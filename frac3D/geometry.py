import numpy as np
import math
from scipy.spatial import ConvexHull, Delaunay
import numba
from sklearn.decomposition import PCA


def quadrante(valor):
    a = np.where(valor > 270, valor,(np.where(valor > 180 , valor -180, (np.where(valor>90,valor+180,valor)))))
    return a

def plane_intersect(a, b):
    """
    a, b   4-tuples/lists
           Ax + By +Cz + D = 0
           A,B,C,D in order  

    output: 2 points on line of intersection, np.arrays, shape (3,)
    """
    a_vec, b_vec = np.array(a[:3]), np.array(b[:3])

    aXb_vec = np.cross(a_vec, b_vec)

    A = np.array([a_vec, b_vec, aXb_vec])
    d = np.array([-a[3], -b[3], 0.]).reshape(3,1)

    p_inter = np.linalg.solve(A, d).T

    return p_inter[0], (p_inter + aXb_vec)[0]


def vnorm(vector):
    vector = np.asarray(vector, dtype=np.float64)
    return math.sqrt(np.sum(np.power(vector, 2)))

def plane_rotate(points, normal, centroid):
    normal0 = np.cross(points[2], points[1])
    costheta = np.dot(normal0, np.reshape(normal, (-1,1)))/(vnorm(normal0) *vnorm(normal))
    axis = np.cross(normal0, normal)/vnorm(np.cross(normal0, normal))

    _c = costheta
    s = np.sqrt(1 - _c * _c)
    C = 1 - _c

    x, y, z = axis
    rmat = np.asarray([[x*x*C+_c, x*y*C-z*s, x*z*C+y*s],
                      [y*x*C+z*s, y*y*C+_c, y*z*C-x*s],
                      [z*x*C-y*s, z*y*C+x*s, z*z*C+_c]])

    rmat = np.reshape(rmat, (3,3))

    new_points = []
    for point in points:
        new_points.append(np.dot(rmat, point) + centroid)

    new_points = np.asarray(new_points)

    return new_points

def point_distance(a, b):
    x0, y0, z0 = a
    x1, y1, z1 = b
    d = math.sqrt(math.pow((x0 - x1), 2) + math.pow((y0 - y1), 2) + math.pow((z0 - z1), 2))
    return d

def triangle_area(points):
    distA = point_distance(points[0], points[1])
    distB = point_distance(points[1], points[2])
    distC = point_distance(points[2], points[0])

    p = (distA + distB + distC)/2

    A = math.sqrt( p*(p-distA)*(p-distB)*(p-distC))

    return A

def surface_area(points):
    tri = Delaunay(points[:,0:2], qhull_options='QJ')
    
    sum_area = 0
    for triangle in range(0, np.shape(tri.simplices)[0]):
        sum_area += triangle_area(points[tri.simplices[triangle]])
    return sum_area
	
	
class Cube():
    def __init__(self, cube_size, offset):
        # vertices do cubo
        V1 = [-1,1,1]
        V2 = [1,1,1]
        V3 = [-1,1,-1]
        V4 = [-1,-1,1]
        V5 = [1,-1,1]
        V6 = [1,1,-1]
        V7 = [-1,-1,-1]
        V8 = [1,-1,-1]

        # Face do cubo
        F1 = [V1, V2, V6, V3]
        F2 = [V2, V5, V8, V6]
        F3 = [V5, V4, V7, V8]
        F4 = [V4, V1, V3, V7]
        F5 = [V5, V2, V1, V4]
        F6 = [V6, V3, V7, V8]

        # Arestas do cubo
        E1 = [V1, V2]
        E2 = [V1, V3]
        E3 = [V2, V6]
        E4 = [V2, V5]
        E5 = [V1, V4]
        E6 = [V4, V5]
        E7 = [V4, V7]
        E8 = [V7, V8]
        E9 = [V7, V3]
        E10 = [V3, V6]
        E11 = [V6, V8]
        E12 = [V8, V5]

        edges = np.asarray([E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12])
        edges = edges*cube_size
        self.edges = edges + offset

        faces = np.asarray([F1, F2, F3, F4, F5, F6])
        faces = faces*cube_size
        self.faces = faces + offset

        vertices = np.asarray([V1, V2, V3, V4, V5, V6, V7, V8])
        vertices = vertices*cube_size
        self.vertices = vertices + offset
        
        
def in_poly_hull_single(poly, point):
    hull = ConvexHull(poly, qhull_options='QJ')
    new_hull = ConvexHull(np.concatenate((poly, [point])), qhull_options='QJ')
    return np.array_equal(new_hull.vertices, hull.vertices)


def in_poly_hull_multi(poly, points):
    hull = ConvexHull(poly)
    res = []
    for p in points:
        new_hull = ConvexHull(np.concatenate((poly, [p])))
        res.append(np.array_equal(new_hull.vertices, hull.vertices))
    return res

def compute_normal(points):
    pca = PCA(n_components=3)
    pca.fit(points)
    eig_vec = pca.components_

    normal = eig_vec[2, :]
    return normal

def point_in_segment(termination_A, termination_B, point):
    x1, y1, z1 = termination_A
    x2, y2, z2 = termination_B
    x, y, z = point
    
    AB = np.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
    AP = np.sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))
    PB = np.sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y)+(z2-z)*(z2-z))
    if AB == (AP + PB):
        return True
    else:
        return False

@numba.jit(nopython=True)
def intersect_line_triangle(q1,q2,p1,p2,p3):
    # https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
    def signed_tetra_volume(a,b,c,d):
        return np.sign(np.dot(np.cross(b-a,c-a),d-a)/6.0)

    s1 = signed_tetra_volume(q1,p1,p2,p3)
    s2 = signed_tetra_volume(q2,p1,p2,p3)

    if s1 != s2:
        s3 = signed_tetra_volume(q1,q2,p1,p2)
        s4 = signed_tetra_volume(q1,q2,p2,p3)
        s5 = signed_tetra_volume(q1,q2,p3,p1)
        if s3 == s4 and s4 == s5:
            n = np.cross(p2-p1,p3-p1)
            t = np.dot(p1-q1,n) / np.dot(q2-q1,n)
            return q1 + t * (q2-q1)
    return None