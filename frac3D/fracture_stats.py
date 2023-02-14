import numpy as np
from scipy.spatial import Delaunay
import mplstereonet as mpl
import math
import open3d as o3d
import concurrent.futures
from .geometry import Cube, surface_area, triangle_area, plane_rotate, intersect_line_triangle


def fisher_stats(plunge, trend):

    trend_rad = trend*np.pi/180 #imprime usando a função e trend em rad
    plunge = plunge*np.pi/180 #plung em rad

    #########CRONIN2008
    l = np.cos(plunge) * np.cos(trend_rad)
    m = np.cos(plunge) * np.sin(trend_rad)
    n = np.sin(plunge)
        
    #mean dip vector l,m,n   (soma dos valores)
    l_soma = l.sum()
    m_soma = m.sum()
    n_soma = n.sum()

    #R = vetor resultante
    r = np.sqrt(l_soma**2 + m_soma**2 + n_soma**2)

    #unit vector l,m,n
    l_unit = l_soma/r
    m_unit = m_soma/r
    n_unit = n_soma/r

    #delta (ângulo plunge of mean dip vector) em graus
    delta_rad = np.arcsin(n_unit)
    delta_graus= np.arcsin(n_unit)*180/np.pi

    #trend (direção média do conjunto) em graus
    trend_m = np.where(m_unit < 0, 360- np.arccos(l_unit/np.cos(delta_rad)), np.arccos(l_unit/np.cos(delta_rad)))*180/np.pi

    #k (dispersão dos dados)
    qt = len(trend)
    k = (qt-1)/(qt-r)

    #O raio do cone do intervalo de confiança de 95% (α95) em graus
    # alpha95_rad = np.arccos(1-((qt-r)/r)*(((1/0.05)**(1/(qt-1)))-1))
    # alpha95_graus = alpha95_rad*180/np.pi

    #Erro total: soma o erro aleatório (do equipamento) - para Brunton 0.5° a 1°
    # r_total = np.sqrt((alpha95_rad**2) + ((1*np.pi/180)**2))
    # r_total_graus = r_total*180/np.pi

    return (trend_m, delta_graus), (qt, r, k)

def compare_centers(center_a, center_b):
    
    strike_a, dip_a = center_a
    strike_b, dip_b = center_b

    lat_a, lon_a = mpl.analysis._convert_measurements((strike_a, dip_a), measurement='poles')
    lat_b, lon_b = mpl.analysis._convert_measurements((strike_b, dip_b), measurement='poles')


    angle = mpl.stereonet_math.angular_distance((lat_a, lon_a), (lat_b, lon_b))
    
    return angle

# def normals_to_stikedip(normals):
    
#     strike = []
#     dip = []
    
#     # print(np.shape(normals)[0])
#     for i in range(0, np.shape(normals)[0]):
        
#         normal = normals[i]
#         # centroid = centroids[i]
    
#         # d = -centroid.dot(normal)

#         a, b, c = normal
#         alpha = math.acos(abs(c/(math.sqrt((math.pow(a, 2) + math.pow(b, 2)
#                                             + math.pow(c, 2))))))
#         beta = math.acos(a/(math.sqrt(math.pow(a, 2) + math.pow(b, 2))))
    
#         dip.append(math.degrees(alpha))
#         # a_value.append(a)
#         # b_value.append(b)
#         # c_value.append(c)
#         # d_value.append(d)
    
#         if a > 0 and c < 0: # quadrante 2 #*
#             strike.append(360-math.degrees(beta))
#         elif a > 0 and c > 0: # quadrante 3
#             strike.append(math.degrees(beta)) #* 180+
#         elif a < 0 and c < 0: # quadrante 4
#             strike.append(180+math.degrees(beta)) #*
#         else: # quadrante 1
#             strike.append(180-math.degrees(beta))
            
#     strike = np.asarray(strike)
#     dip = np.asarray(dip)
            
#     return strike, dip


def compute_fracture_sets_spacing(labels, n_clusters, a, b, c, d):
    # Compute spacing statistics. Returns mean distance, standard deviation, and median
    spacing_stats = []
    for cluster in range(0, n_clusters):
        #np.asarray(a)[np.where(labels==cluster)]
        a_mean = np.mean(np.asarray(a)[np.where(labels==cluster)])
        b_mean = np.mean(np.asarray(b)[np.where(labels==cluster)])
        c_mean = np.mean(np.asarray(c)[np.where(labels==cluster)])
    
        d_temp = np.copy(np.asarray(d)[np.where(labels==cluster)])
    
        d_temp.sort()
    
        spacing = []
    
        for i in range(np.size(d_temp)-1):
            SP = abs(d_temp[i] -d_temp[i+1])/math.sqrt(math.pow(a_mean,2) + math.pow(b_mean,2) + math.pow(c_mean,2))
            spacing.append(SP)
    
        if np.size(spacing) > 1:
            spacing_stats.append((np.mean(spacing), np.std(spacing), np.median(spacing)))
        else:
            print(spacing)
    spacing_stats = np.reshape(spacing_stats, (np.shape(spacing_stats)[0], 3))
    
    return spacing_stats 
    
    
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def compute_p32_statistics(new_normals, new_centroids, model_dimension, cube_size, plane_size, extrapolate):

    # Format inputs to     
    normals = np.asarray([item for sublist in new_normals for item in sublist])
    centroids = np.asarray([item for sublist in new_centroids for item in sublist])
    
  
    # Data range
    x_range = np.asarray([0, model_dimension[0]]) + np.asarray([-extrapolate, extrapolate])
    y_range = np.asarray([0, model_dimension[1]]) + np.asarray([-extrapolate, extrapolate])
    z_range = np.asarray([0, model_dimension[2]]) + np.asarray([-extrapolate, extrapolate])
    
    # Maximum quantity of cubes in each dimension
    maxX = (x_range[1] - x_range[0])/(cube_size*2)
    maxY = (y_range[1] - y_range[0])/(cube_size*2)
    maxZ = (z_range[1] - z_range[0])/(cube_size*2)
    
    offsets = []
    for i in range(0, math.ceil(maxX)):
        for j in range(0, math.ceil(maxY)):
            for k in range(0, math.ceil(maxZ)):
                x0 = float(i*cube_size*2) + cube_size + int(x_range[0])
                y0 = float(j*cube_size*2) + cube_size + int(y_range[0])
                z0 = float(k*cube_size*2) + cube_size + int(z_range[0])
                offsets.append((x0, y0, z0))
    
    
    
    progress_count = 0
    progress_total = np.shape(offsets)[0]*np.shape(normals)[0]
    printProgressBar(0, progress_total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    

    intensity_values = []
    for offset in offsets:
        cube_definitions = Cube(cube_size, offset)
        edges = cube_definitions.edges
        vertices = cube_definitions.vertices
        sum_area = 0


        for p in range(0, np.shape(normals)[0]):            
            frac_points = np.asarray([[-1, -1,0],[1,-1,0],[1,1,0], [-1,1,0]]) # base horizontal plane
            frac_points = frac_points*plane_size
            frac_points = plane_rotate(frac_points, normals[p], centroids[p])
            
            points = []
   
            # Cube faces
            faces = [[0,2,3], [2,6,3], [0,5,1], [0,2,5], [1,4,5], [4,7,5],
                      [4,3,6], [4,7,6], [0,1,4], [0,3,4], [2,5,7], [2,6,7]]
            np_triangles = np.array(faces).astype(np.int32)
            cube_mesh = o3d.geometry.TriangleMesh()
            cube_mesh.vertices = o3d.utility.Vector3dVector(vertices)
            cube_mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
            
            
            face_mesh = o3d.geometry.TriangleMesh()
            np_vertices = np.array(frac_points)
            np_triangles = np.array([[0, 1, 2], [0, 2, 3]]).astype(np.int32)
            
            face_mesh.vertices = o3d.utility.Vector3dVector(np_vertices)
            face_mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
            face_mesh.compute_vertex_normals()
            
            # o3d.visualization.draw_geometries([face_mesh, cube_mesh], mesh_show_back_face=True, mesh_show_wireframe=True)
            
            
            if face_mesh.is_intersecting(cube_mesh):
                
                
                # check if points from fracture plane are inside cube
                for i in range(0, 4):
                    # res = in_poly_hull_single(vertices, frac_points[i])
                    res = Delaunay(vertices).find_simplex(frac_points[i])
                    if res >=0:
                        points.append(frac_points[i])

                
                # check if fracture edges cross the cube faces
                for i in range(0, 4): 
                    if i == 3:
                        point_a = frac_points[i]
                        point_b = frac_points[0]
                    else:
                        point_a = frac_points[i]
                        point_b = frac_points[i+1]
                
                
                    for j in range(0, 6): # for each cube face
                        face = cube_definitions.faces[j]
                        point = intersect_line_triangle(point_a, point_b, 
                                                        face[0],face[1],face[2])
                        
                        if isinstance(point, np.ndarray):
                            points.append(point)
                        
                        point = intersect_line_triangle(point_a, point_b, 
                                                        face[0],face[2],face[3])
                        
                        if isinstance(point, np.ndarray):
                            points.append(point)

                
                # check cube edges cross the fracture area
                for edge in edges:
                    point = intersect_line_triangle(edge[0], edge[1],
                                                    frac_points[0], frac_points[1], frac_points[2])
                    
                    if isinstance(point, np.ndarray):
                        points.append(point)
                    
                    point = intersect_line_triangle(edge[0], edge[1],
                                                    frac_points[0], frac_points[2], frac_points[3])
                    
                    if isinstance(point, np.ndarray):
                        points.append(point)


            if np.shape(points)[0] != 0:
                points = np.asarray(points)
                if np.shape(points)[0] > 3:
                    sum_area += surface_area(points)
                else:
                    sum_area += triangle_area(points)

            progress_count+=1
            printProgressBar(progress_count, progress_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

        P32 = sum_area/((cube_size*2)*(cube_size*2)*(cube_size*2))
        intensity_values.append(P32)

    return intensity_values, offsets



def compute_p32_statistics_thread(new_normals, new_centroids, model_dimension, cube_size=4, plane_size=4, extrapolate=0):

    # Format inputs to     
    normals = np.asarray([item for sublist in new_normals for item in sublist])
    centroids = np.asarray([item for sublist in new_centroids for item in sublist])
    
  
    # Data range
    x_range = np.asarray([0, model_dimension[0]]) + np.asarray([-extrapolate, extrapolate])
    y_range = np.asarray([0, model_dimension[1]]) + np.asarray([-extrapolate, extrapolate])
    z_range = np.asarray([0, model_dimension[2]]) + np.asarray([-extrapolate, extrapolate])
    
    # Maximum quantity of cubes in each dimension
    maxX = (x_range[1] - x_range[0])/(cube_size*2)
    maxY = (y_range[1] - y_range[0])/(cube_size*2)
    maxZ = (z_range[1] - z_range[0])/(cube_size*2)
    
    offsets = []
    for i in range(0, math.ceil(maxX)):
        for j in range(0, math.ceil(maxY)):
            for k in range(0, math.ceil(maxZ)):
                x0 = float(i*cube_size*2) + cube_size + int(x_range[0])
                y0 = float(j*cube_size*2) + cube_size + int(y_range[0])
                z0 = float(k*cube_size*2) + cube_size + int(z_range[0])
                offsets.append((x0, y0, z0))
    
    
    
    # progress_count = 0
    # progress_total = np.shape(offsets)[0]*np.shape(normals)[0]
    # printProgressBar(0, progress_total, prefix = 'Progress:', suffix = 'Complete', length = 50)
    

    intensity_values = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        future = executor.submit(compute_box_intensity, offsets, normals, centroids, plane_size, cube_size)
        intensity_values.append(future.result())
        # print(future.result())
        
        
    # intensity_values = []
    # thread_pool = []
    # with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
    #     for i in range(0,10):
    #         future = executor.submit(compute_box_intensity, offsets[np.int64(np.shape(offsets)[0]/10*i):np.int64((i+1)*np.shape(offsets)[0]/10)], normals, centroids, plane_size, cube_size)
    #         thread_pool.append(future)
            
    #     for thread in thread_pool:
    #         intensity_values.append(thread.result())
    


    return intensity_values[0], np.asarray(offsets)


def compute_box_intensity(offsets, normals, centroids, plane_size, cube_size):
    
    intensity_values = []
    for offset in offsets:
        cube_definitions = Cube(cube_size, offset)
        edges = cube_definitions.edges
        vertices = cube_definitions.vertices
        sum_area = 0


        for p in range(0, np.shape(normals)[0]):            
            frac_points = np.asarray([[-1, -1,0],[1,-1,0],[1,1,0], [-1,1,0]]) # base horizontal plane
            frac_points = frac_points*plane_size
            frac_points = plane_rotate(frac_points, normals[p], centroids[p])
            
            points = []
   
            # Cube faces
            faces = [[0,2,3], [2,6,3], [0,5,1], [0,2,5], [1,4,5], [4,7,5],
                      [4,3,6], [4,7,6], [0,1,4], [0,3,4], [2,5,7], [2,6,7]]
            np_triangles = np.array(faces).astype(np.int32)
            cube_mesh = o3d.geometry.TriangleMesh()
            cube_mesh.vertices = o3d.utility.Vector3dVector(vertices)
            cube_mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
            
            
            face_mesh = o3d.geometry.TriangleMesh()
            np_vertices = np.array(frac_points)
            np_triangles = np.array([[0, 1, 2], [0, 2, 3]]).astype(np.int32)
            
            face_mesh.vertices = o3d.utility.Vector3dVector(np_vertices)
            face_mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
            face_mesh.compute_vertex_normals()
            
            # o3d.visualization.draw_geometries([face_mesh, cube_mesh], mesh_show_back_face=True, mesh_show_wireframe=True)
            
            
            if face_mesh.is_intersecting(cube_mesh):
                
                
                # check if points from fracture plane are inside cube
                for i in range(0, 4):
                    # res = in_poly_hull_single(vertices, frac_points[i])
                    res = Delaunay(vertices).find_simplex(frac_points[i])
                    if res >=0:
                        points.append(frac_points[i])

                
                # check if fracture edges cross the cube faces
                for i in range(0, 4): 
                    if i == 3:
                        point_a = frac_points[i]
                        point_b = frac_points[0]
                    else:
                        point_a = frac_points[i]
                        point_b = frac_points[i+1]
                
                
                    for j in range(0, 6): # for each cube face
                        face = cube_definitions.faces[j]
                        point = intersect_line_triangle(point_a, point_b, 
                                                        face[0],face[1],face[2])
                        
                        if isinstance(point, np.ndarray):
                            points.append(point)
                        
                        point = intersect_line_triangle(point_a, point_b, 
                                                        face[0],face[2],face[3])
                        
                        if isinstance(point, np.ndarray):
                            points.append(point)

                
                # check cube edges cross the fracture area
                for edge in edges:
                    point = intersect_line_triangle(edge[0], edge[1],
                                                    frac_points[0], frac_points[1], frac_points[2])
                    
                    if isinstance(point, np.ndarray):
                        points.append(point)
                    
                    point = intersect_line_triangle(edge[0], edge[1],
                                                    frac_points[0], frac_points[2], frac_points[3])
                    
                    if isinstance(point, np.ndarray):
                        points.append(point)


            if np.shape(points)[0] != 0:
                points = np.asarray(points)
                if np.shape(points)[0] > 3:
                    sum_area += surface_area(points)
                else:
                    sum_area += triangle_area(points)

            # progress_count+=1
            # printProgressBar(progress_count, progress_total, prefix = 'Progress:', suffix = 'Complete', length = 50)

        P32 = sum_area/((cube_size*2)*(cube_size*2)*(cube_size*2))
        intensity_values.append(P32)
        
    return intensity_values
