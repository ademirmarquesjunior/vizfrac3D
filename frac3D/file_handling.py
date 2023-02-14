import numpy as np
# import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
# from sklearn.preprocessing import normalize

import math
#from sympy import Point3D, Line3D, Plane
# from scipy.spatial import Delaunay
import csv
import shapefile
#import matplotlib.colors as mcolors
# %matplotlib inline

#!pip install mplstereonet
import mplstereonet as mpl
from .plot import normals_to_stikedip

# from plotly.offline import plot
# import plotly.graph_objects as go
# import plotly.express as px
# import pandas as pd

# #!pip install spherical_kde
# import spherical_kde
# import pycircstat as circ



def load_planes_from_file(File):
    #File = 'Soledade_VR.data'
    strike = []
    dip = []
    centroids = []
    normals = []
    a_value = []
    b_value = []
    c_value = []
    d_value = []
    with open(File) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')

        row_count = 0
        for row in csv_reader:
          for cell in row:
            if cell == 'Plane':
              planes_start = int(row[1])
              planes_end = int(row[2]) + int(row[1])
              break

          if row_count >= planes_start and row_count < planes_end:
            points = []
            for i in range(np.size(row)):
              if i > 2 and i < np.size(row)-1:
                points.append(float(row[i]))
            points = np.reshape(points, (int(np.size(points)/3), 3))
            points = np.asarray([points[:,0], points[:,2], -points[:,1]]).T

            pca = PCA(n_components=3)
            pca.fit(points)
            eig_vec = pca.components_

            normal = eig_vec[2, :]  # (a, b, c)
            normals.append(normal)
            centroid = np.mean(points, axis=0)
            centroids.append(centroid)

            d = -centroid.dot(normal)

            a, b, c = normal
            alpha = math.acos(abs(c/(math.sqrt((math.pow(a, 2) + math.pow(b, 2)
                                                + math.pow(c, 2))))))
            beta = math.acos(a/(math.sqrt(math.pow(a, 2) + math.pow(b, 2))))

            dip.append(math.degrees(alpha))
            a_value.append(a)
            b_value.append(b)
            c_value.append(c)
            d_value.append(d)

            if a > 0 and b > 0: # quadrante 2 #*
                strike.append(360-math.degrees(beta))
            elif a > 0 and b < 0: # quadrante 3
                strike.append(180+math.degrees(beta)) #* 180+
            elif a < 0 and b < 0: # quadrante 4
                strike.append(math.degrees(beta)) #*
            else: # quadrante 1
                strike.append(math.degrees(beta))

          row_count += 1
          
    strike = np.asarray(strike)
    dip = np.asarray(dip)
    centroids = np.asarray(centroids)
    normals = np.asarray(normals)

    return strike, dip, normals, centroids, a_value, b_value, c_value, d_value


def load_planes_from_file_xpp(File):
    #File = 'SoledadeRavinaLocal.xpp'
    strike = []
    dip = []
    centroids = []
    normals = []
    a_value = []
    b_value = []
    c_value = []
    d_value = []
    with open(File) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='|')

        row_count = 0
        for row in csv_reader:
          for cell in row:
            if cell == 'Plane':
              planes_start = int(row[2])
              planes_end = int(row[2]) + int(row[1])
              break
          if 'planes_start' in locals():
              if row_count >= planes_start and row_count < planes_end:
                points = []
                row_values = row[6].split()
                for i in range(np.size(row_values)):
                  if i > 2 and i < np.size(row_values):
                    points.append(float(row_values[i]))
                points = np.reshape(points, (int(np.size(points)/3), 3))
                points = np.asarray([points[:,0], points[:,2], -points[:,1]]).T
    
                pca = PCA(n_components=3)
                pca.fit(points)
                eig_vec = pca.components_
    
                normal = eig_vec[2, :]  # (a, b, c)
                normals.append(normal)
                centroid = np.mean(points, axis=0)
                centroids.append(centroid)
    
                d = -centroid.dot(normal)
    
                a, b, c = normal
                alpha = math.acos(abs(c/(math.sqrt((math.pow(a, 2) + math.pow(b, 2)
                                                    + math.pow(c, 2))))))
                beta = math.acos(a/(math.sqrt(math.pow(a, 2) + math.pow(b, 2))))
    
                dip.append(math.degrees(alpha))
                a_value.append(a)
                b_value.append(b)
                c_value.append(c)
                d_value.append(d)
    
                if a > 0 and b > 0: # quadrante 2 #*
                    strike.append(360-math.degrees(beta))
                elif a > 0 and b < 0: # quadrante 3
                    strike.append(180+math.degrees(beta)) #* 180+
                elif a < 0 and b < 0: # quadrante 4
                    strike.append(math.degrees(beta)) #*
                else: # quadrante 1
                    strike.append(math.degrees(beta))

          row_count += 1
          
    strike = np.asarray(strike)
    dip = np.asarray(dip)
    centroids = np.asarray(centroids)
    normals = np.asarray(normals)

    return strike, dip, normals, centroids, a_value, b_value, c_value, d_value



def load_planes_from_shapefile(file_address):
    # file_address = 'data/GaivotaLocal_planeA'
    sf = shapefile.Reader(file_address)
    
    bbox = sf.bbox
    
    model_center = [np.mean([bbox[0], bbox[2]]), np.mean([bbox[1], bbox[3]])]
    
    
    
    shapes = sf.shapes()   
    # fields = sf.fields
    rec = sf.records()
    
    # data = np.asarray(rec) # Corrigir problemas com decimais
    
    strike = []
    dip = []
    centroids = []
    normals = []
    a_value = []
    b_value = []
    c_value = []
    d_value = []
    
    for i in range(0, len(shapes)):
        centroid = np.asarray((rec[i][2], rec[i][4], rec[i][3]))
        centroids.append(centroid)
        strike.append(rec[i][6])
        dip.append(rec[i][5])
        
        plunge, bearing = mpl.pole2plunge_bearing(rec[i][6], rec[i][5])
        
        trend_rad = bearing*np.pi/180 #imprime usando a função e trend em rad
        plunge = plunge*np.pi/180 #plung em rad

        #########CRONIN2008
        l = np.cos(plunge) * np.cos(trend_rad)
        m = np.cos(plunge) * np.sin(trend_rad)
        n = np.sin(plunge)
            
        #mean dip vector l,m,n   (soma dos valores)
        l_soma, m_soma, n_soma = l.sum(), m.sum(), n.sum()

        #R = vetor resultante
        r = np.sqrt(l_soma**2 + m_soma**2 + n_soma**2)

        #unit vector l,m,n
        normal = np.asarray((l_soma/r, m_soma/r, n_soma/r))
        a, b, c = normal
        d = -centroid.dot(normal)
        
        a_value.append(a)
        b_value.append(b)
        c_value.append(c)
        d_value.append(d)
        normals.append(normal)
        
    strike = np.asarray(strike)
    dip = np.asarray(dip)
    normals = np.asarray(normals)
    centroids = np.asarray(centroids)
 

    return strike, dip, normals, centroids, a_value, b_value, c_value, d_value, model_center




def load_lines_from_shapefile(file_address):
    # file_address = 'data/GaivotaLocal_line'
    
    sf = shapefile.Reader(file_address)
    
    shapes = sf.shapes()   
    # fields = sf.fields
    rec = sf.records()
    
    
    
    lines = []
    for i in range(0, len(shapes)):
        data =  sf.shapeRecords()[i].shape.__geo_interface__
        coordinates = data['coordinates']
        coordinates = np.asarray((coordinates[0], coordinates[-1]), dtype=np.float64)
        lines.append(coordinates[0])
        
    lines = np.asarray(lines)
    
    lines[:, 0] = np.min(lines[:, 0])
    
    lines[:, 1] = np.min(lines[:, 1])
    
    maxX = np.max(lines[:, 0])
    maxY = np.man(lines[:, 1])
    
        
    return lines, maxX. maxY



def save_shapefile(file_address, regression_lines, segm_group_angles, dataset):
    '''
    Save line segments to shapefiles.

    Parameters
    ----------
    file_address : String
        DESCRIPTION.
    regression_lines : Array of int64
        Array of line segments.
    segm_group_angles : Array of int64
        Array of angles and line lengths.
    dataset : Object
        Rasterio object.

    Returns
    -------
    None.

    '''
    w = shapefile.Writer(file_address, shapeType=3)
    w.field('fracture_i', 'N')
    w.field('fractdir', 'N')
    w.field('fractlength', 'N')

    for i in range(0, np.shape(regression_lines)[0]):
        point0 = dataset.xy(regression_lines[i, 1], regression_lines[i, 0])
        point1 = dataset.xy(regression_lines[i, 3], regression_lines[i, 2])
        # Add record
        w.record(i, segm_group_angles[i, 0], 0)
        # Add geometry
        w.line([[[point0[0], point0[1]], [point1[0], point1[1]]]])

    w.close()
    return

def export_dfn_csv(file_address, offsets, p_32_statistics):
    # file_address = 'dfn.csv'
    
    with open(file_address, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ', quoting=csv.QUOTE_MINIMAL)
        
        spamwriter.writerow(['x', 'y', 'z', 'value'])
        

        for i in range(0, np.shape(offsets)[0]):
            position = offsets[i]
            spamwriter.writerow([position[0], position[1], position[2], p_32_statistics[i]])
            
            
            
def export_planes_from_normals(file_address, new_normals, new_centroids, model_center, z_offset=30):
    
    if '_dip' in dir():
        del(_dip)
    
    with open(file_address, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        
        spamwriter.writerow(['x', 'y', 'z', 'strike', 'dip'])
    

        for i in range(0, np.shape(new_normals)[0]):
            positions = np.copy(new_centroids[i])
            positions[:,0] = model_center[0] + positions[:,0]
            positions[:,1] = model_center[1] - positions[:,1]
            positions[:,2] = positions[:,2]

            strike, dip = normals_to_stikedip(new_normals[i]) # mpl.vector2pole(new_normals[group][:,0], new_normals[group][:,1], new_normals[group][:,2])
            strike = np.where(strike > 180, strike - 180, strike)
            plunge, bearing = mpl.pole2plunge_bearing(strike, dip)
            plunge = np.asarray(plunge)
            bearing = np.asarray(bearing)
            
            if '_dip' not in dir():
                _strike = strike
                _dip = dip
                _positions = positions
            else:
                _strike = np.concatenate((_strike, strike), axis=0)
                _dip = np.concatenate((_dip, dip), axis=0)
                _positions = np.concatenate((_positions, positions), axis=0)
                
            
            
    
        for j in range(0, np.shape(_strike)[0]):
            x, y, z = _positions[j]
            z = z + z_offset
            spamwriter.writerow([x, y, z, _strike[j], _dip[j]])