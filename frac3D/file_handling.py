import numpy as np
# import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
# from sklearn.preprocessing import normalize

import math
#from sympy import Point3D, Line3D, Plane
# from scipy.spatial import Delaunay
import csv
#import matplotlib.colors as mcolors
# %matplotlib inline

#!pip install mplstereonet
# import mplstereonet as mpl

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