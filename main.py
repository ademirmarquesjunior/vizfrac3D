# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 13:48:12 2022

@author: adeju
"""


import numpy as np
import matplotlib.pyplot as plt
import mplstereonet as mpl
from time import time
from frac3D.file_handling import load_planes_from_file
from frac3D.plot import *
from frac3D.fracture_clustering import compute_kmeans_runs, optimal_number_of_clusters, get_fracture_sets
from frac3D.fracture_generator import generator
from frac3D.fracture_stats import  compute_fracture_sets_spacing # compute_fracture_sets_fisher,


# Loading Plane data from Mosis XP
strike, dip, normals, centroids, a, b, c, d = load_planes_from_file('data/SoledadeRavinaVR.data') # Soledade_VR.data
centroids = centroids-[np.min(centroids[:,0]), np.min(centroids[:,1]), np.min(centroids[:,2])]

strike = np.where(strike > 180, strike - 180, strike)
plunge, bearing = mpl.pole2plunge_bearing(strike, dip)
plunge = np.asarray(plunge)
bearing = np.asarray(bearing)
rake = mpl.project_onto_plane(strike, dip, plunge, bearing)
rake = np.where(np.isnan(rake), 89, rake)


# Plot deterministic fracture planes
fracture_planes_plot(n_clusters=1, new_normals=np.reshape(normals, (1,-1, 3)), new_centroids=np.reshape(centroids, (1,-1, 3)), model_dimension = [100,100,100], plane_size=4)


# Compute fracture intensity and plot deterministics DFN 
model_dimension = [100, 60, 5]
fracture_plane_size = 4
dfn_cell_size = 2
start = time()
p_32_statistics, offsets = compute_p32_statistics(np.reshape(normals, (1,-1, 3)), np.reshape(centroids, (1,-1, 3)), model_dimension, dfn_cell_size, fracture_plane_size, 0)
end = time()
print(end-start)
dfn_plot(p_32_statistics, offsets, dfn_cell_size, -1, 200, model_dimension, only_walls=False)


# Compute k-means runs for the modified elbow method
max_clusters=10
max_iter=50
start = time()
kmeans_data, elbow = compute_kmeans_runs(strike, dip, plunge, bearing, max_clusters, max_iter)
# np.save('kmeans_data.npy', kmeans_data)
# np.save('elbow_runs.npy', elbow)
end = time()
print(end-start)


# Compute the ideal number of clusters 
# kmeans_data = np.load('kmeans_data.npy')
# elbow = np.load('elbow_runs.npy')
elbow = np.reshape(elbow, (int(np.shape(elbow)[0]/max_clusters), max_clusters))
elbow = np.where(elbow == np.inf, np.nan, elbow) 
max_K2 = np.nanmax(np.power(np.nanmean(elbow, axis=1), 2))
kmeans_data2 = max_K2-np.power(np.nanmean(elbow, axis=1), 2)
n_clusters = optimal_number_of_clusters(kmeans_data, 1)


# Plot Elbow method
fig = plt.figure(10, figsize=(8,4))
ax1 = fig.add_subplot()
ax1.plot(list(range(1, 10)), kmeans_data[0:9] , 'go-')
# ax1.set_yscale('log')
# ax1.set_title('K-means Elbow method with Fisher k sum'),
ax1.set_xlabel('Number of clusters', fontsize=14)
ax1.set_ylabel('Sum of Fisher k', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig.show()


# Generate Fisher and spacing statistics
n_clusters = 7
labels, centers = get_fracture_sets(strike, dip, n_clusters)
# np.save('cluster_labels.npy', labels)
# labels = np.load('cluster_labels.npy')

plane_stats = compute_fracture_sets_fisher(strike, dip, plunge, bearing, n_clusters, labels)
spacing_stats = compute_fracture_sets_spacing(labels, n_clusters, a, b, c, d)

plot_stereonet(strike, dip, labels, n_clusters, plane_stats, mode='pole')
# plot_stereonet(strike, dip, labels, n_clusters, plane_stats, mode='plane')



# Generate new DFN model
model_dimension = [40, 40, 40] # 40. 40, 40
fracture_plane_size = 4 
dfn_cell_size = 0.5

# Gerar novos planos de fraturas e computar estatísticas
new_normals, new_centroids = generator(n_clusters, plane_stats, spacing_stats, model_dimension, fracture_plane_size)

start = time()
p_32_statistics, offsets = compute_p32_statistics(new_normals, new_centroids, model_dimension, dfn_cell_size, fracture_plane_size, 0)
end = time()
print(end-start)


# Plot fracture intensity model and generated fractures
fracture_planes_plot(n_clusters, new_normals2, new_centroids, model_dimension, fracture_plane_size)
dfn_plot(p_32_statistics, offsets, dfn_cell_size, 0, 200, model_dimension, only_walls=True)

# Plot statistics
mean_p_32 = np.mean(p_32_statistics)
sd_mean_p32 = np.std(p_32_statistics)
print(mean_p_32, sd_mean_p32)


# # Save model
# with open('model_data_40model_1_cell.npy', 'wb') as f:
#     np.save(f, new_normals)
#     np.save(f, new_centroids)
#     np.save(f, model_dimension)
#     np.save(f, fracture_plane_size)
#     np.save(f, dfn_cell_size)
#     np.save(f, p_32_statistics)
#     np.save(f, offsets)

# # Load model
# with open('model_data_40model_1_cell.npy', 'rb') as f:
#     new_normals = np.load(f, allow_pickle=True)
#     new_centroids = np.load(f, allow_pickle=True)
#     model_dimension = np.load(f, allow_pickle=True)
#     fracture_plane_size = np.load(f, allow_pickle=True)
#     dfn_cell_size = np.load(f, allow_pickle=True)
#     p_32_statistics = np.load(f, allow_pickle=True)
#     offsets = np.load(f, allow_pickle=True)