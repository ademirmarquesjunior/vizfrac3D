# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 13:48:12 2022

@author: adeju
"""


import numpy as np
import matplotlib.pyplot as plt
import mplstereonet as mpl
from time import time
from frac3D.file_handling import load_planes_from_file, load_planes_from_file_xpp
from frac3D.plot import compute_fracture_sets_fisher, dfn_plot, fracture_planes_plot, plot_stereonet, compute_fisher_from_normals
from frac3D.fracture_clustering import compute_kmeans_runs, optimal_number_of_clusters, get_fracture_sets
from frac3D.fracture_generator import generator
from frac3D.fracture_stats import  compute_fracture_sets_spacing, compute_p32_statistics, compute_p32_statistics_thread, compare_centers


# Loading Plane data from Mosis XP
# strike, dip, normals, centroids, a, b, c, d = load_planes_from_file('data/SoledadeRavinaVR.data') # Soledade_VR.data
strike, dip, normals, centroids, a, b, c, d = load_planes_from_file_xpp('data/SoledadeRavinaLocal.xpp') 
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
model_dimension = [100, 50, 5]
model_center = [np.max(centroids[:, 0]), np.max(centroids[:, 0]), np.max(centroids[:, 0])]
fracture_plane_size = 4
dfn_cell_size = 2
start = time()
p_32_statistics, offsets = compute_p32_statistics(np.reshape(normals, (1,-1, 3)), np.reshape(centroids-model_center, (1,-1, 3)), model_dimension, dfn_cell_size, fracture_plane_size, 0)
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
kmeans_data2 = np.power(np.nanmean(elbow, axis=1), 2)/max_K2
n_clusters = optimal_number_of_clusters(kmeans_data2[0:10], 0)

# Plot Elbow method
fig = plt.figure(10, figsize=(12,4))
ax1 = fig.add_subplot()
ax1.plot(list(range(1, 11)), kmeans_data2[0:10] , 'go-')
# ax1.set_yscale('log')
# ax1.set_title('K-means Elbow method with Fisher k sum'),
ax1.set_xlabel('Number of clusters', fontsize=14)
ax1.set_ylabel('Normalized squared \n mean of Fisher k', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig.show()

# Number of clusters defined
n_clusters = 7

# Evaluate k-means clustering with and without the stabilization check
k_stats_balanced = []
k_stats_not_balanced = []
for i in range(0, 10):
    labels_temp, centers_temp = get_fracture_sets(strike, dip, n_clusters, balance=False)
    plane_stats_temp = compute_fracture_sets_fisher(strike, dip, plunge, bearing, n_clusters, labels_temp, plot=False)
    k_stats_not_balanced.append(plane_stats_temp[:,4])
    
    labels_temp, centers_temp = get_fracture_sets(strike, dip, n_clusters, balance=True)
    plane_stats_temp = compute_fracture_sets_fisher(strike, dip, plunge, bearing, n_clusters, labels_temp, plot=False)
    k_stats_balanced.append(plane_stats_temp[:,4])
    
k_stats_balanced = np.asarray(k_stats_balanced)
k_stats_not_balanced = np.asarray(k_stats_not_balanced)


# Generate Fisher and spacing statistics
labels, centers = get_fracture_sets(strike, dip, n_clusters, balance=True)
# np.save('cluster_labels.npy', labels)
# labels = np.load('cluster_labels.npy')

plane_stats = compute_fracture_sets_fisher(strike, dip, plunge, bearing, n_clusters, labels, plot=True)
spacing_stats = compute_fracture_sets_spacing(labels, n_clusters, a, b, c, d)

plot_stereonet(strike, dip, labels, n_clusters, plane_stats, mode='pole')
# plot_stereonet(strike, dip, labels, n_clusters, plane_stats, mode='plane')



# Generate new DFN model
model_dimension = [40, 40, 40] # 40. 40, 40
fracture_plane_size = 4 
dfn_cell_size = 0.5 # 0.5

# Gerar novos planos de fraturas e computar estatísticas
new_normals, new_centroids = generator(n_clusters, plane_stats, spacing_stats, model_dimension, fracture_plane_size)

start = time()
# p_32_statistics, offsets = compute_p32_statistics(new_normals, new_centroids, model_dimension, dfn_cell_size, fracture_plane_size, 0)
p_32_statistics, offsets = compute_p32_statistics_thread(new_normals, new_centroids, model_dimension, dfn_cell_size, fracture_plane_size, 0)
end = time()
print(end-start)


# Plot fracture intensity model and generated fractures
fracture_planes_plot(n_clusters, new_normals, new_centroids, model_dimension, fracture_plane_size)
dfn_plot(p_32_statistics, offsets, dfn_cell_size, -1, 200, model_dimension, only_walls=True)

# Plot statistics
mean_p_32 = np.mean(p_32_statistics)
sd_mean_p32 = np.std(p_32_statistics)
print(mean_p_32, sd_mean_p32)


plane_stats2 = compute_fisher_from_normals(new_normals)
plot_stereonet_from_normals(new_normals, mode='pole')



# Validation of the stochastic generation
model_dimension = [40, 40, 40] # 40, 40, 40
fracture_plane_size = 4 


distances = []
directions = []
for i in range(0, 30):
    new_normals_temp, new_centroids_temp = generator(n_clusters, plane_stats, spacing_stats, model_dimension, fracture_plane_size)
    plane_stats_temp_stochastic = compute_fisher_from_normals(new_normals_temp, plot=False)
    for i in range(0, n_clusters):
        distances.append(compare_centers(plane_stats[i, 0:2], plane_stats_temp_stochastic[i, 0:2])[0])
        directions.append(plane_stats_temp_stochastic[i,0:2])
distances = np.reshape(distances, (30, 7))
directions = np.reshape(directions, (30, 7, 2))

np.mean(distances, axis=0)
np.std(distances, axis=0)

for i in range(0, n_clusters):
    plunge_temp, bearing_temp = mpl.pole2plunge_bearing(directions[:, i, 0], directions[:, i, 1])
    vector, stats = mpl.find_fisher_stats(plunge_temp, bearing_temp, conf=95)
    mean_strike, mean_dip = mpl.plunge_bearing2pole(vector[0], vector[1])
    print(mean_strike, mean_dip)
    


# # Save model
# with open('model_data_40model_1_cell.npy', 'wb') as f:
#     np.save(f, new_normals)
#     np.save(f, new_centroids)
#     np.save(f, model_dimension)
#     np.save(f, fracture_plane_size)
#     np.save(f, dfn_cell_size)
#     np.save(f, p_32_statistics)
#     np.save(f, offsets)

# Load model
# with open('model_data_40model_1_cell.npy', 'rb') as f:
#     new_normals = np.load(f, allow_pickle=True)
#     new_centroids = np.load(f, allow_pickle=True)
#     model_dimension = np.load(f, allow_pickle=True)
#     fracture_plane_size = np.load(f, allow_pickle=True)
#     dfn_cell_size = np.load(f, allow_pickle=True)
#     p_32_statistics = np.load(f, allow_pickle=True)
#     offsets = np.load(f, allow_pickle=True)
