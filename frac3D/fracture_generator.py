import numpy as np
import mplstereonet as mpl
import spherical_kde


def generator(n_clusters, direction_stats, spacing_stats, model_dimension, plane_size=10):
    
    max_x, max_y, max_z = model_dimension
    
    new_normals = []
    # new_distances = []
    new_centroids = []
    
    n_levels = int(max_z/plane_size)
    if n_levels == 0: n_levels = 1
    
    
    for group in range(0, n_clusters):
        normals_temp = []
        centroids_temp = []
        
        for i in range(0, n_levels):            
            # Compute plane quantities
            n_fractures = int((max_x*max_y)/spacing_stats[group,0])
            if n_fractures == 0:
                n_fractures = 1

            # Compute plane directions
            mean_strike, mean_dip, vector, r_fisher, k_fisher = direction_stats[group]
            lon, lat = mpl.pole(mean_strike, mean_dip)
            new_lat, new_lon = spherical_kde.distributions.VonMisesFisher_sample(lat[0], lon[0], np.deg2rad(r_fisher), size=n_fractures)
            #newstrike, new_dip = mpl.geographic2pole(new_lon, new_lat)
            normal = mpl.stereonet2xyz(new_lon, new_lat)
            normal = [normal[0], normal[1], normal[2]]
            
            normals_temp.append(np.rot90(normal))

            centroids = np.rot90((np.random.uniform(0, max_x, size=n_fractures), np.random.uniform(0, max_y, size=n_fractures), np.full(n_fractures, i*plane_size) + np.random.random(n_fractures)*plane_size))
            centroids_temp.append(centroids)
            

        
        normals_temp = np.asarray(normals_temp)
        normals_temp = np.reshape(normals_temp, (int(np.size(normals_temp)/3), 3))
        new_normals.append(normals_temp)
        
        centroids_temp = np.asarray(centroids_temp)
        centroids_temp = np.reshape(centroids_temp, (int(np.size(centroids_temp)/3), 3))
        new_centroids.append(centroids_temp)
        
    
    return new_normals, new_centroids
