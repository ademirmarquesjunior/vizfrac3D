import numpy as np
import math
import mplstereonet as mpl


def optimal_number_of_clusters(kmeans_data, start=1):
    x1, y1 = start, kmeans_data[0]
    x2, y2 = len(kmeans_data) + start, kmeans_data[-1]

    distances = []
    for i in range(len(kmeans_data)):
        x0 = i+2
        y0 = kmeans_data[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = math.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        distances.append(numerator/denominator)
    
    return distances.index(min(distances)) + start



def compute_kmeans_runs(strike, dip, plunge, bearing, max_clusters=10, max_iter=50):
    elbow_metric = []

    
    for n_clusters in range(1, max_clusters+1):
      for iter in range(0,max_iter):
          
          for i in range(0,100):
            while True:
                try:
                    centers = mpl.kmeans(strike, dip, num=n_clusters, measurement='poles', tolerance = 0.00005) # Função instável a ser substituída
                except Exception:
                    continue
                break
          #centers = mpl.kmeans(strike, dip, num=n_clusters, measurement='poles', tolerance = 0.00005) # Função instável a ser substituída
          strike_cent, dip_cent = mpl.geographic2pole(*zip(*centers))
    
          labels = np.full((np.size(strike)), 0, dtype=np.int)
          distances = np.full((np.size(strike), 1), np.inf)
    
          for i in range(np.shape(centers)[0]):
              center = centers[i]
              for j in range(0, np.shape(strike)[0]):
                lat, lon = mpl.analysis._convert_measurements((strike[j], dip[j]), measurement='poles')
                angle = mpl.stereonet_math.angular_distance(center, (lat, lon))
                if angle < distances[j]:
                    distances[j] = angle
                    labels[j] = i
    
          plane_stats = []
          for group in range(0, n_clusters):
              if np.size(plunge[np.where(labels == group)]) > 1:
                  vector, stats = mpl.find_fisher_stats(plunge[np.where(labels == group)],
                                                        bearing[np.where(labels == group)], conf=95)
                  mean_strike, mean_dip = mpl.plunge_bearing2pole(vector[0], vector[1])
                  
                  plane_stats.append((mean_strike[0], mean_dip[0], stats[0], stats[1], stats[2]))
                  #print(mpl.find_fisher_stats(strike[np.where(labels == group)], dip[np.where(labels == group)], measurements = 'poles', conf=95))
                  #print(fisher_stats(plunge[np.where(labels == group)], np.asarray(bearing)[np.where(labels == group)]))
              else:
                  plane_stats.append((strike[np.where(labels == group)],
                          dip[np.where(labels == group)], 1, 0, np.inf))
          
          plane_stats = np.reshape(plane_stats, (np.shape(plane_stats)[0], 5))
          if np.shape(plane_stats)[0] > 1:
              #elbow_metric.append(np.rad2deg(circ.mean(np.deg2rad(plane_stats[:,3]))))
              elbow_metric.append(np.sum(plane_stats[:,4]))
          else:
              #elbow_metric.append(plane_stats[0,3])
              elbow_metric.append(plane_stats[0,4])
    
    #plt.plot(elbow_metric, list(range(1,8)))
    elbow_metric = np.reshape(elbow_metric, (int(np.shape(elbow_metric)[0]/max_iter), max_iter))
    elbow_metric = np.where(elbow_metric == np.inf, np.nan, elbow_metric)
    
    max_K2 = np.nanmax(np.power(np.nanmean(elbow_metric, axis=1), 2))
    kmeans_data = max_K2-np.power(np.nanmean(elbow_metric, axis=1), 2)
    
    return kmeans_data, elbow_metric


def get_fracture_sets(strike, dip, n_clusters, balance=True):

    #ax.rake(strike, dip, marker='.', color='black')
    
    count = 0
    cond = True
    
    # Find the two modes
    while cond:
        for i in range(0,100):
          while True:
              try:
                  centers = mpl.kmeans(strike, dip, num=n_clusters, measurement='poles', tolerance = 0.00001) # Função instável a ser substituída
              except Exception:
                  continue
              break
        #centers = mpl.kmeans(strike, dip, num=n_clusters, measurement='poles') # Função instável a ser substituída
        
        
        strike_cent, dip_cent = mpl.geographic2pole(*zip(*centers))
        
        labels = np.full((np.size(strike)), 0, dtype=np.int64)
        distances = np.full((np.size(strike), 1), np.inf)
        
        for i in range(np.shape(centers)[0]):
            center = centers[i]
            for j in range(0, np.shape(strike)[0]):
              lat, lon = mpl.analysis._convert_measurements((strike[j], dip[j]), measurement='poles')
              angle = mpl.stereonet_math.angular_distance(center, (lat, lon))
              if angle < distances[j]:
                  distances[j] = angle
                  labels[j] = i
                  
                  
        # Check if clusters with only one item exists. If true retry       
        for i in range(n_clusters):
            if np.size(strike[np.where(labels == i)]) < 2:
                print("exit condition")
                count += 1
                cond = True
                break
                  

        if balance == True:
            k_fisher = []
            for i in range(n_clusters):
                if np.size(strike[np.where(labels == i)]) > 1:
                    plunge, bearing = mpl.pole2plunge_bearing(strike[np.where(labels == i)], dip[np.where(labels == i)])
                    vector, stats = mpl.find_fisher_stats(plunge, bearing, conf=95)
                    k_fisher.append(stats[2])

        
            if np.size(k_fisher) < n_clusters:
                count += 1
                cond = True
            # if np.sum(np.where(k_fisher <= np.mean(k_fisher)-np.std(k_fisher), k_fisher, 0)) >= 1:
            elif np.sum(np.where(k_fisher < np.max(k_fisher)*0.25, k_fisher, 0)) >= 1:
                count += 1
                cond = True
            else:
                cond = False
                print(str(count+1)+" balancing tries.")
        else:
            cond = False
    
    return labels, centers
