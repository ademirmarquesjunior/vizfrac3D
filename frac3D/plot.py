import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import colorsys
import mplstereonet as mpl
import open3d as o3d
from .geometry import plane_rotate, Cube



def rosechart(angles):

    # angles = segm_group_angles[:, 0]
    bin_edges = np.arange(-5, 366, 10)
    number_of_strikes, bin_edges = np.histogram(angles, bin_edges)

    number_of_strikes[0] += number_of_strikes[-1]

    half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
    two_halves = np.concatenate([half, half])

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='polar')

    ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves,
           width=np.deg2rad(10), bottom=0.0, color='red', edgecolor='k')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
    # ax.set_rgrids(np.arange(1, two_halves.max()+1, 2),angle=0,weight='black')
    ax.set_title("N = "+str(np.size(angles)), y=1.10, fontsize=25)
    # plt.show()
    fig.tight_layout()
    fig.canvas.draw()

    # Now we can save it to a numpy array.
    graph_plot = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    graph_plot = graph_plot.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    plt.show()

    return
	
	
def fisher_plot(plunge, bearing, strike, dip):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='stereonet')
    #ax.line(plunge, bearing, color="black", markersize=2)
    ax.pole(strike, dip, color="black", markersize=2)
    confidence = 95

    # plunge = dip / trend = bearing = strike
    vector, stats = mpl.find_fisher_stats(plunge, bearing, conf=confidence)

    template = (u"Mean Strike/Dip: {strike:0.0f}\u00B0/{dip:0.0f}\u00B0\n"
                "Confidence: {conf}%\n"
                u"Fisher Angle: {fisher:0.2f}\u00B0\n"
                u"R-Value {r:0.3f}\n"
                "K-Value: {k:0.2f}")
    
    mean_strike, mean_dip = mpl.plunge_bearing2pole(vector[0], vector[1])

    #mean_strike = np.where(mean_strike > 180, mean_strike - 180, mean_strike)

    label = template.format(strike=mean_strike[0], dip=mean_dip[0], conf=confidence,
                            r=stats[0], fisher=stats[1], k=stats[2])

    ax.pole(mean_strike, mean_dip, color="red", label=label)
    ax.cone(vector[0], vector[1], stats[1], facecolor="None", edgecolor="red")

    # ax.legend(bbox_to_anchor=(1.1, 1.1), numpoints=1)
    plt.show()
    
    return mean_strike[0], mean_dip[0], stats

    
   
def plot_stereonet(strike, dip, labels, n_clusters, plane_stats, mode='plane'):

    colors = [(12, 12, 255), (24, 150, 33), (255,0,0), (127, 0, 110), (255, 0, 220), (255, 248, 56), (178, 0, 255)]
    colors = np.asarray(colors)/255
    
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111, projection='stereonet')
    if mode == 'pole':
        ax.density_contour(strike, dip, cmap='cool',  measurement='poles', sigma=2)
    
    # for i in range(0, n_clusters):
    #     ax.pole(strike_cent[i], dip_cent[i], colors[i]+'s', markersize=1) # marker = '.'
    
    
    for i in range(0, np.size(dip)):
        if mode == 'plane':
            ax.plane(strike[i], dip[i], colors[labels[i]]+'-', linewidth=2)
        if mode == 'pole':
            ax.pole(strike[i], dip[i], color=tuple(colors[labels[i]]), marker='.', markersize=7)


    ax.grid()
    
    
    legend_patches = []
    for i in range(0, n_clusters):
      legend_patches.append(mpatches.Patch(color=colors[i], label="" + str(i+1)))
    ax.legend(handles=legend_patches, loc="upper right", bbox_to_anchor=(1.7, 0.8), title="Sets", fontsize=14, title_fontsize=14)
    
    plt.show()
    


def fracture_planes_plot(n_clusters, new_normals, new_centroids, model_dimension, plane_size=10):
    
    #model_center = model_dimension/2
    # points = np.asarray([[-1, 0,-1],[1,0,-1],[1,0,1], [-1,0,1]]) # base vertical square
    points = np.asarray([[-1, -1,0],[1,-1,0],[1,1,0], [-1,1,0]]) # base horizontal square
    points = points*plane_size
    
    rgb_colors = [[12, 12, 255], [24, 150, 33], [255,0,0], [127, 0, 110], [255, 0, 220], [255, 248, 56], [178, 0, 255]]
    rgb_colors = np.asarray(rgb_colors)/255    

    
    meshes = []
    
    for i in range(0, n_clusters):
        for j in range(0, np.shape(new_centroids[i])[0]):
            normal = new_normals[i][j]
            centroid = new_centroids[i][j]
            
            new_points = plane_rotate(points, normal, centroid)
            
            
            mesh = o3d.geometry.TriangleMesh()
            np_vertices = np.array(new_points)
            np_triangles = np.array([[0, 1, 2], [0, 2, 3]]).astype(np.int32)
            
            mesh.vertices = o3d.utility.Vector3dVector(np_vertices)
            mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
            mesh.paint_uniform_color(rgb_colors[i])
            mesh.compute_vertex_normals()
            meshes.append(mesh)
            
    mesh_frame = o3d.geometry.TriangleMesh.create_coordinate_frame(
    size=4, origin=[-4, -4, -4])
    
    meshes.append(mesh_frame)
            
    o3d.visualization.draw_geometries(meshes, mesh_show_back_face=True)
    
    return



def dfn_plot(intensity_values, offsets, cube_size, cut_out, reference_angle_color, model_dimension, only_walls=False):
    
    
    if only_walls==True:
        dim_x, dim_y, dim_z = model_dimension
        positions = np.reshape(offsets, (dim_x, dim_y, dim_z, 3))
        intensity = np.reshape(intensity_values, (dim_x, dim_y, dim_z))
        
        temp_intensity = []
        temp_positions = []
        for i in range(0, dim_x):
            for j in range(0, dim_y):
                for k in range(0, dim_z):
                    if i == 0 or i == dim_x-1:
                        temp_positions.append(positions[i, j, k])
                        temp_intensity.append(intensity[i, j, k])
                    else:
                        if j == 0 or j == dim_y-1:
                            temp_positions.append(positions[i, j, k])
                            temp_intensity.append(intensity[i, j, k])
                        if k == 0 or k == dim_z-1:
                            temp_positions.append(positions[i, j, k])
                            temp_intensity.append(intensity[i, j, k])                

        intensity_values = np.asarray(temp_intensity)
        offsets = np.asarray(temp_positions)
    

    voxel_colors = []
    # reference_angle_color = 180
    if np.size(intensity_values) > 1:
        hue = reference_angle_color-(intensity_values-np.min(intensity_values))/(np.max(intensity_values)- np.min(intensity_values))*reference_angle_color
    else:
        hue = [reference_angle_color-(intensity_values[0]*reference_angle_color)]
        
        
        
    hsv_colors =  np.rot90((hue, np.full(np.size(hue), 1), np.full(np.size(hue), 0.8)))

    voxel_colors = np.zeros(np.shape(hsv_colors))
    for i in range(0, np.shape(hsv_colors)[0]):
        voxel_colors[i] = colorsys.hsv_to_rgb(hsv_colors[i, 0]/360, hsv_colors[i, 1], hsv_colors[i, 2])
    
    
    # cube_list = []
    mesh_list = []
    i = 0
    for offset in offsets:
        
        if intensity_values[i] > cut_out:
        
            cube_definitions = Cube(cube_size, offset)
            vertices = cube_definitions.vertices
            
            # Cube faces
            faces = [[0,2,3], [2,6,3], [0,5,1], [0,2,5], [1,4,5], [4,7,5],
                      [4,3,6], [4,7,6], [0,1,4], [0,3,4], [2,5,7], [2,6,7]]
            
            np_triangles = np.array(faces).astype(np.int32)
            
            
            mesh = o3d.geometry.TriangleMesh()
            mesh.vertices = o3d.utility.Vector3dVector(vertices)
            mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
            mesh.paint_uniform_color(voxel_colors[i])
            # mesh.compute_vertex_normals()
            mesh_list.append(mesh)
            
        i+=1
        
    mesh_frame = o3d.geometry.TriangleMesh.create_coordinate_frame(
    size=4, origin=[-4, -4, -4])
    
    mesh_list.append(mesh_frame)
            
    o3d.visualization.draw_geometries(mesh_list, mesh_show_back_face=True)
    
    return



def compute_fracture_sets_fisher(strike, dip, plunge, bearing, n_clusters, labels):
    
    # Compute Fisher statistics for each set
    plane_stats = []
    for group in range(0, n_clusters):
        if np.size(plunge[np.where(labels == group)]) > 1:
            mean_strike, mean_dip, stats = fisher_plot(plunge[np.where(labels == group)],
                                                       bearing[np.where(labels == group)],
                                                       strike[np.where(labels == group)],
                                                       dip[np.where(labels == group)])
            
            plane_stats.append((mean_strike, mean_dip, stats[0], stats[1], stats[2]))
            #print(mpl.find_fisher_stats(strike[np.where(labels == group)], dip[np.where(labels == group)], measurements = 'poles', conf=95))
            #print(fisher_stats(plunge[np.where(labels == group)], np.asarray(bearing)[np.where(labels == group)]))
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='stereonet')
            ax.pole(strike[np.where(labels == group)],
                    dip[np.where(labels == group)], color="red")
            plt.show()
            plane_stats.append((strike[np.where(labels == group)],
                    dip[np.where(labels == group)], 1, 0, np.inf))
            print(strike[np.where(labels == group)], dip[np.where(labels == group)])
        print("-----")
    
    plane_stats = np.reshape(plane_stats, (np.shape(plane_stats)[0], 5))
    
    # Fisher statistics results
    # print('Fisher statistics')
    # print('Strike / Dip / r_value, Fisher_angle, K')
    # for i in range(0, np.shape(plane_stats)[0]):
    #     print(plane_stats[i])
    # print('\n K-means strike/dip centers')
    # print(np.asarray((strike_cent, dip_cent)).T)
        
    print("\\begin{table}[h!]\centering\\begin{tabular}{ccccc} \\toprule Set & Mean Strike/Dip & Fisher Angle & R-Value & K-Value \\\ \\midrule")
    for i in range(0, np.shape(plane_stats)[0]):
        
        print( str(i+1) + " & " + str(round(plane_stats[i][0])) + "°/" + str(round(plane_stats[i][1]))  + "° & " + 
              str("{:.2f}".format(plane_stats[i][3])) + "° & " + str("{:.2f}".format(plane_stats[i][2]))  + " & " + str("{:.2f}".format(plane_stats[i][4]))  + " \\\  ")  
    
    
    print("\\bottonrule \\end{tabular}  \caption{Directional statistics for each set showing mean direction/mean dip and computed dispersion values and Fisher angles.} \label{tab:clusters} \end{table}")

    return plane_stats
