import pandas as pd
import numpy as np
import os

def find_20x_fov(mpps):
    # find correct patchsize for a 20x field of view. Assumes FOV is a square.
    patchsize_20xfov_in_px = int(1120/mpps) # we found a fov of 1120 um for a fov of 20x in our microscopes
    return patchsize_20xfov_in_px


from CCA_nuclei import bounding_box
def make_lst_of_centorids(components):
    #for each component, if there are more than 3 nuclei find the bbox and calcualte the centroid of the bbox. 
    # test = [[[0,0],[0,1],[1,0],[1,1]]]. centroid_bbox = lst_of_centorids(test)
    lst_of_centroids_of_components = []
    for component in components:
        if len(component) >=3 :
            bbox_info, bbox_minx, bbox_miny, w, h = bounding_box(component)
            centroid_bbox = [(bbox_minx + w/2), (bbox_miny + h/2)]
            lst_of_centroids_of_components.append(centroid_bbox)
    return lst_of_centroids_of_components



def find_foci_per_patch(lst_of_centroids, xmin, xmax, ymin, ymax):
   
    r = []
    for x, y in lst_of_centroids:
        if ((xmin <= x) and (x <= xmax)) and ((ymin <= y ) and( y <= ymax)):
            #print("found a match")
            r.append((x, y))

    return r
  
    
def run_Infl_by_foci_for_all_imgs(path, processed, experiment_name, save_loc):
    
    for img in os.listdir(path):
        print(img)
        if img.endswith(".pkl"):
            ## find file name
            if img.endswith(f"{experiment_name}_global_CCA.pkl"):
                name = img.replace(f"{experiment_name}_global_CCA.pkl","")
            
            if not(name in processed):
                print(f"working on {name}")
                processed.append(name)
                global_CCA = pd.read_pickle(path + name + f"{experiment_name}_global_CCA.pkl")
                clusters = pd.read_pickle(path.replace("CCA_dfs", "cluster_info") +  name + f"{experiment_name}_clusters.pkl")
                
                cluster_centroids = make_lst_of_centorids(clusters)
                
                
                mpps = global_CCA["mpps"].iloc[0] #always the same
                subpatch_size = global_CCA["subpatch_size(px)"].iloc[0]#always the same
                start_x = int(min(global_CCA["global_x_coords(px)"]))
                start_y = int(min(global_CCA["global_y_coords(px)"]))
                end_x = int(max(global_CCA["global_x_coords(px)"])) + subpatch_size
                end_y = int(max(global_CCA["global_y_coords(px)"])) + subpatch_size
              
                
                num_centroids_in_patch = []
                coords_x_l = []
                coords_y_l = []
                
                patchsize_20xfov = find_20x_fov(mpps) ### replace with a different function for another fov. Remember to change naming of final pkl file (last line of function)
                print("working with patchsize in px", patchsize_20xfov)
                for coords_y in range(start_y, end_y, patchsize_20xfov):
                    for coords_x in range(start_x, end_x, patchsize_20xfov):
                        
                        centroids_in_patch = find_foci_per_patch(lst_of_centroids = cluster_centroids, xmin=coords_x, xmax = coords_x + patchsize_20xfov, ymin=coords_y, ymax=coords_y + patchsize_20xfov)
                        num_centroids_in_patch.append(len(centroids_in_patch))
                        coords_x_l.append(coords_x)
                        coords_y_l.append(coords_y)
                
                df_of_num_centroids_per_fov = pd.DataFrame(data={"coords_x(px)": coords_x, "coords_y(px)": coords_y, "num_centroids_in_patch":num_centroids_in_patch })
                print(f"saving {name}")    
                df_of_num_centroids_per_fov.to_pickle(f"{save_loc}/{experiment_name}/num_foki_per_20fov/{name}_foki_per_20xfov.pkl") ## change name if not using 20x fov
                
                
                
                
                