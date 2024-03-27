import pandas as pd
import numpy as np
import os
from scipy.spatial import KDTree
from scipy.spatial import ConvexHull, convex_hull_plot_2d,  Delaunay
import networkx as nx
import openslide
import matplotlib.pyplot as plt
import math
import PIL
from PIL import Image, ImageDraw
PIL.Image.MAX_IMAGE_PIXELS = 100000 * 10000
import openslide
from hovernet_functions import load_preprocessed_img
import pickle
import math
import pdb
from pathlib import Path
import skimage
from skimage import measure




def extract_areas_convert_to_um2(df):
    areas_in_um = []
    indexes = [] #returns the df index for each area (which subpatch the nucleus was in)
    mpps = df["mpps"].iloc[0]
    for index, row in df.iterrows():
        areas = row["nuclei_area(pxs)"]
        if len(areas) > 0:
            for nuc_area in areas:
                indexes.append(index)
                areas_in_um.append(nuc_area * (mpps**2)) #bc it is an area measurement, px*px
    return areas_in_um, indexes


def find_contour_info(contour):
    xs = [i[0] for i in contour]
    ys = [i[1] for i in contour]
    width = max(xs) - min(xs)
    height = max(ys) - min(ys)
    contour_mapped = []
    for i in range(len(xs)):
        contour_mapped.append((xs[i]-min(xs), ys[i]-min(ys)))
    shape_mask = make_array_out_of_vertices(contour_mapped, int(width), int(height))
    regions = skimage.measure.regionprops(shape_mask)
    assert (len(regions) ==1)
    region = regions[0]
    r_ext = region["extent"]
    r_circ = (4 * math.pi * region["Area"]) / ((region["Perimeter"]) ** 2)
    axis_ratio = region["major_axis_length"] / region["minor_axis_length"]
    r_eccen = region["eccentricity"]    
    return [axis_ratio, r_ext, r_eccen, r_circ ]

def extract_contour_return_circularity(df):
    axis_ratio = []
    circ = []
    eccen = []
    extent = []
    mpps = df["mpps"].iloc[0]
    indexes = []
    for index, row in df.iterrows():
        contours = row["nuclei_contour_global(px)"]
        if len(contours) > 0:
            for contour in contours:
                indexes.append(index)
                contour_info  = find_contour_info(contour)
                axis_ratio.append(contour_info[0]) 
                extent.append(contour_info[1])
                eccen.append(contour_info[2])
                circ.append(contour_info[3])
                
    return axis_ratio, circ, eccen, extent, indexes



def combine_dataframes(subdf, nf):
    """adds information about nuclei to the fat algorithm
    MAKE SURE nucdf ONLY takes in "relevant nuclei!"
    returns a dataframe, which can be used as input for CCA_analysis"""

    assert [i == True for i in nf["is_a_relevant_nucleus"].to_list()] #MAKE SURE nucdf ONLY takes in "relevant nuclei!
    
    nucdf = nf.copy()
    centroids = nucdf["centroid_global(px)"].to_list()
    xs = [i[0] for i in centroids]
    ys = [i[1] for i in centroids]
    nucdf["centroid_global_x(px)"] = xs
    nucdf["centroid_global_y(px)"] = ys
    nucdf["id"] = np.arange(0, len(nucdf))
    
    #collects nuclei data
    nuclei_coords = []
    nuclei_area = []
    nuclei_contour = []
    nuclei_bbox = []
    nuclei_type = []
    nuclei_type_probability = []
    nuclei_number = []
    subdf = subdf.copy()
    
    for index, row in subdf.iterrows(): #each subpatch is a row
        x_coord = row["global_x_coords(px)"]
        y_coord = row["global_y_coords(px)"]
        subpatch_size = row["subpatch_size(px)"]
        mask_x, mask_y = create_mask(nucdf, x_coord, y_coord, subpatch_size)
        nucdf_slice = nucdf[mask_x & mask_y] #nuclei that are only in this location.
        nuclei_coords.append([nucdf_slice["centroid_global(px)"].to_list()])
        nuclei_bbox.append(nucdf_slice["bbox_global(px)"].to_list())
        nuclei_contour.append(nucdf_slice["contour_global(px)"].to_list())
        nuclei_type.append(nucdf_slice["type"].to_list())
        nuclei_type_probability.append(nucdf_slice["type_probability"].to_list())
        nuclei_area.append(nucdf_slice["area_nuc(pxs)"].to_list())
        nuclei_number.append(len(nucdf_slice["centroid_global(px)"].to_list()))
        assert (len(nucdf_slice["type"].to_list()) == len(nucdf_slice["area_nuc(pxs)"].to_list())) & (len(nucdf_slice["bbox_global(px)"].to_list()) == len(nucdf_slice["centroid_global(px)"].to_list())) #make sure same number of nuclei infos being added everywhere
        
        """some nuclei were being added doubled"""
        ids = nucdf_slice["id"].to_list()
        for i in ids:
            original = len(nucdf)
            nucdf = nucdf[nucdf.id != i]
            assert len(nucdf) == (original-1)
        
        
    # add nuclei info to the subdf. Nuclei in a subpatch are stored in lists
    #### that means nucleus#2 info is stored in nuclei_coords[2], nuclei_bbox[2] ...
    subdf["nuclei_coords_global(px)"] = nuclei_coords
    subdf["nuclei_area(pxs)"] = nuclei_area
    subdf["nuclei_contour_global(px)"] = nuclei_contour
    subdf["nuclei_bbox_global(px)"] = nuclei_bbox
    subdf["nuclei_type"] = nuclei_type
    subdf["nuclei_type_probability"] = nuclei_type_probability
    subdf["number_of_nuclei_in_subpatch"] = nuclei_number
    print("num nucleis ", subdf["number_of_nuclei_in_subpatch"].sum(), len(nf))
    assert subdf["number_of_nuclei_in_subpatch"].sum() <=  len(nf)
    return subdf

def create_mask(f, x, y, size):
    mask_x = (f["centroid_global_x(px)"] > x) & (f["centroid_global_x(px)"] < x + size)
    mask_y = (f["centroid_global_y(px)"] > y) & (f["centroid_global_y(px)"] < y + size)
    return mask_x, mask_y


"""PART 1: Global CCA"""


def global_connected_component_anlysis(df, max_distance, centroids): 
    """clusters nuclei together globally irrespective of patch boundaries"""
    
    area_removed_px = []
    cluster_l = []
    hulls_l = []
    
    print(f"         |_> building KDTree with {len(centroids)} centroids")
    tree = KDTree(data = centroids)
    #d_to_nearest_neighbor = distance_analysis(tree=tree, centroids=centroids, mpps=mpps) #uncomment for neighbor data
    components, graph_stats = cluster_tree(tree, max_distance, centroids)
    print(f"         |_> completed clustring, found {len(components)} clusters [not accounting for number of nuclei/cluster]")
    original_components = components#.copy() #commented out to save memory for now
    
    return components, tree, graph_stats


def extract_centroids(df):
    centroids = []
    for index, row in df.iterrows():
        nuclei_centroids = row["nuclei_coords_global(px)"]
        assert len(nuclei_centroids) == 1 #they are alwyas stored within a list
        nuclei_centroids = nuclei_centroids[0]
        if len(nuclei_centroids) > 0:
            for pt in nuclei_centroids:
                centroids.append(pt)
    assert df["number_of_nuclei_in_subpatch"].sum() == len(centroids)
    return centroids



def extract_centroids_by_type(df, nuc_type, dic):
    print("matching type: ", nuc_type)
    centroids = []
    probabilities = []
    bboxs_l = []
    contours_l = []
    for index, row in df.iterrows():
        nuclei_centroids = row["nuclei_coords_global(px)"]
        nuc_types = row["nuclei_type"]
        probs = row["nuclei_type_probability"]
        bboxs = row['nuclei_bbox_global(px)']
        contours = row['nuclei_contour_global(px)']
        
        assert len(nuclei_centroids) == 1 #they are alwyas stored within a list
        nuclei_centroids = nuclei_centroids[0]
        if len(nuclei_centroids) > 0:
            for i in range(len(nuclei_centroids)):
                pt = nuclei_centroids[i]
                pt_type = dic[nuc_types[i]]
                if pt_type == nuc_type: #if is the correct type
                    centroids.append(pt)
                    probabilities.append(probs[i])
                    bboxs_l.append(bboxs[i])
                    contours_l.append(contours[i])
    return centroids, probabilities, bboxs_l, contours_l



    
def cluster_tree(tree, max_distance, centroids):
    neighbor_indices = tree.query_ball_tree(tree, r=max_distance) ## Find all pairs of points whose distance is at most r
    G = nx.Graph() #create a graph with the edges being the connection to each neighbor
    for i in range(len(neighbor_indices)):
        for j in neighbor_indices[i]:
            G.add_edge(i, j)
    components = list(nx.connected_components(G)) ## list of indices that "belong" together 
    num_components = nx.number_connected_components(G)
    clustered_nucs = [] # [ [(x,y), (x,y)], [(x,y)], ... ]. len list = total num of components, 
    for component in components:
        nucs_in_component = []
        for nuc in component: #nuc is the index of the nucleus in respect to the original list
            nucs_in_component.append(centroids[nuc])
        clustered_nucs.append(nucs_in_component)
    assert num_components == len(clustered_nucs)
    
    avg_clustering_coeffs, clustering_coeffs, num_neigh, non_neigh_num, avg_dis_to_neighbor, avg_dis_to_nonneighbor, num_comm_neighs, degree, is_in_cluster = graph_statistics(G, centroids, clustered_nucs)
    
    return clustered_nucs,  [avg_clustering_coeffs, clustering_coeffs, num_neigh, non_neigh_num, avg_dis_to_neighbor, avg_dis_to_nonneighbor, num_comm_neighs, degree, is_in_cluster]

def graph_statistics(G, rel_nuclei_centroids, components):
    neighbors = [] #lst with [node, [(neighbor_coordsx,neighbor_coordsy)]]
    num_neigh = []
    neigh_indexes = [] #lst with [node, [neighbor_indexes]]
    non_neigh_num =[]
    non_neighbors =[] #lst with [node, [(neighbor_coordsx,neighbor_coordsy)]]
    degree = []
    
    for node in list(G.nodes):
        neigh = [rel_nuclei_centroids[n] for n in nx.all_neighbors(G, node)]
        num_neigh.append(len(neigh) -1 ) ## itself is always considered a neighbor!!!
        neighbors.append([rel_nuclei_centroids[node], neigh])
        neigh_indexes.append([node, [n for n in nx.all_neighbors(G, node)]])
        non_neigh = [rel_nuclei_centroids[n] for n in nx.non_neighbors(G, node)]
        non_neigh_num.append(len(non_neigh)) #itself is not included
        non_neighbors.append([rel_nuclei_centroids[node], non_neigh])
        degree.append(G.degree(node))
    
    avg_dis_to_neighbor = return_avg_distance_between_node_and_pts(neighbors)
    avg_dis_to_nonneighbor = return_avg_distance_between_node_and_pts(non_neighbors)

    #local_con = local_connectivity(G) # returns {node, {connectivity to every other node}} #takes forever to run

    num_comm_neighs = []
    for n in neigh_indexes:
        node = n[0]
        neighs = n[1]
        avg_num_common_neighbors = []
        for neigh in neighs:
            if not(node == neigh):
                avg_num_common_neighbors.append(len([n for n in nx.common_neighbors(G, node, neigh)]))
        num_comm_neighs.append(np.mean(avg_num_common_neighbors))
    
    clustering_coeffs = nx.clustering(G)
    avg_clustering_coeffs = nx.average_clustering(G)
    
    is_in_cluster = []
    for point in rel_nuclei_centroids:
        is_in_c = find_pt_in_a_cluster(point, components)
        is_in_cluster.append(is_in_c)

    return avg_clustering_coeffs, clustering_coeffs, num_neigh, non_neigh_num, avg_dis_to_neighbor, avg_dis_to_nonneighbor, num_comm_neighs, degree, is_in_cluster

def find_pt_in_a_cluster(pt, clusters):
    for cluster in clusters:
        for nuc in cluster:
            if (nuc == pt) & (len(cluster) > 2):
                return True
    return False

def return_avg_distance_between_node_and_pts(lst):
    avg_dd = []
    for node in lst:
        n = node[0]
        pts = node[1]
        dd = []
        for i in pts:
            d = math.dist(n,i)
            dd.append(math.dist(n, i))
        if len(dd) == 0: #no points
            avg_dd.append(np.nan)
        else:
            avg_dd.append(np.mean(dd))
    assert len(avg_dd) == len(lst)
    return avg_dd


"""PART 2: Process components"""

def find_patchsize(df):
    keys = np.unique(df["original_key"].to_list())
    coords = df[df["original_key"]==keys[0]]["global_x_coords(px)"].to_list()#has lots of subpatches
    return max(coords) - min(coords) # was always process as a square, which is why this applies to y too


def run_CCA(base_path, fat_path, max_distance, experiment_name, removed_df_label):
    """run CCA_analysis on imgs that have been processed by hovernet and algorithm"""
    processed = []
    
    #make folders to save information
    Path(f"{base_path}{experiment_name}/cluster_info/").mkdir(parents=True, exist_ok=True)
    Path(f"{base_path}{experiment_name}/CCA_dfs/").mkdir(parents=True, exist_ok=True)
    print(len(os.listdir(base_path)))
    for img in os.listdir(base_path):
        
        name = img.replace("_fat_and_nucs.pkl","")
        print(name)
        if img == experiment_name:
            continue
        if ((os.path.exists(f"{base_path}{img}")) & (not(os.path.exists(f"{base_path}{experiment_name}/CCA_dfs/{name}_{experiment_name}_global_CCA.pkl")))):
            
            processed.append(img)
            
            print("working on img:", img)
            print("---> experiment: ", experiment_name)
            
            

            #algo_hvnet_df = pd.read_pickle(f"{base_path}{img}_fat_and_nucs.pkl")
            algo_hvnet_df = pd.read_pickle(f"{base_path}{img}")
            print("---> loaded hovernet_and_algo_df")

            fat_loc =  f"{fat_path}{name}/saved_patches/fat_patches/"
            clusters_loc =  f"{base_path}{experiment_name}/cluster_info/"

            max_distance_px = int(max_distance / algo_hvnet_df["mpps"].iloc[0]) ## convert from_um to px
      
            patch_size_px = algo_hvnet_df["original_patch_size(px)"].iloc[0]
        
    
            print("---> max_distance in um: ", max_distance, "max distance in px: ", max_distance_px)
            print("---> patch_size in px: ", patch_size_px)
            
            
            centroids =  extract_centroids(df=algo_hvnet_df) # globally cluster nuclei!
            
            components, tree, graph_stats = global_connected_component_anlysis(df=algo_hvnet_df, max_distance = max_distance_px, centroids=centroids)
            
            print("finished finding components")
            with open(f"{clusters_loc}{name}_{experiment_name}_clusters.pkl", 'wb') as j:
                pickle.dump(components, j)
            
            with open(f"{clusters_loc}{name}_{experiment_name}_graph_stats.pkl", 'wb') as j:
                pickle.dump(graph_stats, j)
            
            
            
            cca_df, hull_df = process_components(components=components, algo_hvnet_df=algo_hvnet_df, name=img, path_to_fat_masks=fat_loc, clusters_loc=clusters_loc, data_saver_mode=False, patch_size_px=patch_size_px, removed_df_label=removed_df_label, experiment_name=experiment_name)
            print("saving final dfs")
            cca_df.to_pickle(f"{base_path}{experiment_name}/CCA_dfs/{name}_{experiment_name}_global_CCA.pkl")
            
                                        
    

def get_patch_data(path, patch_size):
    r = {}
    for i in os.listdir(path):
        patch_key = int(i.split("_")[0]) #check this :D
        coords_x = int(i.split("_")[1])
        coords_y = int(i.split("_")[2].replace(".tif",""))
        if not(patch_key in r):
                r[patch_key] = {"global_x_coords(px)":[], "global_y_coords(px)":[],"patch_xmax":[], "patch_ymax":[]}
        patch_r = r[patch_key]
        patch_r["global_x_coords(px)"].append(coords_x)
        patch_r["global_y_coords(px)"].append(coords_y)
        patch_r["patch_xmax"].append(coords_x + patch_size)
        patch_r["patch_ymax"].append(coords_y + patch_size)
    return r
    
        
        

def process_components(components, algo_hvnet_df, name, path_to_fat_masks, clusters_loc, data_saver_mode, patch_size_px, removed_df_label, experiment_name):
    
    """1. create hull masks for each component & return data about the hull
       2. check which hulls overlap over multiple subpatches (irregardless of original patch)
       3. load the fat patch
       3. select hulls in that patch and make one mask out of all the little hull masks
       4. overlay the mask with the fat_patch to remove fat from the area. 
    
    """

    
    hull_df = make_hulls(components)
    #hull_df.to_pickle(f"{clusters_loc}{name}_{experiment_name}_hull_df.pkl")
   
    
    print("adding overlap data")

    patch_data = get_patch_data(path_to_fat_masks, patch_size_px)
            
    hulls_in_patch = {}
    hulls_patches = []
    for index, hrow in hull_df.iterrows():
        bbox = hrow["hull_bbox_global(px)"]
        bbox_x = int(bbox[0][0]) #is smallest x
        bbox_y = int(bbox[0][1] )#is smallest y
        bbox_w = bbox[1] #width
        bbox_h = bbox[2] #height
        bbox_xmax = bbox_x + bbox_w
        bbox_ymax = bbox_y + bbox_h
        #print("bbox_info: ",bbox_x, bbox_y, bbox_xmax, bbox_ymax)
        hull_patches = set()
        for original_key, patch_hulls in patch_data.items():
            #print(original_key, patch_hulls)
            assert len(patch_hulls["global_x_coords(px)"]) ==1
            assert len(patch_hulls["global_y_coords(px)"]) ==1
            assert len(patch_hulls["patch_xmax"]) ==1
            assert len(patch_hulls["patch_ymax"]) ==1
            
            subpatchx = patch_hulls["global_x_coords(px)"][0]
            subpatchy = patch_hulls["global_y_coords(px)"][0]
            subpatchxmax = patch_hulls["patch_xmax"][0]
            subpatchymax = patch_hulls["patch_ymax"][0]
            
            if ((subpatchx < bbox_xmax) & (subpatchxmax > bbox_x)) & (subpatchy < bbox_ymax) & (subpatchymax > bbox_y):
                hull_patches.add(original_key)
        hulls_patches.append(list(hull_patches)) 
        
        for patch_key in hull_patches:
            if not(patch_key in hulls_in_patch):
                hulls_in_patch[patch_key] = {"masks":[], "bboxs":[], "areas":[], "patch_x":[], "patch_y":[],"patch_xmax":[], "patch_ymax":[]}
            patch_hulls = hulls_in_patch[patch_key]
            patch_hulls["masks"].append(hrow["hull_masks"])
            patch_hulls["bboxs"].append(hrow["hull_bbox_global(px)"])
            patch_hulls["areas"].append(hrow["hull_area(pxs)"])
            patch_hulls["patch_x"].append(patch_data[patch_key]["global_x_coords(px)"][0])
            patch_hulls["patch_y"].append(patch_data[patch_key]["global_y_coords(px)"][0])
            patch_hulls["patch_xmax"].append(patch_data[patch_key]["patch_xmax"])
            patch_hulls["patch_ymax"].append(patch_data[patch_key]["patch_ymax"])
            
    hull_df["patch_keys_for_hull"] = hulls_patches  
    
    
    #hull_df.to_pickle(f"{clusters_loc}{name}_{experiment_name}_hull_df.pkl") #overwrites previous save, has more data :)
    """
    with open(f"{clusters_loc}{name}_hulls_in_patch.pkl", 'wb') as j:
        pickle.dump(hulls_in_patch, j)
    
    """
    print("processing each key")
  

    removed_area = []
    patch_key_l = []
    
    for original_key, patch_hulls in hulls_in_patch.items():
        if original_key >= 0: #keep for testing
            patch, patch_x_global, patch_y_global = load_preprocessed_img(loc=path_to_fat_masks, wanted_key=original_key)
            print("original_key", original_key)
            print("patch_shape", patch.shape)
            masks = patch_hulls["masks"]
            bboxs = patch_hulls["bboxs"]
            patch_xmax = patch_hulls["patch_xmax"][0][0]
            patch_ymax = patch_hulls["patch_ymax"][0][0]

            """
            patch_info = algo_hvnet_df[algo_hvnet_df["original_key"]== original_key]

            patch_x_global_2 = min(patch_info["global_x_coords(px)"].to_list())
            patch_y_global_2 = min(patch_info["global_y_coords(px)"].to_list())
            patch_xmax2 = max(patch_info["global_x_coords(px)"].to_list()) + patch_info["subpatch_size(px)"].iloc[0]
            patch_ymax2 = max(patch_info["global_y_coords(px)"].to_list()) + patch_info["subpatch_size(px)"].iloc[0]

            patch_xmax = patch_x_global + patch_size_px
            patch_ymax = patch_y_global + patch_size_px

            assert patch_xmax2 == patch_xmax
            assert patch_ymax2 == patch_ymax

            assert patch_x_global == patch_x_global_2
            assert patch_y_global == patch_y_global_2

            """
            #print(f"key: {original_key}, patch_x_global {patch_x_global}, patch_y_global {patch_y_global}, patch_xmax {patch_xmax}, patch_ymax {patch_ymax}")

            for bbox in bboxs: #check that each bbox really is in patch
                bbox_x = int(bbox[0][0]) #is smallest x
                bbox_y = int(bbox[0][1] )#is smallest y
                bbox_w = bbox[1] #width
                bbox_h = bbox[2] #height
                bbox_xmax = bbox_x + bbox_w
                bbox_ymax = bbox_y + bbox_h
                bboxxs = np.arange(bbox_x, (bbox_x + bbox_w +1))
                bboxys = np.arange(bbox_y, (bbox_y + bbox_h +1))
                patch_xs = np.arange(patch_x_global, patch_xmax + 1)
                patch_ys = np.arange(patch_y_global, patch_ymax + 1)
                overlap_x = return_overlap(patch_xs, bboxxs)
                overlap_y = return_overlap(patch_ys, bboxys)
                if len(overlap_x) <= 0 :
                    print("bbox info:" , bbox_x, bbox_y, bbox_w, bbox_h, bbox_xmax, bbox_ymax)
                    print("patch info:", patch_xs, patch_ys)
                assert len(overlap_x) > 0
                assert len(overlap_y) > 0
            #print("making big hull mask for patch: ", original_key, patch_x_global, patch_y_global)


            hulls_mask = make_template(masks=masks, bboxs=bboxs, glob_x=patch_x_global, glob_y=patch_y_global, max_globx=patch_xmax, max_globy=patch_ymax, size=patch_size_px, canvas= np.ones_like(patch))

            assert (hulls_mask.shape == patch.shape)

            tm = patch.copy()
            hm = hulls_mask.copy()
            remove = np.logical_or(hm, np.invert(tm)) #overlapping fat_areas removed

            """
            if np.count_nonzero(tm==0) > 0 :                          
                fig, ax = plt.subplots(ncols = 3)
                ax[0].imshow(tm)
                ax[1].imshow(hm)
                ax[2].imshow(remove)
                ax[1].axis("off")
                ax[0].axis("off")
                ax[2].axis("off")
                ax[0].set_title(f"fat, key {original_key}")
                ax[1].set_title("hull_mask")
                ax[2].set_title("removed area")
                plt.show()
            """


            area_removed_by_subpatch = process_subpatches(remove, algo_hvnet_df["subpatch_size(px)"].iloc[0])
            removed_area.append([original_key, area_removed_by_subpatch])
            patch_key_l.append(original_key)
        
    CCA_df = add_data_to_df(algo_hvnet_df, patch_key_l, removed_area, removed_df_label)
    return CCA_df, hull_df


def add_data_to_df(df, patches, removed_area_l, label):
    area_removed = []
    for patch in np.unique(np.array(df["original_key"].to_list())):
        if patch in patches: #make sure same type
            for i in removed_area_l:
                if i[0] == patch:
                    for j in i[1]:
                        area_removed.append(j[1])
        if not(patch in patches): #for ex bc no nuclei in patch
            df_slice = df[df["original_key"] == patch]
            for index, row in df_slice.iterrows():
                area_removed.append(0)
    df[f"{label}"] = area_removed
    return df



def return_overlap(l1,l2):
    l1 = set(l1)
    l2 = set(l2)
    return sorted(list(l1.intersection(l2)))



def make_template(masks, bboxs, glob_x, glob_y, max_globx, max_globy, size, canvas):
    
    assert (canvas.shape[0]== size) & (canvas.shape[1]== size) #is a square
    assert len(masks) == len(bboxs)
    
    
    for i in range(len(masks)):
        mask = masks[i]
        bbox = bboxs[i]
        bx = int(bbox[0][0])
        by = int(bbox[0][1])
        bw = int(bbox[1])
        bh = int(bbox[2])
    
        canvas_minx = int(glob_x)
        canvas_miny = int(glob_y)
        canvas_maxx = max_globx
        canvas_maxy = max_globy
        canvas_x = np.arange(canvas_minx, canvas_maxx +1 ) #lst x_coords global
        canvas_y = np.arange(canvas_miny, canvas_maxy  +1) #lst y_coords global

        mask_minx = bx
        mask_miny = by 
        mask_maxx = bx + bw
        mask_maxy = by + bh
        mask_x = np.arange(mask_minx, mask_maxx + 1) #lst x_coords global
        mask_y = np.arange(mask_miny, mask_maxy + 1) #lst y_coords global


        assert (canvas_maxx - canvas_minx) == size
        assert (canvas_maxy - canvas_miny) == size
        
        xs = return_overlap(canvas_x, mask_x) #overlapping x coords
        ys = return_overlap(canvas_y, mask_y) #overlapping y coords

        if (len(xs) ==0) or len(ys)==0:
            print(f"mask num: {i}")
            print(f"canvas minx: {canvas_minx}, canvas max x {canvas_maxx}")
            print(f"canvas miny: {canvas_miny}, canvas maxy {canvas_maxy}")

            print(f"mask minx: {mask_minx}, mask maxx {mask_maxx}")
            print(f"mask miny: {mask_miny}, mask maxy {mask_maxy}")

            print(f"overlap x {xs}")
            print(f"overlap y {ys}")
            
        
        assert len(xs) > 0 
        assert len(ys) > 0
        
        # map overlaps to canvas scale
        min_lcx = min(xs) - glob_x 
        max_lcx = max(xs) - glob_x
        min_lcy = min(ys) - glob_y
        max_lcy = max(ys) - glob_y
        
        # map overlaps to mask scale
        min_lmx = min(xs) - bx # should be greater than 0
        max_lmx = max(xs) - bx
        min_lmy = min(ys) - by # should be greater than 0
        max_lmy = max(ys) - by

        #print(min_lcx,max_lcy, min_lcx,max_lcx, min_lmy,max_lmy, min_lmx,max_lmx)

        if len(mask.shape)>0: 
            canvas[min_lcy:max_lcy, min_lcx:max_lcx] = np.invert(np.logical_and(canvas[min_lcy:max_lcy, min_lcx:max_lcx], mask[min_lmy:max_lmy, min_lmx:max_lmx]))
    return canvas
            
            
def process_subpatches(array, subpatch_size):
    sp = 0
    l = []
    array = np.invert(array) ## areas that should be removed ==1 now
    for y in range(0, array.shape[1], subpatch_size): 
        for x in range(0, array.shape[0], subpatch_size):
            subpatch = array[y:y+subpatch_size, x:x+subpatch_size]
            area_removed = np.sum(subpatch)
            l.append((sp, area_removed))
            sp += 1
    return l



def make_hulls(components):
    """ 
    Components is a list of list of list(nuc_x, nuc_y) global position
    
    returns a dataframe for each hull:
      -lst_of_hull_vertices (global)
      -area_of_clusters (px2)
      -hull_centroids_l (global)
      -masks_of_hulls_l (centered to (0,0) in top left corner)
      -bounding_box_l (global)
      -which subpatches & original patches does hull overlap with
      
      - a component must have at least 2 distinct points! #hovernet sometimes records multiple nuclei in the same location. Alternatively "clean" hovernet dataframe by removing duplicates...
    """
    print("making hulls")
    verts = [] ## for plotting
    areas = []
    hull_centroids_l = []
    masks = []
    bboxs = []
    cluster_num = 0
    for cluster in components:
        #print("cluster number: ", cluster_num, " length: ", len(cluster))
        cluster_num +=1
        
        distinct_pts = list(set(tuple(sorted(pt)) for pt in cluster)) 
        
        """
        #just comment this part out
        if len(distinct_pts) <= 2: 
            areas.append(np.nan) # not a cluster, b/c too few nuclei
            hull_centroids_l.append(np.nan)
            bboxs.append(np.nan)
            masks.append(np.nan)
            verts.append(np.nan)
         """
        
            
        if len(distinct_pts) > 2:
            hull =  ConvexHull(cluster)
            areas.append(hull.volume) # see docs, here hull.volume = area
            vert = hull.vertices #is an index
            
            mapped_verts = [cluster[v] for v in vert] #lst of global coords of vertices
            verts.append(mapped_verts)
            
            hull_centroid = centroid_poly(np.array(mapped_verts)) #GET RID OF THIS
            hull_centroids_l.append(hull_centroid)
            
            bbox, min_x, min_y, w, h = bounding_box(mapped_verts) #min_x, min_y=top left corner of bb
            bboxs.append(bbox)
            verts_local = scaled_down(mapped_verts, min_x, min_y) #scale to 0 to make mask
            mask = make_array_out_of_vertices(verts_local, width=w, height=h)
            masks.append(mask)
            
        data = {"global_hull_verts(px)":verts, "hull_area(pxs)":areas, "hull_masks":masks, "hull_bbox_global(px)":bboxs, "hull_centroids": hull_centroids_l}
        df = pd.DataFrame(data)
    return df

def scaled_down(lst, min_x, min_y):
    r = []
    for i in lst:
        r.append(( i[0] - min_x, i[1] - min_y))
    return r

def bounding_box(lst):
    """takes a lst of a lst of x and y coordinates.
    returns min(x), min(y), width and height of bbox"""
    xs = []
    ys = []
    for i in lst:#separate into xs and ys GLOBAL
        xs.append(i[0])
        ys.append(i[1])
    w = int(np.max(xs)) - int(np.min(xs)) #width
    h = int(np.max(ys)) - int(np.min(ys)) #height
    return [(np.min(xs), np.min(ys)), w, h], np.min(xs), np.min(ys), w, h
        
def make_array_out_of_vertices(lst, width, height):
    """takes a list of coordinates and draws a polygon using those points
    https://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask
    Isaac Sutherland, Sep 17, 2010
    """
    polygon = lst
    img = Image.new('L', (width, height), 0)
    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
    mask = np.array(img)
    return mask

def centroid_poly(poly): # from example on stackoverflow
    """https://stackoverflow.com/questions/31562534/scipy-centroid-of-convex-hull
    Answered by Alec Day, Sept 11 2019
    Edited Jul 18, 2021 by klwire
    """
    T = Delaunay(poly).simplices
    n = T.shape[0]
    W = np.zeros(n)
    C = 0
    for m in range(n):
        sp = poly[T[m, :], :]
        W[m] = ConvexHull(sp).volume
        C += W[m] * np.mean(sp, axis=0)
    return C / np.sum(W)

    
    
    
    
    
    
    
    
    
    
    
    


    
    
    

    
### Cluster by group



def overlay_masks(keys, loc, hull_df1, hull_df2, patch_size_px, algo_hvnet_df):
    removed_area = []
    patch_key_l = []
    
    for original_key in keys:
        patch, patch_x_global, patch_y_global = load_preprocessed_img(loc=loc, wanted_key=original_key) #load the img
        masks = [] #collect relevant masks
        bboxs=[]
        for index, row in hull_df1.iterrows():
            key_lst = row["hull_original_key"]
            if original_key in key_lst:
                masks.append(row["hull_masks"])
                bboxs.append(row["hull_bbox_global(px)"])
                
        for index, row in hull_df2.iterrows():
            key_lst = row["hull_original_key"]
            if original_key in key_lst:
                masks.append(row["hull_masks"])
                bboxs.append(row["hull_bbox_global(px)"])

        
        hulls_mask = make_template(masks, bboxs, patch_x_global, patch_y_global, patch_size_px, canvas= np.ones_like(patch))
        assert (hulls_mask.shape == patch.shape)
        tm = patch.copy()
        hm = hulls_mask.copy()
        remove = np.logical_or(hm, np.invert(tm)) #overlapping fat_areas removed
        area_removed_by_subpatch = process_subpatches(remove, algo_hvnet_df["subpatch_size(px)"].iloc[0])
        removed_area.append([original_key, area_removed_by_subpatch])
        patch_key_l.append(original_key)
  
    return removed_area, patch_key_l


def add_hull_data(df, hull_df, nuc_type):
    r = []
    for patch in np.unique(np.array(df["original_key"].to_list())):
        print(patch)
        print(type(patch))
        sliced = hull_df[hull_df["hull_original_key"] == patch]
        for index, row in hull_df.iterrows():
            if row["hull_original_key"] == patch:
                r.append(row["hull_area(pxs)"])
    df[nuc_type] = r
    return df


def cluster_by_group (base_path, fat_path, patch_size_um, max_distance_1, max_distance_2, experiment_name, type_dic, nuc_type_1, nuc_type_2, processed):
    to_process = [x for x in os.listdir(base_path) if x not in processed]
    for img in to_process:
        if (os.path.isdir(base_path + img)) & (os.path.exists(f"{base_path}{img}/dataframes/{img}_fat_and_nucs.pkl")):
            print("working on img: ", img)
            print("experiment: ", experiment_name)
            fat_loc =  f"{fat_path}{img}/saved_patches/fat_patches/"
            clusters_loc =  f"{base_path}{img}/cluster_info/clusters/"
            overlap_loc = f"{base_path}{img}/cluster_info/overlap_lsts/"

            algo_hvnet_df = pd.read_pickle(f"{base_path}{img}/dataframes/{img}_fat_and_nucs.pkl")

            max_distance_1 = int(max_distance_1 / algo_hvnet_df["mpps"].iloc[0]) ## convert from_um to px
            max_distance_2 = int(max_distance_2 / algo_hvnet_df["mpps"].iloc[0]) ## convert from_um to px
            patch_size_px = int(patch_size_um / algo_hvnet_df["mpps"].iloc[0]) ## convert from_um to px
            #print("max_distance in um: ", max_distance, "max distance in px: ", max_distance_px)
            print("patch_size_um: ", patch_size_um, "patch_size in px: ", patch_size_px)
            print("loaded hovernet_and_algo_df")
          
            centroids_1, probabilities,bboxes, contours = extract_centroids_by_type(algo_hvnet_df, nuc_type_1, type_dic)
            centroids_2, probabilities, bboxes, contours= extract_centroids_by_type(algo_hvnet_df, nuc_type_2, type_dic)
            print(f"found {len(centroids_1)} nuclei of type {nuc_type_1}")
            print(f"found {len(centroids_2)} nuclei of type {nuc_type_2}")

            components_1, tree_1, graph_stats_1 = global_connected_component_anlysis(df=algo_hvnet_df, max_distance = max_distance_1, centroids=centroids_1)
            components_2, tree_2, graph_stats_2 = global_connected_component_anlysis(df=algo_hvnet_df, max_distance = max_distance_2, centroids=centroids_2)

            print("making hulls and adding data about overlapping patches")

            CCA_df1, hull_df1 = process_components(components_1, algo_hvnet_df, img, fat_loc, clusters_loc, overlap_loc, False, patch_size_px, removed_df_label = f"area_{nuc_type_1}_removed(px)")
            CCA_df2, hull_df2 = process_components(components_2, algo_hvnet_df, img, fat_loc, clusters_loc, overlap_loc, False, patch_size_px, removed_df_label = f"area_{nuc_type_2}_removed(px)")
            ## (area1 + area_2) - area_removed should be area of overlapping space 

            rm1 = CCA_df1[f"area_{nuc_type_1}_removed(px)"].to_list()
            rm2 = CCA_df2[f"area_{nuc_type_2}_removed(px)"].to_list()

            #select patches to load
            keys = set() # which original patches contain hulls
            for index, row in hull_df1.iterrows():
                keys.update(row["hull_original_key"])

            for index, row in hull_df2.iterrows():
                keys.update(row["hull_original_key"])

            # the hulls could overlap
            removed_area, patch_key_l = overlay_masks(keys, fat_loc, hull_df1, hull_df2, patch_size_px, algo_hvnet_df)
            CCA_df = add_data_to_df(algo_hvnet_df, patch_key_l, removed_area, label = f"total_area_removed(px)")
            CCA_df[f"area_{nuc_type_1}_removed(px)"] = rm1
            CCA_df[f"area_{nuc_type_2}_removed(px)"] = rm2
            CCA_df[f"overlapping_{nuc_type_1}_{nuc_type_2}_area(px)"] =  (CCA_df[f"area_{nuc_type_1}_removed(px)"] + CCA_df[f"area_{nuc_type_2}_removed(px)"]) - CCA_df["total_area_removed(px)"] 
            CCA_df.to_pickle(f"{base_path}{img}/dataframes/{img}_{experiment_name}_global_CCA.pkl")
