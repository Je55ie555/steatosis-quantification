import pandas as pd
import os
import json
from PIL import Image
Image.MAX_IMAGE_PIXELS = 100000 * 10000
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import subprocess
from scipy.spatial import ConvexHull



def create_mask(nucdf, x, y, size):
    f = nucdf.copy()
    centroids = nucdf["centroid_global(px)"].to_list()
    xs = [i[0] for i in centroids]
    ys = [i[1] for i in centroids]
    f["centroid_global_x(px)"] = xs
    f["centroid_global_y(px)"] = ys
    mask_x = (f["centroid_global_x(px)"] >= x) & (f["centroid_global_x(px)"] <= x + size)
    mask_y = (f["centroid_global_y(px)"] >= y) & (f["centroid_global_y(px)"] <= y + size)
    return mask_x, mask_y, f


def combine_dataframes(subdf, nucdf):
    """adds information about nuclei to the fat algorithm
    MAKE SURE nucdf ONLY takes in "relevant nuclei!"
    returns a dataframe, which can be used as input for CCA_analysis"""

    assert [i == True for i in nucdf["is_a_relevant_nucleus"].to_list()] #MAKE SURE nucdf ONLY takes in "relevant nuclei!
    
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
        mask_x, mask_y, nucdf = create_mask(nucdf, x_coord, y_coord, subpatch_size)
        nucdf_slice = nucdf[mask_x & mask_y] #nuclei that are only in this location.
        nuclei_coords.append([nucdf_slice["centroid_global(px)"].to_list()])
        nuclei_bbox.append(nucdf_slice["bbox_global(px)"].to_list())
        nuclei_contour.append(nucdf_slice["contour_global(px)"].to_list())
        nuclei_type.append(nucdf_slice["type"].to_list())
        nuclei_type_probability.append(nucdf_slice["type_probability"].to_list())
        nuclei_area.append(nucdf_slice["area_nuc(pxs)"].to_list())
        nuclei_number.append(len(nucdf_slice["centroid_global(px)"].to_list()))

    # add nuclei info to the subdf. Nuclei in a subpatch are stored in lists
    #### that means nucleus#2 info is stored in nuclei_coords[2], nuclei_bbox[2] ...
    subdf["nuclei_coords_global(px)"] = nuclei_coords
    subdf["nuclei_area(pxs)"] = nuclei_area
    subdf["nuclei_contour_global(px)"] = nuclei_contour
    subdf["nuclei_bbox_global(px)"] = nuclei_bbox
    subdf["nuclei_type"] = nuclei_type
    subdf["nuclei_type_probability"] = nuclei_type_probability
    subdf["number_of_nuclei_in_subpatch"] = nuclei_number
    return subdf


def start_hovernet(stored_patches, save_loc, hovernet_command, master_loc, img_name, coords_x, coords_y, data_saver_mode):
    """ hovernet is run in sets of patches, because running in wsi mode was buggy for mrxs images"""
    
    call = hovernet_command[4]
    a = hovernet_command[0] # which gpu to use
    b = hovernet_command[1] # number of types of nuclei to predict 
    c = hovernet_command[2] # model mode
    d  = hovernet_command[3] # model path
    
    for i in stored_patches:
        pic = i[0]
        key = i[1]
        input_dir, output_dir = make_folders(save_loc=save_loc, key=key)
        
        #save images as tifs for hovernet to run
        Image.fromarray(pic).save(f"{save_loc}/{key}/original_patches/{key}_{coords_x}_{coords_y}.tif")
        
        ## preparing command to run hovernet
        if int(b) > 0:
            type_info_path={output_dir}
            command = ["python3", "run_infer.py", f"--gpu={a}",f"--nr_types={int(b)}", f"--type_info_path=type_info.json", f"--model_mode={c}",f"--model_path={d}", 'tile', f"--input_dir={input_dir}", f"--output_dir={output_dir}"]    
        else:
            command = ["python3", "run_infer.py",f"--gpu={a}", f"--model_path={d}",f"--model_mode={c}", 'tile', f"--input_dir={input_dir}", f"--output_dir={output_dir}"]
        

        ## start the process
        os.chdir(master_loc) # this works and gets you to the right directory
        proc = subprocess.Popen(command, env=dict(os.environ, KMP_DUPLICATE_LIB_OK = "TRUE"), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait(timeout=None) #wait for patches to finish processing
        #print(proc.stderr.readlines()) #debug
        print("                 |_entering hovernet process // pid:", proc.pid)

        ## data cleanup
        delete_all_tifs_from_directory(f"{save_loc}/{key}/original_patches/") #delete b/c easy to regenerate
        if data_saver_mode:
            delete_mats(f"{save_loc}/{key}/generated_output/mat/")
            delete_pngs(f"{save_loc}/{key}/generated_output/overlay/")
        


def make_folders(save_loc, key):
    input = f"{save_loc}{key}/original_patches/"
    output = f"{save_loc}{key}/generated_output/"
    Path(input).mkdir(parents=True, exist_ok=True) #where to save the patches
    Path(output).mkdir(parents=True, exist_ok=True) #where hovernet should put its output
    return input, output

def delete_all_tifs_from_directory(path):
    """delete saved patches that result from having to run hovernet in patch mode"""
    files = os.listdir(path)
    for file in files:
        if file.endswith(".tif"):
            os.remove(f"{path}{file}")

def delete_mats(path):
    """delete saved patches that result from having to run hovernet in patch mode"""

    files = os.listdir(path)
    for file in files:
        if file.endswith(".mat"):
            os.remove(f"{path}{file}")

def delete_pngs(path):
    """delete saved patches that result from having to run hovernet in patch mode"""

    files = os.listdir(path)
    for file in files:
        if file.endswith(".png"):
            os.remove(f"{path}{file}")

def put_coords_in_global_context(lst_of_coords, global_x, global_y):
    return [lst_of_coords[0]+float(global_x), lst_of_coords[1]+ float(global_y)]
    
def put_lst_lst_in_global_context(lst_of_lst_coords, global_x, global_y):
    r = []
    for i in lst_of_lst_coords:
        r.append([i[0] + float(global_x), i[1] + float(global_y)])
    return r

def calculate_area(contour):
    area = type(contour) #NoneType if not a list
    if type(contour)== list:
        nuc = np.array(contour)
        hull = ConvexHull(nuc)
        area = hull.volume #returns an area since contour is 2D, see documentation
    return area

def crawl_through_data(path):
    """crawls through folder strucsture as defined in run_algo to gather all the hovernet data into one df. path starts at temproray_data"""
    df_lst = []
    for patch in os.listdir(f"{path}"):
        patch_path = f"{path}{patch}/generated_output/"
        df = consolidate_hovernet_data(patch_path)
        df_lst.append(df)
    return pd.concat(df_lst)


def consolidate_hovernet_data(path):
    """put hovernet output into dataframe with global context"""
    df_name = []
    df_patch_key = []
    df_inst_centroid = []
    df_inst_bbox = []
    df_inst_contour = []
    df_inst_type = []
    df_inst_type_prob =[]
    df_area = []
    df_local_coords = []

    json_path = path + "json/"
    img_name = json_path.split("/")[-6]
    assert len(os.listdir(json_path)) <= 1 #sanity check. There should only ever be one patch per folder
    if len(os.listdir(json_path)) == 1:
        json_file = os.listdir(json_path)[0]
        global_key = json_file.split("_")[0]
        coords_x = json_file.split("_")[-2]#global with respect to wsi, in px
        coords_y = json_file.split("_")[-1].replace(".json","")#global with respect to wsi, in px
        with open(json_path + json_file) as file:
            data = json.loads(file.read())
            nuc_info = data['nuc']
            """
            Hovernet data:
                mag_info = data['mag']
                nuc_info = data['nuc']
                nuc info has following information:
                    'bbox','centroid', 'contour','type_prob','type'
            """
            for inst in nuc_info:
                inst_info = nuc_info[inst]
                inst_centroid = inst_info['centroid']
                inst_contour = inst_info['contour']
                inst_bbox = inst_info['bbox']
                type_inst = inst_info["type"] 
                type_prob = inst_info["type_prob"]
                
                #print("check if global_centeroid calcualation is correct!!! seems weird to me")

                global_centroid = put_coords_in_global_context(lst_of_coords=inst_centroid, global_x=coords_x, global_y=coords_y)
                global_contour = put_lst_lst_in_global_context(lst_of_lst_coords=inst_contour, global_x=coords_x, global_y=coords_y)
                global_bbox = put_lst_lst_in_global_context(lst_of_lst_coords=inst_bbox, global_x=coords_x, global_y=coords_y)
                area_px = calculate_area(global_contour)
                
                df_local_coords.append(inst_centroid)
                df_name.append(img_name)
                df_patch_key.append(global_key)
                df_inst_centroid.append(global_centroid)
                df_inst_bbox.append(global_bbox)
                df_inst_contour.append(global_contour)
                df_inst_type.append(type_inst)
                df_inst_type_prob.append(type_prob)
                df_area.append(area_px)

    df_data = {"Name": df_name, "patch_key": df_patch_key, "local_centroid(px)": df_local_coords, "centroid_global(px)": df_inst_centroid, "bbox_global(px)": df_inst_bbox, "contour_global(px)": df_inst_contour, "type": df_inst_type, "type_probability": df_inst_type_prob, "area_nuc(pxs)": df_area}
    df = pd.DataFrame(data=df_data)
    return df

def load_preprocessed_img(loc, wanted_key):
    """loads an image by patch key"""
    for i in os.listdir(loc):
        patch_key = i.split("_")[0] #check this :D
        coords_x = i.split("_")[1]
        coords_y = i.split("_")[2].replace(".tif","")
        if patch_key == str(wanted_key):
            img = np.array(Image.open(loc + i))
            return img, int(coords_x), int(coords_y)
        
def patch_in_directory(wanted_key, loc):
    """originally hovernet was run on all patches, even if there was no tissue there. Hovernet recognizes nuclei in unscanned area. to fix this, hovernet is only run on patches with tissue detected in new run. however when using data from old hovernet run, we have to skip patches with no tissue (patches which werent saved and therefore can't be loaded)"""
    for patch in os.listdir(loc):
        patch_key = patch.split("_")[0]
        if patch_key == wanted_key:
            return True
    return False
            

def globalize_hovernet_data(raw_df, binary_loc):
    """to avoid noise and detection of nuclei in 'dirty' areas, only count tissue located in a tissue area"""
    dfs_lst = []
    key_lst = np.unique(raw_df["patch_key"].to_list()) 
    for key in key_lst:
        nucs_in_patch = raw_df[raw_df["patch_key"]==key].copy() #split nuclei by which binary image they are located in
        
        patch_exists = patch_in_directory(key, binary_loc)
    
        if patch_exists:
            preprocessed_img, coords_x, coords_y = load_preprocessed_img(loc=binary_loc, wanted_key = key) #load the appropriate binary image
            to_keep = []
            count = 0 
            for index, row in nucs_in_patch.iterrows():
                coords = row["centroid_global(px)"]
                local_x = coords[0] - coords_x
                local_y = coords[1] - coords_y
                #assert (coords[0] >= coords_x) & (coords[0] <= coords_x + 10000)
                #assert (coords[1] >= coords_y) & (coords[0] <= coords_y + 10000)
                count += 1
                if preprocessed_img[int(local_y), int(local_x)] == 0: #careful, loads [y:x] ## THEN is located in a tissue area
                    to_keep.append(True)
                else:
                    to_keep.append(False)
        else:#patch wasn't saved, since no tissue detected
            print("patch not saved as binary_patch: ", key)
            to_keep = []
            for index, row in nucs_in_patch.iterrows():
                to_keep.append(False)
        nucs_in_patch["is_a_relevant_nucleus"] = to_keep
        dfs_lst.append(nucs_in_patch)
    final_df = pd.concat(dfs_lst) # to access only "real" nuclei, do: df[df["is_a_relevant_nucleus"]==True]
    return final_df


def visualize_nuclei(df, img_loc):
    """visualize the nuclei in a dataframe on top of the original patch."""
    keys = np.unique(df["patch_key"].to_list())
    for key in keys:
        # if key == some_key: ## to visualize only a certain patch
        img, x_coords, y_coords = load_preprocessed_img(loc=img_loc, wanted_key=key) #load patch & return global patch coords
        nucs_in_key = df[df["patch_key"]==key]
        relevant_nuclei = nucs_in_key[nucs_in_key["is_a_relevant_nucleus"] == True]
        irrelevant_nuclei = nucs_in_key[nucs_in_key["is_a_relevant_nucleus"] == False]
       
        centroids = relevant_nuclei["local_centroid(px)"].to_list()
        xs_rel = [i[0] for i in centroids]
        ys_rel = [i[1] for i in centroids]
        
        centroids = irrelevant_nuclei["local_centroid(px)"].to_list()
        xs_irrel = [i[0] for i in centroids]
        ys_irrel = [i[1] for i in centroids]
        
        plt.imshow(img, cmap="binary")
        plt.plot(xs_rel, ys_rel, marker="v", color="green", linestyle='None')
        plt.plot(xs_irrel, ys_irrel, marker="v", color = "red", linestyle='None')
        plt.show()