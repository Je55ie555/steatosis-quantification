import os
import numpy as np
import openslide
import skimage.filters
from skimage.color import rgb2gray
from skimage import measure
import skimage.segmentation
from skimage.feature import peak_local_max
import pandas as pd
import csv
from findmaxima2d import find_maxima # find_local_maxima
from scipy import ndimage as ndi
import math
import cv2 as cv
from PIL import Image
from pathlib import Path
from hovernet_functions import *
import czifile


################# Data Management ####################################

def create_csv(f, csv_files_path):
    """creates a csv file("_data.csv") with information about every selected object.
    for one patch following information stays constant:
        patch_key
        patch_x
        patch_y
        black_area
    for wsi stays constant:
        WSI
        global_threshold
        patch_level
        th_mpp
        mag"""
    data_loc = csv_files_path + f + "_data.csv"
    csv_file = open(data_loc, "w", newline="")
    csv_file.write( "WSI,global_threshold,patch_level,th_mpp,patch_key,patch_x,patch_y,Centroid_x,Centroid_y,Area,Bbox_x,Bbox_y,Bbox_width,Bbox_height,Extent,Axis_ratio,is_fat,Black_area,Circularity,Mpp,Magnification\n"
    ) 
    return csv_file


def write_line_in_csv(
    csv,
    img,
    th,
    img_level,
    th_mpp,
    k,
    x,
    y,
    c_x,
    c_y,
    area,
    bbox_x,
    bbox_y,
    bbox_width,
    bbox_height,
    ext,
    is_fat,
    tissue,
    circ,
    axis_ratio,
    mpp_x,
    mag
):
    
    csv.write(
        "%s, %.5f, %i, %.5f, %i, %i, %i, %.3f, %.3f, %.5f, %.3f, %.3f, %.3f, %.3f, %.5f, %.5f, %i, %.5f, %.5f, %.5f, %.5f\n"
        % (
            img,
            th,
            img_level,
            th_mpp,
            k,
            x,
            y,
            c_x,
            c_y,
            area,
            bbox_x,
            bbox_y,
            bbox_width,
            bbox_height,
            ext,
            axis_ratio,
            is_fat,
            tissue,
            circ,
            mpp_x,
            mag
        )
    )





###################### Opening Images, Metadata Extraction, Parameter Adjustments ####################

def find_target_level(levels, wsi, target):
    mpp = [] # assume x and y mpps are same
    
    for i in range(int(levels)):
        prop = "mirax.LAYER_0_LEVEL_" + str(i) + "_SECTION.MICROMETER_PER_PIXEL_X"
        mppx= float(wsi.properties[prop])
        mpp.append(mppx)
    
    abs_dif = lambda value: abs(value - target)
    closest_value = min(mpp, key=abs_dif)
    
    for i in range(int(levels)):
        prop = "mirax.LAYER_0_LEVEL_" + str(i) + "_SECTION.MICROMETER_PER_PIXEL_X"
        mppx= float(wsi.properties[prop])
        if mppx == closest_value:
            level = i
        else:
            continue

    return level


def open_mrxs(file_location, file, target):
    """takes in file location, returns the opened wsi (as a np.array) at level for calculating global threshold as defined by target. 
    Also returns  image parameters"""

    wsi = openslide.OpenSlide(os.path.join(file_location, file))
    mag = int(wsi.properties['openslide.objective-power'])
    scanned_height = int(wsi.properties["openslide.bounds-height"])
    scanned_width = int(wsi.properties["openslide.bounds-width"])
    start_x = int(wsi.properties["openslide.bounds-x"])
    start_y = int(wsi.properties["openslide.bounds-y"])
    finish_y = start_y + scanned_height
    finish_x = start_x + scanned_width
    levels = wsi.properties[ 'openslide.level-count']
    mpp_x = float(wsi.properties["mirax.LAYER_0_LEVEL_0_SECTION.MICROMETER_PER_PIXEL_X"])
    mpp_y = float(wsi.properties["mirax.LAYER_0_LEVEL_0_SECTION.MICROMETER_PER_PIXEL_Y"])
    assert mpp_x == mpp_y # after this, code assumes x and y scales are the same
    level = find_target_level(levels, wsi, target)
    wsi_img_array = np.array(
        wsi.read_region(
            (0, 0),
            level,
            (wsi.level_dimensions[level][0], wsi.level_dimensions[level][1]),
        ))
    return wsi_img_array, mag, scanned_height, scanned_width, start_x, start_y, finish_x, finish_y, mpp_x, wsi

def open_czi(file_location, file, target):
    """takes in file location, returns the opened wsi (as a np.array) at level for calculating global threshold as defined by target. 
    Also returns image parameters"""

    czi = czifile.CziFile(file_location+file)
    image = czifile.imread(file_location+file)
    img_array = np.zeros_like(image[0,0,0,:,:,:])
    for i in range(image.shape[1]):#for all the scenes
        scene_img = image[0, i, 0 :,:,:]
        scene = scene_img[0,:,:,:]
        img_array [scene> 0 ] = scene[scene > 0] ## stiched image at highest resolution  
    md_string = czi.metadata() #raw=False creates dictionary, none creates string
    czi_parsed = etree.fromstring(md_string)
    mag = int(czi_parsed.xpath("//NominalMagnification")[0].text)
    md = czi.metadata(raw=False)
    scaling_info = md["ImageDocument"]["Metadata"]["Scaling"]["Items"]["Distance"]
    mpps_x = scaling_info[0]["Value"]
    mpps_y = scaling_info[1]["Value"]
    
    assert mpps_x == mpps_y # after this, code assumes x and y scales are the same
    
    mpp = mpps_x *(10**6)  #mpps are in m => convert to micrometers to match mrxs units
    scanned_height = img_array.shape[1]
    scanned_width = img_array.shape[0]
    start_x = 0
    start_y = 0
    finish_y = start_y + scanned_height
    finish_x = start_x + scanned_width
    factor = target // mpp ## calculates approximately same mpp
    wsi_img_array = np.array(Image.fromarray(img_array).resize((int(scanned_height//factor), int(scanned_width//factor))))
    return wsi_img_array, mag, scanned_height, scanned_width, start_x, start_y, finish_x, finish_y, mpp


def open_wsi(file_location, file, target):
    if file.endswith(".mrxs"):
        f=file.replace(".mrxs","") 
        wsi_img_array, mag, scanned_height, scanned_width, start_x, start_y, finish_x, finish_y, mpp, wsi = open_mrxs(file_location, file, target)
        img_type = "mrxs"
        
        
    if file.endswith(".czi"):
        f=file.replace(".czi","")
        wsi_img_array, mag, scanned_height, scanned_width, start_x, start_y, finish_x, finish_y, mpp = open_czi(file_location, file, target)
        img_type = "czi"
        wsi = np.nan
    return wsi_img_array, mag, scanned_height, scanned_width, start_x, start_y, finish_x, finish_y, mpp, f, img_type, wsi

def convert_from_um_to_mpp(is_area, mpp, um):
    if is_area:
        pixels = um/((mpp*mpp))
    else:
        pixels = um/mpp
    return pixels


def convert_and_adjust(mpp, min_fat_area, max_fat_area, filled_bg_hole_size, filled_fg_hole_size, kernel, patch_size):
    """ adjusts default parameters to magnification/scan size of the wsi, based on pixel size.
     Default parameters were calculated at mpps:013913085559427"""

    #convert_from_um_to_mpp (is_area, mpp, um) for each of these parameters (because input is given in micrometers)
    min_fat_area = convert_from_um_to_mpp(is_area=True, mpp=mpp, um=min_fat_area)
    max_fat_area = convert_from_um_to_mpp(is_area=True, mpp=mpp, um=max_fat_area)
    filled_bg_hole_size = convert_from_um_to_mpp(is_area=True, mpp=mpp, um=filled_bg_hole_size)
    filled_fg_hole_size = convert_from_um_to_mpp(is_area=True, mpp=mpp, um=filled_fg_hole_size)
    kernel = convert_from_um_to_mpp(is_area=False, mpp=mpp, um=kernel)
    patch_size = int(convert_from_um_to_mpp (is_area=False, mpp=mpp, um=patch_size))
    
    

    #adjusts default parameters to magnification
    orig_mpp = 0.13913085559427 ## original parameters were calculated at these mpps.
    if mpp < orig_mpp:
        print("READUJST MPP DISTANCE MEASUREMENTS, does not support <0.139 mpps at the moment") #should theoretically work, hasn't been tested yet. Saftey break :D
        assert 1==2
    else: 
        factor = mpp / orig_mpp #float
        int_factor = int(np.round(mpp/orig_mpp)) ## not ideal..... .5> gets rounded up
        area_adjustment = factor **2
        min_fat_area = min_fat_area/area_adjustment
        max_fat_area = max_fat_area/area_adjustment
        filled_fg_hole_size = filled_fg_hole_size/area_adjustment
        filled_bg_hole_size = filled_bg_hole_size/area_adjustment
        kernel = int(kernel/int_factor)
        mag_factor = int_factor
    return min_fat_area, max_fat_area, filled_fg_hole_size, filled_bg_hole_size, kernel, mag_factor, patch_size

def return_scanned_parts(img, img_type): 
    """ returns subsection of img that contains information as an array"""

    if img_type == "czi":
        m1 = img[:, :, 0] != 0
        m2 = img[:, :, 1] != 0
        m3 = img[:, :, 2] != 0
        scanned = (m1 & m2) & m3
    if (img_type == "mrxs"):
        scanned= img[:,:, 3] != 0 
    return scanned

def open_patch(img_type, wsi, wsi_img_array, coords_x, coords_y, patch_size, patch_level):
    """opens a patch of the wsi, returns patch as an array"""

    if (img_type=="mrxs"):
        patch = np.array(
            wsi.read_region(
                (coords_x, coords_y), patch_level, (patch_size, patch_size)
            )
        ) #THE ARRAY IS STORED AS y,x       
    if img_type ==".czi":
        patch = wsi_img_array[coords_x:coords_x+patch_size, coords_y:coords_y+patch_size]
    return patch


################## Image Functions ########################
def global_threshold(wsi, img_type):
    """Uses otsu method to calculate threshold for scanned parts of the wsi, using pixel variance"""

    scanned_sections = return_scanned_parts(wsi, img_type)
    threshold = skimage.filters.threshold_otsu(wsi[scanned_sections ==True][:, 0:3].std(axis=1)) #threshold with variance information
    return threshold

def apply_global_threshold(threshold, patch):
    """applies global threshold on a patchwise basis"""

    variance = patch.std(axis=2)
    mask = np.zeros_like(variance,  dtype=bool)
    mask[variance > threshold] = 1
    return mask

def remove_holes_from_patch(patch, fg_hole_size, bg_hole_size):
    """using input to remove objects too small or too large from the fg and bg (as defined by presets)"""

    patch_filled_bg = skimage.morphology.remove_small_holes(
        np.invert(patch), area_threshold=bg_hole_size
    )
    patch_filled_fg = skimage.morphology.remove_small_holes(
        np.invert(patch_filled_bg), area_threshold=fg_hole_size
    )
    return patch_filled_fg.astype(np.uint8)

def apply_closing(struct_element, patch):
    """uses a rectangular element to perform a closing on an img"""
    kernel = cv.getStructuringElement(cv.MORPH_RECT ,(struct_element, struct_element)) 
    closed_patch = cv.morphologyEx(patch, cv.MORPH_CLOSE, kernel)
    return closed_patch

def preprocess (patch, global_threshold, filled_bg_hole_size, filled_fg_hole_size, kernel):
    """
        1. apply global threshold of pixel variance to patch ## we found that using pixel variance information instead of intensity, or intensity and varaince improved the results
        2. remove holes that are too large or too small to be fat objects from patch (uses preset variables)
        3. preform a closing (uses preset variables)
        4. perform a hysteresis thresholding using gray patch, closed patch and mask of the variance information
        5. remove small holes remaining after the hysteresis
        # note, closed patch is returned in addition to the preprocessed patch as it is needed in the watershedding step
    """

    mask = apply_global_threshold(global_threshold, patch)
    filled_patch = remove_holes_from_patch(mask, filled_fg_hole_size, filled_bg_hole_size)
    closed_patch = apply_closing(kernel, filled_patch)
    gray_patch = rgb2gray(patch[:, :, 0:3])
    hysteresis_patch = skimage.filters.apply_hysteresis_threshold(
        image=gray_patch, low=mask, high=closed_patch
    )
    preprocessed_patch = skimage.morphology.remove_small_holes(
        hysteresis_patch, area_threshold=filled_bg_hole_size
    )
    closed_patch = np.invert(closed_patch.astype(bool))
    return closed_patch, preprocessed_patch, gray_patch


########## object functions ###########
def object_analysis(
    closed_patch,# used for watershedding
    name,
    min_area,
    max_area,
    max_axis_ratio,
    extent,
    x, #of patch, in pixels
    y, #of patch, in pixels
    glob_th, #global threshold
    key,
    csv_file,
    scanned, #parts of scanned img
    img_level,#for csv data
    th_mpp, #for csv data
    mag,
    preprocessed_patch,#used to be called th_patch
    mpp, #for csv data
    magnification
):
    """calculates black area (tissue area without fat area), performs watershed, and then processes objects. Writes data to CSV file"""

    black_area = np.sum(preprocessed_patch[scanned == 1] == 0)
    regions, seg = watershed(closed_patch=closed_patch, max_area=max_area, mag=mag, preprocessed_patch=preprocessed_patch)
    tissue_mask = np.ones_like(closed_patch)
    # write in csv #
    if len(regions) == 0:  # no objects found in image  
        write_line_in_csv(
                csv=csv_file,
                img=name,
                th=glob_th,
                img_level=img_level,
                th_mpp=th_mpp,
                k=key,
                x=x,
                y=y,
                c_x=np.nan,
                c_y=np.nan,
                area=np.nan,
                bbox_x=np.nan,
                bbox_y=np.nan,
                bbox_width=np.nan,
                bbox_height=np.nan,
                ext=np.nan,
                is_fat=0,
                tissue=black_area,
                circ=np.nan,
                axis_ratio=np.nan,
                mpp_x=mpp,
                mag=magnification
            )
        
        #total_tissue = np.invert(tissue_mask & closed_patch)
        return tissue_mask, black_area #, total_tissue, black_area 
    else:
        for region in regions:
            is_fat = filters(
                region=region,
                min_area=min_area,
                max_area=max_area,
                extent=extent,
                max_axis_ratio=max_axis_ratio,
            )
            r_area = region["Area"]
            r_bbox_x = region["bbox"][0]
            r_bbox_y = region["bbox"][1]
            r_bbox_width = region["bbox"][2]
            r_bbox_height = region["bbox"][3]
            r_ext = region["extent"]
            r_circ = (4 * math.pi * region["Area"]) / ((region["Perimeter"]) ** 2)
            r_centroid_x = region["centroid"][1]
            r_centroid_y = region["centroid"][0]
            if (region["major_axis_length"] ==0) or (region["minor_axis_length"] == 0 ):
                axis_ratio = np.nan
            else:
                axis_ratio = region["major_axis_length"] / region["minor_axis_length"]
            write_line_in_csv(
                csv=csv_file,
                img=name,
                th=glob_th,
                img_level=img_level,
                th_mpp=th_mpp,
                k=key,
                x=x,
                y=y,
                c_x=r_centroid_x,
                c_y=r_centroid_y,
                area=r_area,
                bbox_x=r_bbox_x,
                bbox_y=r_bbox_y,
                bbox_width=r_bbox_width,
                bbox_height=r_bbox_height,
                ext=r_ext,
                is_fat=is_fat,
                tissue=black_area,
                circ=r_circ,
                axis_ratio=axis_ratio,
                mpp_x=mpp,
                mag=magnification
            )
            if not(is_fat):
                continue
            if is_fat:
                tissue_mask[seg ==region["label"]] = 0 
        #total_tissue = np.invert(tissue_mask & closed_patch)
        #total_tissue[scanned==0]=0 #removes unscanned regions
        return tissue_mask, black_area #, total_tissue, black_area
    

def seed_watershed(patch, max_area, magnification):
    """uses the closed image to place the seeds for the watershed. """

    filled_patch = fill_large_holes(patch=patch, max_fat_area=max_area)  #makes sure no seeds are placed in background. Only objs that are potentially fat remain
    distance_map = skimage.filters.gaussian(ndi.distance_transform_edt(filled_patch), sigma=5//magnification) 
    peak_max_coords = peak_local_max(distance_map, min_distance=50//magnification, threshold_abs=8//magnification) # returns coordinates of local peaks in img
    peak_mask = np.zeros(distance_map.shape, dtype=bool) #create mask
    peak_mask[tuple(peak_max_coords.T)] = True #fill in local maxes
    h, w, _ = find_maxima(distance_map, peak_mask, 10//magnification) #find the maxima of the local maxes. 
    local_max_mask = np.zeros(distance_map.shape, dtype=bool)
    for m_x, m_y in zip(
        w, h  
    ): 
        local_max_mask[m_y, m_x] = True # set local_max_mask to true at the coordinates of the maxima
    markers = measure.label(local_max_mask)  # seeds for watershed area non-zero values
    return markers

def watershed(closed_patch, max_area, mag, preprocessed_patch):
    """watershed to separate connected objects into two separate objects"""

    seeds = seed_watershed(patch=closed_patch, max_area=max_area, magnification=mag)
    filled_patch = fill_large_holes(patch = preprocessed_patch, max_fat_area = max_area) # Only objs that are potentially fat remain         
    distance_map = skimage.filters.gaussian(ndi.distance_transform_edt(filled_patch), sigma=5//mag)
    seg = skimage.segmentation.watershed(-distance_map, seeds, mask=filled_patch, watershed_line=False) #only points where mask == True are labeled.
    regions = skimage.measure.regionprops(seg)
    return regions, seg

def fill_large_holes(patch, max_fat_area):
    p = np.copy(patch)  #prevents overwriting
    label_img = skimage.measure.label(p)
    regions = skimage.measure.regionprops(label_img)
    for region in regions:
        if region["Area"] > 5 * max_fat_area: 
            label = region["label"]
            mask = label_img == label
            p[mask] = 0
    return p

def filters(region, min_area, max_area, extent, max_axis_ratio):
    """uses features to determine if object is likely fat.
    returns a boolean"""

    is_fat = True
    if region["Area"] < min_area or region["Area"] > max_area:
        is_fat = False
    if region["Extent"] < extent:
        is_fat = False
    if (region["major_axis_length"] ==0) or (region["minor_axis_length"] == 0 ):
        print(f"problem fat spot:{region['centroid'][1]},{region['centroid'][0]}")
        is_fat = False
    else:
        if (region["major_axis_length"] / region["minor_axis_length"]) > max_axis_ratio:
            is_fat = False
    return is_fat

def save_patches(preprocessed_patch, patch_with_fat_objs, binary_loc, fat_loc, key, coords_x, coords_y):
    """saves patches into a given location"""
    preprocessed_patch_img = Image.fromarray(preprocessed_patch)
    fat_obj_img =Image.fromarray(patch_with_fat_objs)        
    preprocessed_patch_img.save(f"{binary_loc}{key}_{coords_x}_{coords_y}.tif")
    fat_obj_img.save(f"{fat_loc}{key}_{coords_x}_{coords_y}.tif")         


def process_subpatches(preprocessed_patch, size, patch_size, subpatch_analysis, scanned, patch_with_fat, coords_x, coords_y, key, mpp):
    """this way, image can be processed in a smaller resolution without compromising object detection due to overlap with borders of patch"""

    subpatch = []
    subpatch_key = 0
    total_foreground = 0
    total_fat = 0 
    
    for y in range(0, patch_size, size): 
        for x in range(0, patch_size, size):
            b_patch =  preprocessed_patch[y:y+size, x:x + size]
            fat_patch =  patch_with_fat[y:y+size, x:x + size]
            scanned_subpatch = scanned[y:y+size, x:x + size]
            foreground_area = np.sum(b_patch[scanned_subpatch==1]==0) # I don't need to do the scanned==0 part, could delete. 
            fat_area = np.count_nonzero(fat_patch==0)
            subpatch.append([size, key, subpatch_key, x + coords_x, y + coords_y, fat_area, foreground_area, mpp])
            subpatch_key +=1 
            total_foreground += foreground_area
            total_fat += fat_area
    subpatch_analysis.append([subpatch])
    return subpatch_analysis, total_foreground, total_fat 


def save_subpatch_info_into_df(subpatch_analysis, output_loc, name_for_file, patch_size):
    """adds information from subpatch proceessing to the whole dataframe"""
    subpatch_size = []
    subpatch_key = []
    original_key = []
    global_x_coords = []
    global_y_coords = []
    fa = []
    ba = []
    mpps = []
    for patch in subpatch_analysis:
        for sp_list in patch:
            for subpatch in sp_list:
                subpatch_size.append(subpatch[0])
                original_key.append(subpatch[1])
                subpatch_key.append(subpatch[2])
                global_x_coords.append(subpatch[3])
                global_y_coords.append(subpatch[4])
                fa.append(subpatch[5])
                ba.append(subpatch[6])
                mpps.append(subpatch[7])
    data = {"original_key": original_key, "subpatch_size(px)": subpatch_size,"subpatch_key":subpatch_key, "global_x_coords(px)":global_x_coords, "global_y_coords(px)": global_y_coords, "fat_area(px)": fa, "black_area(px)":ba, "mpps":mpps }
    df = pd.DataFrame(data)
    df["original_patch_size(px)"] = [patch_size] * len(df)
    df.to_pickle(f"{output_loc}{name_for_file.replace('.','')}_subpatch_df.pkl")
    return df


##### MAIN ################

def run_workflow(
    file_location,
    run_path,
    PATCH_SIZE,
    patch_level, 
    MIN_FAT_AREA,
    MAX_FAT_AREA,
    min_extent,
    max_axis_ratio,
    FILLED_BG_HOLE_SIZE,
    FILLED_FG_HOLE_SIZE,
    file_type,
    KERNEL,
    SUBPATCH_SIZES_FACTOR,
    hovernet_master_loc,
    run_hovernet,
    hovernet_command,
    data_saver_mode,
    imgs_to_process
):
    
    print("LOGS 'run_workflow': \n")
    counter = 0
    print(file_location)
    for file in os.listdir(file_location):
    
        if (any([file in i for i in imgs_to_process]) & (file.endswith(file_type))) or ((len(imgs_to_process)==0) & (file.endswith(file_type))):
            
            print(f"- Processing file: {file}")
            
           

            """open images. Currently supported datatypes mrxs and czi """
            print(f"---------- opening wsi to determining global threshold")
            target= 4.45 # target mpp, corresponds to dev pic level 5 mpp. Used for development
            wsi_img_array, mag, scanned_height, scanned_width, start_x, start_y, finish_x, finish_y, mpp, f, img_type, wsi = open_wsi (file_location, file, target)
            
            """Data managment"""
            print("---------- data directory: ")
            output_loc = f"{run_path}/{f}/dataframes/"
            binary_loc = f"{run_path}/{f}/saved_patches/binary_patches/"
            fat_loc = f"{run_path}/{f}/saved_patches/fat_patches/"
            
           
            print(f"             |_ location dataframe & csv storage :{output_loc}")
            print(f"             |_ location binary patch storage :{binary_loc}")
            print(f"             |_ location fat patch storage :{fat_loc}")
            

            if run_hovernet:
                hovernet_loc = f"{run_path}/{f}/temporary_data/"
                print(f"             |_ location hovernet storage (if hovernet enabled):{hovernet_loc}")
        
            Path(f"{output_loc}").mkdir(parents=True, exist_ok=True)
            Path(f"{fat_loc}").mkdir(parents=True, exist_ok=True)
            Path(f"{binary_loc}").mkdir(parents=True, exist_ok=True)
           
            csv = create_csv(f=f, csv_files_path=output_loc)
        
            """adjust parameters to mpp"""
            print(f"---------- adjusting default parameters to wsi scan")

            min_fat_area, max_fat_area, filled_fg_hole_size, filled_bg_hole_size, kernel, mag_factor, patch_size = convert_and_adjust(mpp, MIN_FAT_AREA, MAX_FAT_AREA, FILLED_BG_HOLE_SIZE, FILLED_FG_HOLE_SIZE, KERNEL, PATCH_SIZE)
            
            
            ## subpatch_size must evenly divide into patch_size! 
            print(f"------- original patch_size: {patch_size}")
            patch_size = int(math.ceil(patch_size/SUBPATCH_SIZES_FACTOR) * SUBPATCH_SIZES_FACTOR)
            print(f"-------- patch_size adjusted for even fitting of subpatches: {patch_size}")
            print(f"------- number of subpatches in a patch: {SUBPATCH_SIZES_FACTOR}")
            subpatch_sizes = [int(patch_size/SUBPATCH_SIZES_FACTOR)]
            
            print("----------  mapping patch size to a suitable value for subpatch processing")
            print(f"---------- calculating global threshold")
            g_threshold = global_threshold(wsi_img_array, img_type)
            
            if run_hovernet:
                stored_patches = []
                patches_stored = 0

            key = 0
            subpatch_analysis = []
            print(f"---------- total number of patches in wsi (at patch size {patch_size} (px)): { math.ceil(float(scanned_width / patch_size)) * math.ceil(float(scanned_height / patch_size))}")
            print("---------- processing patch number:")
            for coords_y in range(start_y, finish_y, patch_size):
                for coords_x in range(start_x, finish_x, patch_size):
                    #print("             |_ ", key)
                    if (key >= 0):#keep in case a particular patch key is of interest :)

                        ### if you want to process only one patch at a given location....
                        # coords_x = start_x + int(6474/0.139)
                        # coords_y = start_y + int(3170/0.139)
                        # print("FIXED X AND Y COORDS!")
                        # patch_to_save = Image.fromarray(patch)
                        # patch_to_save.save(f"{fat_loc}{key}_{coords_x}_{coords_y}_sanity_test.tif")
                        # print("done saving imgae")
                        # print(f"{fat_loc}{key}_{coords_x}_{coords_y}_sanity_test.tif")
                        

                        patch = open_patch(img_type, wsi, wsi_img_array, coords_x, coords_y, patch_size, patch_level) 

                        scanned_patch = return_scanned_parts(img=patch, img_type=img_type)
                        closed_patch, preprocessed_patch, gray_patch = (
                            preprocess(
                                patch=patch,
                                global_threshold = g_threshold,
                                filled_bg_hole_size=filled_bg_hole_size,
                                filled_fg_hole_size=filled_fg_hole_size,
                                kernel=kernel
                            ))
                      
                        
                        patch_with_fat_objects, forground_area = object_analysis( 
                            closed_patch=closed_patch,
                            name=file,
                            min_area=min_fat_area,
                            max_area=max_fat_area,
                            extent=min_extent,
                            max_axis_ratio=max_axis_ratio,
                            x=coords_x,
                            y=coords_y,
                            glob_th=g_threshold,
                            key=key,
                            csv_file=csv,
                            scanned=scanned_patch,
                            img_level=patch_level,
                            th_mpp=target,
                            mag= mag_factor,
                            preprocessed_patch=preprocessed_patch,
                            mpp=mpp,
                            magnification=mag
                        )
                        

                        
                        preprocessed_patch[scanned_patch==0] = 1 # fill in unscanned regions as background
                        
                        if forground_area > 0: #only patches with tissue on it will be processed further
                            
                            save_patches(preprocessed_patch=preprocessed_patch, patch_with_fat_objs=patch_with_fat_objects, binary_loc=binary_loc, fat_loc=fat_loc, key=key, coords_x=coords_x, coords_y=coords_y)#patches need to be saved for connected component anlysis (fat overlap) and for nuclei anaysis (check if nuclei are located in an area with tissue) 
                            
                            if run_hovernet:
                                if patches_stored % 10 == 0: #change to mod 10 after testing. 
                                    stored_patches.append([patch, key])#stored_patches.append([Image.fromarray(patch), key])
                                    start_hovernet(stored_patches=stored_patches, save_loc=hovernet_loc, hovernet_command=hovernet_command, master_loc=hovernet_master_loc, img_name=f, coords_x=coords_x, coords_y=coords_y, data_saver_mode=data_saver_mode)
                                    stored_patches = []
                                    patches_stored = 0
                                    
                                else:
                                    stored_patches.append([patch, key]) 
                                    patches_stored += 1

                            if len(subpatch_sizes) > 0:
                                
                                total_foreground_area = 0
                                total_fat_area = 0
                               
                                

                                for size in subpatch_sizes:
                                    
                                    subpatch_analysis, ba, fa = process_subpatches(preprocessed_patch=preprocessed_patch, size=size, patch_size=patch_size, subpatch_analysis=subpatch_analysis, scanned=scanned_patch, patch_with_fat=patch_with_fat_objects, coords_x=coords_x, coords_y=coords_y, key=key, mpp=mpp)
                                    total_foreground_area += ba
                                    total_fat_area += fa
                                assert (total_foreground_area == forground_area) & (total_fat_area == np.count_nonzero(patch_with_fat_objects==0))
                                subpatch_df = save_subpatch_info_into_df(subpatch_analysis, output_loc, file.replace(img_type,""), patch_size)

                        key +=1

                    else:
                        key +=1
             
            if  (img_type=="mrxs"):           
                wsi.close()
            csv.close()

            if run_hovernet:
                print("running hovernet")
                if patches_stored > 0: ## run hovernet for the last of the pictures in the stored_patch_que
                    start_hovernet(stored_patches=stored_patches, save_loc=hovernet_loc, hovernet_command=hovernet_command, master_loc=hovernet_master_loc, img_name = f, coords_x=coords_x, coords_y=coords_y, data_saver_mode=data_saver_mode)

                print("---------- processing hovernet data")
                
               
                consolidated_hovernet_df = crawl_through_data(hovernet_loc) # CAVE contains nuclei in unscanned regions!
                assert len(consolidated_hovernet_df) > 0 # makes sure that nuclei were detected
                raw_hovernet_df = globalize_hovernet_data(consolidated_hovernet_df, binary_loc)
                relevant_nuclei_df = raw_hovernet_df[raw_hovernet_df["is_a_relevant_nucleus"] == True]
                raw_hovernet_df.to_pickle(f"{output_loc}{f}_raw_global.pkl")
                relevant_nuclei_df.to_pickle(f"{output_loc}{f}_relevant_nuclei.pkl")
                
                #after this in theory json data can be deleted. consider saving since it shouldn't take too much memory and is only needed 1/wsi
                
                algo_hvnet_df = combine_dataframes (subdf=subpatch_df, nucdf=relevant_nuclei_df) #lists all the nuclei (with info) in each subpatch                
                algo_hvnet_df["patch_size(px)"] = [patch_size] * len(algo_hvnet_df)
                algo_hvnet_df["subpatch_size(px)"] = [subpatch_sizes[0]] * len(algo_hvnet_df)
                algo_hvnet_df.to_pickle(f"{output_loc}{f}_fat_and_nucs.pkl")
        
                if data_saver_mode: 
                    delete_all_tifs_from_directory(binary_loc) 
                    delete_all_tifs_from_directory(fat_loc)
            
            print(f"----------Completed file: {file}----------")
            
        else:
            continue

