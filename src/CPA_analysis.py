from fat_detection_script import *
import cv2
from skimage import morphology
import pandas as pd
from skimage.color import label2rgb
from sklearn.cluster import KMeans
from PIL import Image, ImageFilter
from skimage.segmentation import mark_boundaries
from math import cos, sin, pi
import skimage
from skimage import filters
from skimage.morphology import disk

def distance_to_blue(pts):
    blue = np.array([0,0,255])
    dis = [math.sqrt(( blue[0] -i[0] )**2 + ( blue[1] -i[1])**2 + (  blue[2]- i[2])**2) for i in pts]
    return dis

def distance_to_white(pts):
    white = np.array([255,255,255])
    dis = [math.sqrt(( white[0] -i[0] )**2 + (white[1] -i[1])**2 + (  white[2]- i[2])**2) for i in pts]
    return dis



def find_avg_num_pix_where_blue_is_max(patch, mask):
 
    
    patch_with_pixels_selected = patch[mask==False]
    blues = patch_with_pixels_selected[:,2]
    reds = patch_with_pixels_selected[:,0] #high blue and red channels -> purple --> only pixels where blue is the "strongest" should be selected
    greens = patch_with_pixels_selected[:,1]
    template = np.zeros_like(blues)
    template[(blues > reds) & (blues > greens)] = 1
    
    mean_pixels_with_max_blue = np.mean(template) #mean of blue channel
    return mean_pixels_with_max_blue


def convert_to_cartesian_space(lst_of_degrees):
    """transform HSV colorspace to non-cyclical rep.
     Based on  https://stackoverflow.com/questions/35402609/point-on-circle-base-on-given-angle/35402676
     Quentin Pradet, Feb 15, 2016"""
    center = [0,0]
    radius = 1
    hue_radians = [math.radians(i *2) for i in lst_of_degrees]
    xs = [center[0]+ (radius * cos(angle)) for angle in hue_radians]
    ys = [center[1]+ (radius * sin(angle)) for angle in hue_radians]
    return xs,ys
    
                            
    
    
    
    
def preprocess_masson(mask, patch, filled_bg_hole_size, filled_fg_hole_size, kernel):
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
    return closed_patch, preprocessed_patch





def detect_CPA(
    file_location,
    run_path,
    PATCH_SIZE,
    patch_level, 
    MIN_FAT_AREA,
    MAX_FAT_AREA,
    file_type,
    min_extent,
    max_axis_ratio,
    FILLED_BG_HOLE_SIZE,
    FILLED_FG_HOLE_SIZE,
    KERNEL,
    imgs_to_process
):
    print("LOGS 'run_workflow': \n")
    counter = 0
    print(file_location)
    dfl = []
    for file in os.listdir(file_location):


        if (any([file in i for i in imgs_to_process]) & (file.endswith(file_type))) or ((len(imgs_to_process)==0) & (file.endswith(file_type))):
            
            
            print(f"- Processing file: {file}")
            
            """open images. Currently supported datatypes mrxs and czi"""
            print(f"----------opening wsi to determining global threshold")
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
            

            Path(f"{output_loc}").mkdir(parents=True, exist_ok=True)
            Path(f"{fat_loc}").mkdir(parents=True, exist_ok=True)
            Path(f"{binary_loc}").mkdir(parents=True, exist_ok=True)
            Path(f"{output_loc}/patches/").mkdir(parents=True, exist_ok=True)

            
            csv = create_csv(f=f, csv_files_path=output_loc)
        
            """adjust parameters to mpp"""
            print(f"----------adjusting default parameters to wsi scan")
            
            min_fat_area, max_fat_area, filled_fg_hole_size, filled_bg_hole_size, kernel, mag_factor, patch_size = convert_and_adjust(mpp, MIN_FAT_AREA, MAX_FAT_AREA, FILLED_BG_HOLE_SIZE, FILLED_FG_HOLE_SIZE, KERNEL, PATCH_SIZE)
            
            area_threshold_of_tissue_to_remove = (PATCH_SIZE * PATCH_SIZE *.05)/(mpp*2) #tissue object must take up at least 5% of area of patch
            
            print(f"----------calculating global threshold")
            
            
            #find regions with tissue based on saturation value
            scanned_sections = return_scanned_parts(wsi_img_array, img_type)
            wsi_hsv = cv2.cvtColor(wsi_img_array, cv2.COLOR_RGB2HSV)
            wsi_hsv_scanned_parts = wsi_hsv[scanned_sections ==True]
            global_threshold = skimage.filters.threshold_otsu(wsi_hsv_scanned_parts[:,1])
            
            

            key = 0
            key_l = []
            blue_pixels = []
            blue_pixels_without_edges = []
            tissue_without_fat_pixels =[]
            fat_pixels = []
            total_tissue_pixels = []
            mpps = []
            coordsxl = []
            coordsyl = []
            name = []
            stain_vector_channel_selected = []

            
            subpatch_analysis = []
            print(f"----------total number of patches in wsi (at patch size {patch_size} (px)): { math.ceil(float(scanned_width / patch_size)) * math.ceil(float(scanned_height / patch_size))}")
            for coords_y in range(start_y, finish_y, patch_size):
                for coords_x in range(start_x, finish_x, patch_size):
                    print("             |_ ", key)
                    
                    if (key>=0) : #keep in case a particular patch key is of interest :)
                        
                        patch = open_patch(img_type, wsi, wsi_img_array,coords_x , coords_y, patch_size, patch_level)
                        
                      
                        
                        scanned_patch = return_scanned_parts(img=patch, img_type=img_type)
                        
                        
                        #convert patch into hsv and apply global threshold to make binary
                        img_hsv = cv2.cvtColor(patch, cv2.COLOR_RGB2HSV)
                        tissue_mask = img_hsv[:,:,1] > global_threshold
                        #remove small pieces of tissue that could be dirt (any isolated object that is less than 5% of patch area)
                        tissue_mask_cleaned = np.invert(morphology.remove_small_holes(np.invert(tissue_mask), area_threshold = area_threshold_of_tissue_to_remove))
                        
                        liver_tissue_without_white_areas = np.sum(tissue_mask_cleaned)
                        
                        percent_patch_with_liver_tissue = (liver_tissue_without_white_areas/(patch_size**2)) * 100
                        
                       
                        if percent_patch_with_liver_tissue > 1: #if more than 1% tissue was found
                          
                            #Find Fat
                            closed_patch, preprocessed_patch = preprocess_masson(
                                mask = tissue_mask_cleaned,
                                patch=patch,
                                filled_bg_hole_size=filled_bg_hole_size,
                                filled_fg_hole_size=filled_fg_hole_size,
                                kernel=kernel
                            )
                          
                            patch_with_fat_objects, forground_area = object_analysis( 
                            closed_patch=closed_patch,
                            name=file,
                            min_area=min_fat_area,
                            max_area=max_fat_area,
                            extent=min_extent,
                            max_axis_ratio=max_axis_ratio,
                            x=coords_x,
                            y=coords_y,
                            glob_th=global_threshold,
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
                            print("finished fat analysis")
                           
                            #kmeans cluster idea
                            patch_no_alpha= patch[:,:,0:3]

                            #patch_no_alpha = skimage.filters.gaussian(patch_no_alpha)
                            original_shape = patch_no_alpha.shape


                          
                           
                            #blue range is 180-300. Have to be divided by two b/c of csv hue scale (0-180)
                            xs_blue_ideal, ys_blue_ideal= convert_to_cartesian_space(lst_of_degrees=[(-240)])
                            
                            #red is located around 0 or 360
                            xs_red_ideal, ys_red_ideal = convert_to_cartesian_space(lst_of_degrees=[(360)])
                            
                            
                            # create initial centroids for kmeans
                            blue_array = np.array([xs_blue_ideal[0], ys_blue_ideal[0]])
                            red_array = np.array([xs_red_ideal[0], ys_red_ideal[0]])
                            init = np.stack([blue_array, red_array])
                       
                          
                            
                            #prepare pixels for fitting kmeans
                            hue_degrees = img_hsv[tissue_mask_cleaned==True]### cluster only pixels in detected tissue
                            xs_detected_tissue, ys_detected_tissue = convert_to_cartesian_space(lst_of_degrees=hue_degrees[:,0])
                            transformed_pixels_detected_tissue = np.column_stack([xs_detected_tissue, ys_detected_tissue])
                            
                            
                            
                            
                            
                            
                            #fit kmenas
                            kmeans = KMeans(n_clusters=2, random_state=42, init=init, n_init=1)
                            labels_original = kmeans.fit_predict(transformed_pixels_detected_tissue)
                            
                           
                            
                            #predict on whole patch
                            #print("transforming HSV color space")
                            hsv_image = img_hsv
                            hsv_shape = img_hsv.shape
                            hue_degrees_whole_patch = hsv_image[:,:,0].reshape(hsv_shape[0]* hsv_shape[1],1)
                            
                            xs_whole_patch, ys_whole_patch = convert_to_cartesian_space(lst_of_degrees=hue_degrees_whole_patch)
                            transformed_pixels_whole_patch = np.column_stack([xs_whole_patch, ys_whole_patch])
                            
                            labels = kmeans.predict(transformed_pixels_whole_patch)
                            labels = np.array([i + 1 for i in labels]) # we dont want a zero label :D for masking purposes
                            clustered_image = labels.reshape(hsv_shape[0], hsv_shape[1])#mask of the patch with labels
                            centers = kmeans.cluster_centers_
                            
                            
                            mask1 = np.logical_and((clustered_image == 1), tissue_mask_cleaned)
                            mask2 = np.logical_and((clustered_image == 2), tissue_mask_cleaned)
                            
                            #find average hue for each mask
                            hue_mask1_degrees = [np.mean(img_hsv[mask1==True][:,0])]
                            hue_mask2_degrees = [np.mean(img_hsv[mask2==True][:,0])]
                        
                            
                            xs_m1, ys_m1 = convert_to_cartesian_space(lst_of_degrees=hue_mask1_degrees)
                            xs_m2, ys_m2 = convert_to_cartesian_space(lst_of_degrees=hue_mask2_degrees)
                            
                            
                     
                            #for each mask find the distance to blue and red
                            mask1_distance_to_blue = math.dist((xs_m1[0], ys_m1[0]), (xs_blue_ideal[0], ys_blue_ideal[0]))
                            mask1_distance_to_red = math.dist((xs_m1[0], ys_m1[0]), (xs_red_ideal[0], ys_red_ideal[0]))

                            mask2_distance_to_blue = math.dist((xs_m2[0], ys_m2[0]), (xs_blue_ideal[0], ys_blue_ideal[0]))
                            mask2_distance_to_red = math.dist((xs_m2[0], ys_m2[0]), (xs_red_ideal[0], ys_red_ideal[0]))

                            # print("mask1", mask1_distance_to_blue, mask1_distance_to_red)
                            # print("mask2", mask2_distance_to_blue, mask2_distance_to_red)
                            
                            
                            

                            red_distances = [mask1_distance_to_red, mask2_distance_to_red]
                            blue_distances = [mask1_distance_to_blue, mask2_distance_to_blue]

                            selected_red = np.argmin(red_distances) +1 #b/c not zero index
                            selected_blue = np.argmin(blue_distances) +1 #b/c not zero index
                            #print("selected red", selected_red)
                            #print("selected blue", selected_blue)
                            
                            if selected_red == selected_blue:
                                """
                                fig,ax=plt.subplots(ncols = 4, figsize=(40,10))
                                ax[0].imshow(patch_no_alpha)
                                ax[1].imshow(mask1)
                                ax[2].imshow(mask2)
                                ax[3].imshow(patch_no_alpha)
                                ax[3].imshow(selected_mask, cmap="gray", alpha=0.5)
                                plt.show()      
                                """
                                mask_originally_selected = selected_red
                                print("SELECTED SAME CLUSTER FOR RED AND BLUE WHAT DO")
                                #this happens when there are very few points in a particular cluster, for example.
                                #for that mask, check if in all it's closer to blue or to red
                                #if it is closer, then make that the blue mask
                                # if not, make that the red mask
                                
                                if selected_red == 1:
                                    distances = [mask1_distance_to_blue, mask1_distance_to_red]
                                if selected_red == 2:
                                    distances = [mask2_distance_to_blue, mask2_distance_to_red]
                                
                                min_dist = np.argmin(distances) #which color is it closest to?
                                if min_dist == 0:
                                    color_selected = "blue"
                                    selected_blue = mask_originally_selected
                                    #red has to be assigned the other mask
                                    if mask_originally_selected == 1:
                                        selected_red = 2 
                                    if mask_originally_selected == 2:
                                        selected_red = 1
                                    
                                if min_dist ==1:
                                    color_selected = "red"
                                    selected_red = mask_originally_selected
                                    #blue has to be assigned the other mask
                                    if mask_originally_selected == 1:
                                        selected_blue = 2 
                                    if mask_originally_selected == 2:
                                        selected_blue = 1
                                        
                                        
                            #print("finally selected red", selected_red)
                            #print("finally selected blue", selected_blue)
                            
                            
                            selected_mask = np.logical_and((clustered_image == selected_blue), tissue_mask_cleaned)
                            red_mask = np.logical_and((clustered_image == selected_red), tissue_mask_cleaned)
                       
                            


                           
                            
                            num_blue_pixels = np.sum(selected_mask)
                            tissue_pixels_no_fat = np.sum(tissue_mask_cleaned)
                            fat_tissues = np.sum(np.invert(patch_with_fat_objects))
                            
                            cpa_percent = (num_blue_pixels/(tissue_pixels_no_fat + fat_tissues)) *100
                            
                            

                            name.append(file)
                            key_l.append(key)
                            coordsxl.append(coords_x)
                            coordsyl.append(coords_y)
                            blue_pixels.append(num_blue_pixels)
                            total_tissue_pixels.append(tissue_pixels_no_fat + fat_tissues)
                            tissue_without_fat_pixels.append(tissue_pixels_no_fat)
                            fat_pixels.append(fat_tissues)
                            mpps.append(mpp)
                            
                            
                            """
                            #### PLOTTING
                            fig,ax=plt.subplots(ncols = 5, figsize=(40,10))
                            ax[0].imshow(patch_no_alpha)
                            ax[1].imshow(mask1)
                            ax[1].set_title(np.sum(mask1))
                            ax[2].imshow(mask2)
                            ax[2].set_title(np.sum(mask2))
                            ax[3].set_title(f"finally selected blue {selected_blue}, % cpa = {cpa_percent}, {np.sum(selected_mask)}")
                            ax[3].imshow(patch_no_alpha)
                            ax[3].imshow(selected_mask, cmap="gray", alpha=0.5)
                            ax[4].imshow(patch_with_fat_objects)
                            plt.savefig(f"{run_path}/{f}/flow_chart.svg")
                            plt.show() 
                            
                            
                            #plot the distribution of color in hue cast to a cartesian space
                            colors = []
                            for i in labels_original:
                                if i == 1:
                                    colors.append("C3") #red
                                if i == 0:
                                    colors.append("C0") #blue
                            plt.figure(figsize=(5,5))
                            plt.gca().set_aspect(1.0)
                            plt.scatter(xs_detected_tissue, ys_detected_tissue, color=colors)
                            plt.scatter(xs_blue_ideal, ys_blue_ideal, "X",  label="Blue", color="blue")
                            plt.scatter(xs_red_ideal, ys_red_ideal, "X",  label="Red", color="red")
                            plt.savefig(f"{run_path}/{f}/point_distribution.svg")
                            plt.show()
                            """

                       
                            
                            key +=1
                        else:
                            key +=1



                    
                        
                    else:
                        
                        key +=1     
                       
             
            if  (img_type=="mrxs"):           
                wsi.close()
            csv.close()
            
            
            data = {"name": name, "original_patch_key": key_l, "global_x_coords(pxs)": coordsxl, "global_y_coords(pxs)": coordsyl, "mpps": mpps,  "num_blue_pixels":blue_pixels, "total_tissue_pixels(with_fat)": total_tissue_pixels,"num_pixels_tissue_without_fat":tissue_without_fat_pixels, "num_fat_pixels": fat_pixels}
            
            
            df = pd.DataFrame(data)
            df.to_pickle(f"{output_loc}{file.replace(file_type,'')}_CPA.pkl")
        else:
            continue
        



