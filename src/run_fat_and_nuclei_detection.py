import argparse
from datetime import datetime
import os
from os import path
from pathlib import Path
from fat_detection_script import *

parser = argparse.ArgumentParser(
                    prog='Fat and Nuclei Detection',
                    description='This code processes a WSI in patches in order to detect fatty vesicles. Optionally, if analyzing a H&E WSI, nuclei can be detected.',
                    epilog='See ReadMe for more information')

parser.add_argument('--input-dir', required=True, help="path to WSI image data directory") #file_loc

## output organization
parser.add_argument('--output-folder-name', required=True, help="name for the output folder. Date & time are automatically added to prevent overwriting", type=str)
parser.add_argument('--output-dir', required=True, help="output directory", type=str) #run path

## WSI information
parser.add_argument('--WSI-file-type', required=True, help="can be '.mrxs' or '.czi'", type=str)

## determining which WSIs to process
parser.add_argument('--WSI-queue', required=False, help="A list of comma separated file names as a string, ie 'file_a,file_b' [no spaces!] of the WSIs to be processed in the input directory. When not specified, all files with the given file_type in the directory are processed.", type=str, default="") #imgs_to_process

## HoVer-Net parameters (for nucleus detection)
parser.add_argument( '--run-hovernet', required=False, help="must be 'True' to run HoVer-Net", action="store_true") 
parser.add_argument( '--hovernet-dir', required=False, help="path to HoVer-Net code")
parser.add_argument( '--hovernet-weights-file', required=False, help="path to weights (.tar file)", type=str)
parser.add_argument( '--hovernet-gpu', required=False, help="which gpu number to use run HoVer-Net", type=int)
parser.add_argument( '--hovernet-num-nuc-types', required=False, help="how many nuclei types are to be detected. Depends on weights being used.", type=int)
parser.add_argument( '--hovernet-model-mode', required=False, help="can either be 'fast' or 'original'", type=str)
parser.add_argument( '--hovernet-run-function', required=False, help="see HoVer-Net documentation. Default: './run_infer.py'", type=str, default='./run_infer.py') #call
parser.add_argument( '--data-saver-mode', required=False, help="handle HoVer-Net temporary data", action = "store_true")

## Processind parameters
parser.add_argument( '--patch-size', required=False, help="Patch size in micrometers. Patches are always processed in squares. Default is 1500 micrometers", type=int, default=1500) 
parser.add_argument( '--patch-level', required=False, help="level at which patches are processed. Default is 0.", type=int, default=0) 
parser.add_argument( '--min-fat-vesicle-area', required=False, help="minimum macrovesicular fat vesicle area (in micrometers squared). Default is 40 micrometers squared.", type=int, default=40) 
parser.add_argument( '--max-fat-vesicle-area', required=False, help="maximum macrovesicular fat vesicle area (in micrometers squared). Default is 10000 micrometers squared.", type=int, default=10000) 
parser.add_argument( '--min-extent', required=False, help="minimum extent of a fat vesicle. Default is 0.5", type=float, default=0.5) 
parser.add_argument( '--max-axis-ratio', required=False, help=" maximum axis ratio of a fat vesivle. Default: 2", type=int, default=2) 
parser.add_argument( '--subpatch-size-factor', required=False, help=" int representing how many subpatches a patch will be divided into. Default=20.", type=int, default=20) 
parser.add_argument( '--disk', required=False, help="Kernel for image processing. Default=8.3478513356562", type=float, default=8.3478513356562)

args = parser.parse_args()


#check inputs and create output dir with correct naming
run_time = datetime.now().strftime("%Y%m%d-%H%M%S")
assert path.exists(args.input_dir) 
output_path = args.output_dir + args.output_folder_name + run_time # run_path
Path(f"{output_path}").mkdir(parents=True, exist_ok=False)#create the output directory

#check if hovernet parameters are complete
if args.run_hovernet:
    hovernet_args = [args.hovernet_dir,
                     args.hovernet_gpu,
                     args.hovernet_num_nuc_types,
                     args.hovernet_weights_file,
                     args.hovernet_model_mode,
                     args.hovernet_run_function,
                     args.data_saver_mode
    ]

    for arg in hovernet_args:
        if arg is None:
            print("if running HoVer-net, you must provide all the HoVer-net parameters:")
            print("     --hovernet-dir")
            print("     --hovernet-weights_file")
            print("     --hovernet-gpu")
            print("     --num-nuc-types")
            print("     --hovernet-model-mode")
            print("     --hovernet-run-function")
            assert 1==2
    else:
        hovernet_command = [args.hovernet_gpu, args.hovernet_num_nuc_types, args.hovernet_model_mode, args.hovernet_weights_file, args.hovernet_run_function]
if not(args.run_hovernet):
    hovernet_command = []

    

#calculate remaining parameters
FILL_FG_SIZE = args.min_fat_vesicle_area - 0.2 #in micrometers
FILL_BG_SIZE = args.max_fat_vesicle_area + 0.2 #in micrometers



#converts WSI-queue input (str) into a list of strings (file names to process)
if len(args.WSI_queue) == 0: #empty list
    imgs_to_process = []
else:
    imgs_to_process = args.WSI_queue.split(',')


run_workflow(
    file_location = args.input_dir,
    run_path = args.output_dir,
    PATCH_SIZE = args.patch_size,
    patch_level = args.patch_level,
    MIN_FAT_AREA = args.min_fat_vesicle_area,
    MAX_FAT_AREA = args.max_fat_vesicle_area,
    min_extent = args.min_extent,
    max_axis_ratio = args.max_axis_ratio,
    FILLED_BG_HOLE_SIZE=FILL_BG_SIZE,
    FILLED_FG_HOLE_SIZE=FILL_FG_SIZE,
    file_type = args.WSI_file_type,
    KERNEL= args.disk,
    SUBPATCH_SIZES_FACTOR = args.subpatch_size_factor,
    hovernet_master_loc = args.hovernet_dir,
    run_hovernet = args.run_hovernet,
    hovernet_command = hovernet_command,
    data_saver_mode = args.data_saver_mode,
    imgs_to_process = imgs_to_process
)