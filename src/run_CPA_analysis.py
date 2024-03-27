import argparse
from datetime import datetime
import os
from os import path
from pathlib import Path
from CPA_analysis import *


parser = argparse.ArgumentParser(
    prog="CPA Analysis",
    description="This script analyzes the collagen proportionate area (CPA) in Masson's Trichrome stained WSI. ",
    epilog="See ReadMe for more information",
)

parser.add_argument(
    "-i", "--input-dir", required=True, help="path to WSI image data directory"
)  # file_loc

## output organization
parser.add_argument(
    "--output-folder-name",
    required=True,
    help="name for the output folder. Date & time are automatically added to prevent overwriting",
    type=str,
)
parser.add_argument(
    "--output-dir", required=True, help="output directory", type=str
)  # run path

## WSI information
parser.add_argument(
    "--WSI-file-type", required=True, help="can be '.mrxs' or '.czi'", type=str
)

## determining which WSIs to process
parser.add_argument(
    "--WSI-queue",
    required=False,
    help="A list of comma separated file names as a string, ie 'file_a,file_b' [no spaces!] of the WSIs to be processed in the input directory. When not specified, all files with the given file_type in the directory are processed.",
    type=str,
    default="",
)  # imgs_to_process


## Processind parameters
parser.add_argument(
    "--patch-size",
    required=False,
    help="Patch size in micrometers. Patches are always processed in squares. Default is 1500 micrometers",
    type=int,
    default=1500,
)
parser.add_argument(
    "--patch-level",
    required=False,
    help="level at which patches are processed. Default is 0.",
    type=int,
    default=0,
)
parser.add_argument(
    "--min-fat-vesicle-area",
    required=False,
    help="minimum macrovesicular fat vesicle area (in micrometers squared). Default is 40 micrometers squared.",
    type=int,
    default=40,
)
parser.add_argument(
    "--max-fat-vesicle-area",
    required=False,
    help="maximum macrovesicular fat vesicle area (in micrometers squared). Default is 10000 micrometers squared.",
    type=int,
    default=10000,
)
parser.add_argument(
    "--min-extent",
    required=False,
    help="minimum extent of a fat vesicle. Default is 0.5",
    type=float,
    default=0.5,
)
parser.add_argument(
    "--max-axis-ratio",
    required=False,
    help="maximum axis ratio of a fat vesivle. Default: 2",
    type=int,
    default=2,
)
parser.add_argument(
    "--disk",
    required=False,
    help="Kernel for image processing. Default=8.3478513356562",
    type=float,
    default=8.3478513356562,
)

args = parser.parse_args()


# converts WSI-queue input (str) into a list of strings (file names to process)
if len(args.WSI_queue) == 0:  # empty list
    imgs_to_process = []
else:
    imgs_to_process = args.WSI_queue.split(",")

# calculate remaining parameters
FILL_FG_SIZE = args.min_fat_vesicle_area - 0.2  # in micrometers
FILL_BG_SIZE = args.max_fat_vesicle_area + 0.2  # in micrometers


detect_CPA(
    file_location=args.input_dir,
    run_path=args.output_dir,
    PATCH_SIZE=args.patch_size,
    patch_level=args.patch_level,
    MIN_FAT_AREA=args.min_fat_vesicle_area,
    MAX_FAT_AREA=args.max_fat_vesicle_area,
    min_extent=args.min_extent,
    max_axis_ratio=args.max_axis_ratio,
    FILLED_BG_HOLE_SIZE=FILL_BG_SIZE,
    FILLED_FG_HOLE_SIZE=FILL_FG_SIZE,
    file_type=args.WSI_file_type,
    KERNEL=args.disk,
    imgs_to_process=imgs_to_process,
)
