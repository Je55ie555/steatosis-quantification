import argparse
from datetime import datetime
import os
from os import path
from pathlib import Path
from CCA_nuclei import *

parser = argparse.ArgumentParser(
    prog="Nucleus Cluster Analysis",
    description="This code takes information about nuclei location, and analyzes the nuclei distributions.",
    epilog="See ReadMe for more information",
)

parser.add_argument(
    "--data-dir",
    required=True,
    help="path to fat and nucleus detection output ([file_name]_fat_and_nucs.pkl). The output of nucleus cluster analysis will also be stored under this location",
    type=str,
)

parser.add_argument(
    "--output-folder-name",
    required=True,
    help="name for the output folder. Date & time are automatically added to prevent overwriting",
    type=str,
)
parser.add_argument(
    "--max-distance",
    required=True,
    help="How far apart nuclei can be located apart and still be considered neighbors. In micrometers. Default: 15 µm.",
    type=int,
    default=15,
)

args = parser.parse_args()

df_column_label = f"area(px)_taken_up_by_nucleus_clusters_max_dis_{args.max_distance}(micrometers)"  # this column is added to the [file_name]_fat_and_nucs.pkl dataframe

run_CCA(
    base_path=args.data_dir,
    max_distance=args.max_distance,
    experiment_name=args.output_folder_name,
    removed_df_label=df_column_label,
)
