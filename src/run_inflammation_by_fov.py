from inflammation_by_fov import *
import argparse

parser = argparse.ArgumentParser(
    prog="Process nuclei by Field of View",
    description="This script allows an analysis of the number of nuclei clusters in a given (square-shaped) field of view (FOV).",
    epilog="See ReadMe for more information",
)

parser.add_argument(
    "--data-dir",
    required=True,
    help="path to fat and nucleus detection output ([file_name]_fat_and_nucs.pkl). The output of nucleus cluster analysis will also be stored under this location",
    type=str,
)
parser.add_argument(
    "--output-dir", required=True, help="path to output directory", type=str
)
parser.add_argument(
    "--output-folder-name",
    required=True,
    help="name for the output folder. Date & time are automatically added to prevent overwriting",
    type=str,
)
args = parser.parse_args()

path = args.data_dir
experiment_name = args.output_folder_name
save_loc = args.output_dir
processed = os.listdir(
    f"{save_loc}/{experiment_name}/num_foki_per_20fov/"
)  # prevent images from being rerun :D. Remember to change hardcoded part if using another fov.
processed = [
    i.replace("_foki_per_20xfov.pkl", "") for i in processed
]  # creates a list of the names of files that have already been processed.

run_Infl_by_foci_for_all_imgs(
    path=path, processed=processed, experiment_name=experiment_name, save_loc=save_loc
)
