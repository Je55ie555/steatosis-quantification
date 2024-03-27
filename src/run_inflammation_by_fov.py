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
    help="path to fat and nucleus detection output",
    type=str,
)
parser.add_argument(
    "--CCA-folder-name", required=True, help="path to CCA_dfs", type=str
)
parser.add_argument(
    "--FOV_output",
    required=True,
    help="name for the output folder.",
    type=str,
)
args = parser.parse_args()

path = args.data_dir
CCA_experiment_name = args.CCA_folder_name
FOV_experiment_name = args.FOV_output


run_Infl_by_foci_for_all_imgs(
    path=path,
    CCA_experiment_name=CCA_experiment_name,
    FOV_experiment_name=FOV_experiment_name,
)
