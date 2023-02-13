import argparse
from path_managment import list_of_paths, add_paths
parser = argparse.ArgumentParser()
parser.add_argument('--add-paths', nargs='+', type=str, 
                    help="Expected one or multiple paths in string format.")
parser.add_argument('--add-paths-file', nargs='+', type=str,
                    help="Expected one or multiple paths to file in string format.")

args = parser.parse_args()
path_list = list_of_paths(args.add_paths_file, args.add_paths)
add_paths(path_list)