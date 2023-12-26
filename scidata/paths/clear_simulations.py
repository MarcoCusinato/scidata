import argparse
from path_managment import clear_simulation_folder
parser = argparse.ArgumentParser()
parser.add_argument('--simulations', nargs='+', type=str, 
                    help="Expected one or multiple simulations name in string format.")


args = parser.parse_args()
clear_simulation_folder(args.simulations)