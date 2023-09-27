from scidata.quantities.quantities import SimulationAnalysis
import argparse
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--sim-name', required=True, type=str, help="Simulation name")
parser.add_argument('--sim-path', required=False, default=None, help="Optional: path of the simulation")
parser.add_argument('--sure', action='store_true', default=False, help="Bypass the input from user")
args = parser.parse_args()

sim = SimulationAnalysis(simulation_name=args.sim_name, simulation_folder_path=args.sim_path)

if args.sure:
    print("Removing files in " + sim.storage_path + "\nPlease wait...")
    #shutil.rmtree(sim.storage_path)
    #os.mkdir(sim.storage_path)
    print("Done.")
else:
    sure = input("\nThis will delete all the post-processing file in " + sim.storage_path + \
                '. Are you sure to continue? (Y: continue, other: stop)\n')
    if sure == 'Y':
        print("Removing files in " + sim.storage_path + "\nPlease wait...")
        shutil.rmtree(sim.storage_path)
        os.mkdir(sim.storage_path)
        print("Done.")
