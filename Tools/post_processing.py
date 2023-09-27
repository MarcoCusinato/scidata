import scivis
import os
from scidata.quantities.quantities import SimulationAnalysis
import sys
import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("--sim-name", required=True, help="Name of the simulation to analyse")
parser.add_argument("--sim-path", required=False, default=None, help="Path of the selected simulation")
parser.add_argument("--rxlim-GWs", required=False, default=None, type=float, help="Right x-axis limit in ms, if None full GW is plotted")
parser.add_argument("--overwrite-plots", action='store_true', default=False, help='Remakes all the plots even if they already exists.')
parser.add_argument("--GW-movie", action='store_true', help='Generate a movie for the full')
parser.add_argument("--GW-contribution", action='store_true', help='Generate a plot with the radial contribution for the GWs')
parser.add_argument("--remove-post-processing", action='store_true', default=False, help="Removes post processing done.")
args = parser.parse_args()

# Disable print options
def disablePrint():
    sys.stdout = open(os.devnull, 'w')

# Restore print command
def enablePrint():
    sys.stdout = sys.__stdout__

# Disable warnings
warnings.filterwarnings("ignore")

## Get the path of the plotting an movieing (? is this even a word?) scripts
plots_scripts_path = os.path.join(os.path.dirname(scivis.__file__), 'plots')
movies_scripts_path = os.path.join(os.path.dirname(scivis.__file__), 'movies')

## Get the simulation

sim = SimulationAnalysis(args.sim_name, None, args.sim_path)
args.sim_path = os.path.dirname(sim.path)

plot_folder = os.path.expanduser('~/Postprocessing_plots')

if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)
plot_folder = os.path.join(plot_folder, 
                           args.sim_path.split('/')[-1])
if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)
plot_folder = os.path.join(plot_folder,
                           args.sim_name)
if not os.path.exists(plot_folder):
    os.mkdir(plot_folder)

if args.remove_post_processing:
    os.system('python remove_post_processing.py --sim-name ' + args.sim_name + \
              ' --sim-path ' + args.sim_path + ' --sure')

try:
    print('Calculating Brunt-Väisälä frequency...')
    disablePrint()
    sim.BV_frequency_profile()
    enablePrint()
    print('Brunt-Väisälä frequency calculated.')
except:
    enablePrint()
    print('Error in calculating Brunt-Väisälä frequency.')

try:
    print('Calculating GW amplitude from density and velocities...')
    disablePrint()
    sim.AE220()
    enablePrint()
    print('GW amplitude from density and velocities calculated.')
except:
    enablePrint()
    print('Error in calculating GW amplitude from density and velocities.')

try:
    print('Calculating shock radius...')
    disablePrint()
    sim.get_shock_radius()
    enablePrint()
    print('Shock radius calculated.')
except:
    enablePrint()
    print('Error in calculating shock radius.')

try:
    print('Calculating PNS radius...')
    disablePrint()
    sim.get_PNS_radius()
    enablePrint()
    print('PNS radius calculated.')
except:
    enablePrint()
    print('Error in calculating PNS radius.')

try:
    print('Calculating gain radius...')
    disablePrint()
    sim.get_gain_radius()
    enablePrint()
    print('Gain radius calculated.')
except:
    enablePrint()
    print('Error in calculating gain radius.')

try:
    print('Calculating dipole modes...')
    disablePrint()
    sim.get_modes(1)
    enablePrint()
    print('Dipole modes calculated.')
except:
    enablePrint()
    print('Error in calculating dipole modes.')

try:
    print('Calculating masses and energies...')
    disablePrint()
    sim.get_masses_and_energies()
    enablePrint()
    print('Masses and energies calculated.')
except:
    enablePrint()
    print('Error in calculating masses and energies.')

try:
    print('Calculating innercore mass and energies...')
    disablePrint()
    sim.get_innercore_mass_enregies()
    enablePrint()
    print('Innercore mass and energy calculated.')
except:
    enablePrint()
    print('Error in calculating innercore mass and energies.')

try:
    print('Calculating neutrino sphere radii...')
    disablePrint()
    sim.get_neutrino_sphere_radius()
    enablePrint()
    print('Neutrino sphere radii calculated.')
except:
    enablePrint()
    print('Error in calculating neutrino sphere radii.')

try:
    print('Calculating innercore radius...')
    disablePrint()
    sim.get_innercore_radius()
    enablePrint()
    print('Innercore radius calculated.')
except:
    enablePrint()
    print('Error in calculating innercore radius.')

print('Post processing files created.')

# Plot post processing
try:
    print('Plotting entropy profile.')
    name = 'Entropy_profile.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Entropy.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name + ' --vmin 0 --vmax 22')
    print('Entropy profile plotted.')
except:
    print('Error in plotting the entropy profile.')

try:
    print('Plotting pressure profile.')
    name = 'Pressure_profile.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Pressure.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name + ' --vmin 1e18 --vmax 1e34')
    print('Entropy profile plotted.')
except:
    print('Error in plotting the pressure profile.')

try:
    print('Plotting Ye profile.')
    name = 'Ye_profile.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Ye.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name + ' --vmin 0 --vmax 0.5')
    print('Ye profile plotted.')
except:
    print('Error in plotting the Ye profile.')

try:
    print('Plotting temperature profile.')
    name = 'Temperature_profile.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Temperature.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name + ' --vmin 0 --vmax 30')
    print('temperature profile plotted.')
except:
    print('Error in plotting the temperature profile.')

try:
    print('Plotting explosion energy and mass.')
    name = 'E-M_expl.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Explosions.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('Explosion mass and energy plotted.')
except:
    print('Error in plotting the explosion mass and energy.')

try:
    print('Plotting neutrino luminosities and energies.')
    name = 'L-E_neutrino.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Neutrinos.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('Neutrino luminosities and energies plotted.')
except:
    print('Error in plotting the neutrino luminosities and energies.')

try:
    print('Plotting maximum density.')
    name = 'Rho_max.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Rho_max.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('Maximum density plotted.')
except:
    print('Error in plotting the maximum density.')

try:
    print('Plotting shock radius.')
    name = 'Shock_radius.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Shock_radius.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('Shock radius plotted.')
except:
    print('Error in plotting the shock radius.')

try:
    print('Plotting gain radius.')
    name = 'Gain_quantities.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Gain.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('Gain radius plotted.')
except:
    print('Error in plotting the gain radius.')

try:
    print('Plotting neutrino sphere radii.')
    name = 'Nutrino_sphere_radii.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'Neutrino_spheres.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('Neutrino sphere radii plotted.')
except:
    print('Error in plotting the neutrino sphere radii.')

try:
    print('Plotting PNS radius.')
    name = 'PNS_quantities.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'PNS.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('PNS radius plotted.')
except:
    print('Error in plotting the PNS radius.')

try:
    print('Plotting GWs amplitude.')
    name = 'GWs_strain.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'GWs_amplitudes.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name + ' --rxlim ' + str(args.rxlim_GWs))
    print('GWs amplitude plotted.')
except:
    print('Error in plotting the GWs amplitude.')

try:
    print('Plotting GWs spectrogram.')
    name = 'GWs_spectrogram.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'GW_spectrogram.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name)
    print('GWs spectrogram plotted.')
except:
    print('Error in plotting the GWs spectrogram.')

if args.GW_contribution:
    print('Plotting GWs contribution.')
    name = 'AE220_contribution.png'
    if not os.path.exists(os.path.join(plot_folder, name)) or args.overwrite_plots:
        os.system('python ' + os.path.join(plots_scripts_path, 'GWs_contribution.py') + ' --sim-name ' + \
                args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                ' --name ' + name + ' --rxlim ' + str(args.rxlim_GWs))
    print('GWs contribution plotted.')

if args.GW_movie:
    try:
        print('Making GWs movie.')
        name = 'GWs_Ye_S'
        if not os.path.exists(os.path.join(plot_folder, name + '.webm')) or args.overwrite_plots:
            os.system('python ' + os.path.join(movies_scripts_path, 'Entropy_Ye.py') + ' --sim ' + \
                    args.sim_name + ' --sim-path ' + args.sim_path + ' --save-path ' + plot_folder +\
                    ' --name ' + name + ' --rxlim ' + str(args.rxlim_GWs))
        print('GWs movie done.')
    except:
        print('Error in making the GW movie.')
