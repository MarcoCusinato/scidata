import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../'))
from scidata.cell.cell import cell as cl
from scidata.platform.platform import pltf

dim_keys = {1: '1D', 2: '2D', 3: '3D'}

#functions
def windows_mount():
    """
    Creates a dictionary in a file that specifies on which Windows partition
    the server is mounted.
    """
    message = "We detected you are working on windows, but no" + \
              " mounted partition has been assigned."
    print(message)
    mount_path = input("Please insert the mount path (i.e. Z:): ")
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
              'paths/windows.py'), 'w') as f:
        f.write('windows_paths = ' + str({'mount_path': mount_path}))

def standardize_simulation_path(path):
    """
    Returns standardized Linux path (i.e. path starting with '/')
        in: path string
    """
    if pltf() == 'Windows':
        try:
            from scidata.paths.paths.windows import windows_paths
        except:
            windows_mount()
            from scidata.paths.paths.windows import windows_paths
        path = path.replace('\\', '/')
        if path.startswith(windows_paths['mount_path']):
            path = path[path.find('/'):]
    else:
        path = path.replace('\\', '/')
    return path

def read_path_file(path):
    """
    Returns paths from a file in list format
    in:
        path: string, path to file
    """
    if not os.path.exists(path):
        return None
    import numpy
    return list(numpy.loadtxt(path, dtype=str))

def list_of_paths(files_path, single_paths):
    """
    Returns list of all simulation paths to include
    in:
        files_path: list of strings, path to the file containing the simulations paths.
        single_paths: list of the simulation paths to include.
    """
    if single_paths:
        path_list = single_paths
    else:
        path_list = []
    if files_path:
        for p in files_path:
            list_path_from_file = read_path_file(p)
            for pp in list_path_from_file:
                path_list.append(pp)
    while None in path_list:
        path_list.remove(None)
    print("Paths to search: ", path_list)
    return path_list

def remove_not_useful_folders(dirs):
    """
    Removes the folders containing the equations of state and initial models
    from a list of directories
    """
    dirs_to_remove = ['EOS', 'Initial_models', 'presn_models',
                     'Stellarcollapse_EOS', 'CompOSE_EOS']
    for dir in dirs_to_remove:
        while dir in dirs:
            dirs.remove(dir)
    return dirs

def check_subpaths_for_simulations(path, simDict, simulations):
    """
    Recursively checks a folder and its subfolders to ensure the presence
    of a simulation. Since every simulation contains a 'outp-hdf' folder
    we use its presence as a criterion. If a simulation is fuond, it is stored
    in the simulations dictionary.
    """
    try:
        dirs = remove_not_useful_folders(os.listdir(path))
    except:
        return simDict, simulations
    for folder in dirs:
        path_subfolder = os.path.join(path, folder)
        print('\t', path_subfolder)
        try:
            subfolders = remove_not_useful_folders(os.listdir(path_subfolder))
        except:
            continue
        if 'outp-hdf' in subfolders:
            standard_path = standardize_simulation_path(path)
            try:
                key = dim_keys[cl(path_subfolder).dim]
                if folder in simDict[key]:
                    if simDict[key][folder] == standard_path:
                        continue
                    print(folder, ' simulation already exists. The other is in: ',
                          simDict[key][folder], '\nThis one is in:', standard_path)
                    choice = None
                    msg = "Would you like to keep the new one? (y/n)\nWarning: by typing \'y\'" + \
                        "you will overwrite the old one. "
                    while choice not in ['y', 'n']:
                        choice = input(msg)
                    if choice == 'n':
                        continue
                simDict[key][folder] = standard_path
                simulations.append(folder)
            except:
                continue
        else:
            simDict, simulations = check_subpaths_for_simulations(path_subfolder,
                                                            simDict, simulations)
    return simDict, simulations

def add_paths(path_list):
    """
    Creates or updates the simulations dictionary of all the simulations
    found in a list of paths.
    """
    simDict = get_paths_dictionary(True)
    for p in path_list:
        print('Path: ', p)
        print('Checking...')
        simDict, simulations = check_subpaths_for_simulations(p, simDict, [])
        if len(simulations) == 1:
            print('Added 1 simulation.')
        elif len(simulations) == 0:
            print('No simulation added.')
        else:
            print('Added', len(simulations), 'simulations.')
        save_paths_dictionary(simDict)

def get_paths_dictionary(creating = False):
    """
    Returns a simulations dictionary. Simulations are grouped in 1D, 2D and 3D.
    """
    try:
        from scidata.paths.paths.simulations_dictionary import simulations_dictionary
    except:
        if creating:
            simulations_dictionary = {'1D': {}, '2D': {}, '3D': {}}
        else:
            raise TypeError("Dictionary of simulations not found, " + \
                        "please create one or supply simulation's path.")
    return simulations_dictionary
       
def save_paths_dictionary(simulations_dictionary):
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
              'paths/simulations_dictionary.py'), 'w') as f:
        f.write('simulations_dictionary = ' + str(simulations_dictionary))

def return_real_path(folder, path, platform):
    if platform not in ['Linux', 'Windows']:
        raise TypeError("OS not supported. Supported platforms: Linux, Windows")
    if platform == 'Windows':
        try:
            from scidata.paths.paths.windows import windows_paths
        except:
            windows_mount()
            from scidata.paths.paths.windows import windows_paths
        return os.path.join(windows_paths['mount_path'], path, folder)
    else:
        return os.path.join(path, folder)

def find_simulation(folder, platform, sim_path=None):
    """
    Funcion that allows the user to find the path of a specific simulation.
    Parameters:
        folder: (string) the name of the simulation to find
        platform: (string) the OS from which the script is running
        sim_path: (string or None) the path of the chosen simulation
    """
    if sim_path is not None:
        path = os.path.join(sim_path, folder)
        if not os.path.exists(path):
            raise TypeError("Requested simulation (" + folder + ") does not exist in the supplied path.")
        return os.path.join(sim_path, folder)
    
    simDict = get_paths_dictionary()
    mask = [False, False, False]
    sim_possibilities = ['1D', '2D', '3D']
    if folder in simDict['1D']:
        mask[0] = True
    if folder in simDict['2D']:
        mask[1] = True
    if folder in simDict['3D']:
        mask[2] = True
    sim_number = sum(mask)
    if sim_number == 0:
        msg = "Requested simulation (" + folder + ") is not in the search paths. " +\
              "Please supply also a simulation path or update your search paths."
        raise TypeError(msg)
    simDict = list(simDict.values())
    if sim_number == 1:
        import numpy as np
        simDict = np.array(simDict)
        p = simDict[mask][0][folder]
        print("Your simulation (" + folder + ") is in ", p)
        msg = "If this is not the path you have chosen, " +\
              "please turn to the dark side (supply a simulation path)."
        print(msg)
        return return_real_path(folder, p, platform)
    else:
        import numpy as np
        msg = "There are " + str(sim_number) + " simulations with the requested name (" + \
               folder + "). " + \
              "The simulations are in " + str(np.array(sim_possibilities)[mask]) + "."
        print(msg)
        choice = None
        msg = "Please insert now the simulation dimension (1D, 2D or 3D): "
        while choice not in np.array(sim_possibilities)[mask]:
            choice = input(msg)
        p = return_real_path(folder, simDict[sim_possibilities.index(choice)][folder],
                             platform)
        print("Your simulation (" + folder + ") is in ", p)
        msg = "If this is not the path you have chosen, " +\
              "please turn to the dark side (supply a simulation path)."
        print(msg)
        return p

def clear_folder(folder_path):
    """
    Deletes all empty directories in a folder
    """
    dirs = os.listdir(folder_path)
    for folder in dirs:
        path = os.path.join(folder_path, folder)
        try:
            if not os.listdir(path):
                os.rmdir(path)
        except:
            continue

def clear_simulation_folder(folder_list):
    """
    Deletes all empty folders in a list of simulations.
    """
    for folder in folder_list:
        simulation_path = find_simulation(folder, pltf())
        clear_folder(simulation_path)
