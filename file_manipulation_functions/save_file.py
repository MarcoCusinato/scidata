
import h5py

dict_flavour = {'0': 'nue', '1': 'nua', '2': 'nux'}

def save_h5(save_path, data_radius, data_average,
            indices, ghost_cells, corrected_for_tob):
        """
        Save data in h5 format
        """
        file_out = h5py.File(save_path, 'w')
        if 'neutrino' in save_path:
                neutrino_save(file_out, data_radius, data_average, indices,
                              ghost_cells, corrected_for_tob)
        else:
                standard_save(file_out, data_radius, data_average, indices,
                              ghost_cells, corrected_for_tob)
        file_out.close()

def standard_save(file_out, data_radius, data_average, indices, ghost_cells,
                  corrected_for_tob):
        """
        Save radii data:
                Contains the following dataset
                - radii: NumPy array with dimensions (phi, theta, time),
                         theta and phi with ghost cells (default 1),
                         contains the raidius value for every cell
                - time: NumPy array with dimension (time), contains time valuse for each radius value
                - min: NumPy array with dimension (time), contains minimum radius value
                - max: NumPy array with dimension (time), contains maximum radius value
                - average: NumPy array with dimension (time), contains average radius value
                - indices: NumPy array with dimension (time), contains file indices
                - tob_correction: bool, if time value is tob corrected
        """
        file_out.create_dataset('radii', data = data_radius)
        file_out.create_dataset('time', data = data_average[:, 0])
        file_out.create_dataset('min', data = data_average[:, 1])
        file_out.create_dataset('max', data = data_average[:, 2])
        file_out.create_dataset('average', data = data_average[:, 3])
        file_out.create_dataset('indices', data = indices)
        file_out.create_dataset('tob_correction', data = corrected_for_tob)
        g_group = file_out.create_group('ghost_cells')
        for (key, value) in ghost_cells.items():
                g_group.create_dataset(key, data = value)

def neutrino_save(file_out, data_radius, data_average, indices, ghost_cells,
                  corrected_for_tob):
        """
        Save nutrino sphere radii data:
                Contains three groups, one for each neutrino flavour, with
                the following dataset
                - radii: NumPy array with dimensions (phi, theta, time),
                         theta and phi with ghost cells (default 1),
                         contains the raidius value for every cell
                - time: NumPy array with dimension (time), contains time valuse for each radius value
                - min: NumPy array with dimension (time), contains minimum radius value
                - max: NumPy array with dimension (time), contains maximum radius value
                - average: NumPy array with dimension (time), contains average radius value
                - indices: NumPy array with dimension (time), contains file indices
                - tob_correction: bool, if time value is tob corrected
        """
        file_out.create_dataset('time', data = data_average[:, 0])
        nu_group = file_out.create_group('radii')
        max_group = file_out.create_group('max')
        min_group = file_out.create_group('min')
        average_group = file_out.create_group('average')
        for flavour in range(3):
                key = dict_flavour[str(flavour)]
                nu_group.create_dataset(key, data = data_radius[..., flavour, :])
                max_group.create_dataset(key, data = data_average[:, 3 * flavour + 1])
                min_group.create_dataset(key, data = data_average[:, 3 * flavour + 2])
                average_group.create_dataset(key, data = data_average[:, 3 * flavour + 3])
        file_out.create_dataset('indices', data = indices)
        file_out.create_dataset('tob_correction', data = corrected_for_tob)
        g_group = file_out.create_group('ghost_cells')
        for (key, value) in ghost_cells.items():
                g_group.create_dataset(key, data = value)
