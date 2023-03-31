import numpy as np
from typing import Literal

class ghost:
    """
    Class that allows to remove the selected number of ghost cell from
    an array.
    Input parameters:
            ghost_cells: (int) number of ghost cell of the simulation (usually 4)
    Parameters:
            ghost: (int) number of ghost cells
            r_l: (int) ghost cells at the beginning of the radius array
            r_r: (int) ghost cells at the end of the radius array
            t_l: (int) ghost cells at the beginning of the theta angle array
            t_r: (int) ghost cells at the end of the theta angle array
            p_l: (int) ghost cells at the beginning of the phi angle array
            p_r: (int) ghost cells at the end of the phi angle array
    """
    def __init__(self, ghost_cells):
        self.ghost = ghost_cells
        self.__options_default = {'r_l': self.ghost,
                                  'r_r': self.ghost,
                                  't_l': self.ghost,
                                  't_r': self.ghost,
                                  'p_l': self.ghost,
                                  'p_r': self.ghost}
        for key, value in self.__options_default.items():
            self.__setattr__(key, value)
        self.__options_1D = {'radius': [self.r_l, self.r_r],
                           'theta': [self.t_l, self.t_r],
                           'phi': [self.p_l, self.p_r]}

    def restore_default(self):
        """
        Method that allows to restore the default number of ghost cells
        """
        for key, value in self.__options_default.items():
            self.__setattr__(key, value)
        self.__options_1D = {'radius': [self.r_l, self.r_r],
                           'theta': [self.t_l, self.t_r],
                           'phi': [self.p_l, self.p_r]}

    def update_ghost_cells(self, **kwargs):
        """
        Method that allows change the default number of ghost cells.
        Optional parameters allow to select a custom number of ghost cells
        for a specific grid quantity. They must be between 0 (Keep all the cells)
        and ghost_cell (keep only physical cells)
        Parameters:
            r_l: (int) ghost cells at the beginning of the radius array
            r_r: (int) ghost cells at the end of the radius array
            t_l: (int) ghost cells at the beginning of the theta angle array
            t_r: (int) ghost cells at the end of the theta angle array
            p_l: (int) ghost cells at the beginning of the phi angle array
            p_r: (int) ghost cells at the end of the phi angle array
        """
        self.restore_default()
        updated_parameters = self.__options_default.copy()
        updated_parameters.update(kwargs)
        values = list(updated_parameters.values())
        if min(values)<0 or max(values)>self.ghost:
            raise TypeError("The number of ghost cells MUST be between 0 and " + str(self.ghost))
        for key, value in updated_parameters.items():
            self.__setattr__(key, value)
        self.__options_1D = {'radius': [self.r_l, self.r_r],
                           'theta': [self.t_l, self.t_r],
                           'phi': [self.p_l, self.p_r]}
        del updated_parameters

    def return_ghost_dictionary(self):
        """
        Method that returns the current ghost cells options
        """
        return self.__options_1D

    def remove_ghost_cells(self, array, dim, quantity_1D: 
                           Literal['radius', 'theta', 'phi'] = None):
        """
        Class method that allows to remove ghost cells from a given array
        Parameters:
            array: (numpy array) quantity from which you want to remove the ghost cells
            dim: (int) dimension of the supernova simulation
            quantity_1D: (string, optional), only for 1D array, allows to select which grid 
                         parameter to use
        Method result:
            array: (numpy array) quantity with ghost cells removed 
        """
        assert dim in (1, 2, 3), "Simulation MUST be 1, 2 or 3D"
        array_dim = array.ndim
        if dim == 1:
            if array_dim not in (1, 2, 3):
                raise TypeError("Array MUST be 1, 2 or 3D")
            if array_dim == 1:
                return self.__remove_1D_ghost_cells(array, 'radius')
            elif array_dim == 2:
                return self.__remove_ghost_cells_2D_ar_1D_sim(array)
            else:
                return self.__remove_ghost_cells_3D_ar_1D_sim(array)
        elif dim == 2:
            if array_dim not in (1, 2, 3, 4):
                raise TypeError("Array MUST be 1, 2, 3 or 4D")
            if array_dim == 1:
                if not quantity_1D in ['radius', 'theta']:
                    raise TypeError("Quantity type required: " + str(['radius', 'theta']))
                return self.__remove_1D_ghost_cells(array, quantity_1D)
            elif array_dim == 2:
                return self.__remove_2D_ghost_cells(array)
            elif array_dim == 3:
                return self.__remove_ghost_cells_3D_ar_2D_sim(array)
            else:
                return self.__remove_ghost_cells_4D_ar_2D_sim(array)
        else:
            if array_dim not in (1, 3, 4, 5):
                raise TypeError("Array MUST be 1, 3, 4 or 5D")
            if array_dim == 1:
                if not quantity_1D in ['radius', 'theta', 'phi']:
                    raise TypeError("Quantity type required: " + str(['radius', 'theta', 'phi']))
                return self.__remove_1D_ghost_cells(array, quantity_1D)
            elif array_dim == 3:
                return self.__remove_3D_ghost_cells(array)
            elif array_dim == 4:
                return self.__remove_ghost_cells_4D_ar_3D_sim(array)
            else:
                return self.__remove_ghost_cells_5D_ar_3D_sim(array)

    def remove_ghost_cells_radii(self, array, dim, **kwargs):
        assert dim in (1, 2, 3), "Simulation MUST be 1, 2 or 3D"
        
        if kwargs:
            self.update_ghost_cells(**kwargs)
            if dim == 1:
                pass
            elif dim == 2:
                array = self.__remove_2D_ghost_cells_radii(array)
            else:
                array = self.__remove_3D_ghost_cells_radii(array)
            self.restore_default()
            return array
        if dim == 1:
            return array
        if dim == 2:
            return self.__remove_2D_ghost_cells_radii(array)
        else:
            return self.__remove_3D_ghost_cells_radii(array)

    def __remove_1D_ghost_cells(self, array, quantity_1D):
        assert array.ndim == 1, "Array must be 1-dimensional"
        boundaries = self.__options_1D[quantity_1D]
        size = array.shape[0]
        return array[boundaries[0] : size - boundaries[1]]
    
    def __remove_2D_ghost_cells(self, array):
        assert array.ndim == 2, "Array must be 2-dimensional"
        size_y = array.shape[0]
        size_x = array.shape[1]
        return array[self.t_l : size_y - self.t_r,
                     self.r_l : size_x - self.r_r]

    def __remove_3D_ghost_cells(self, array):
        assert array.ndim == 3, "Array must be 3-dimensional"
        size_z = array.shape[0]
        size_y = array.shape[1]
        size_x = array.shape[2]
        return array[self.p_l : size_z - self.p_r,
                     self.t_l : size_y - self.t_r, 
                     self.r_l : size_x - self.r_r]
    
    def __remove_2D_ghost_cells_radii(self, array):
        t_r = abs(self.t_r-self.__options_default['t_r'])
        t_l = abs(self.t_l-self.__options_default['t_l'])
        return array[t_l : array.shape[0] - t_r]

    def __remove_3D_ghost_cells_radii(self, array):
        t_r = abs(self.t_r-self.__options_default['t_r'])
        t_l = abs(self.t_l-self.__options_default['t_l'])
        p_r = abs(self.p_r-self.__options_default['p_r'])
        p_l = abs(self.p_l-self.__options_default['p_l'])
        return array[p_l : array.shape[0] - p_r,
                     t_l : array.shape[1] - t_r]

    def __remove_ghost_cells_2D_ar_1D_sim(self, array):
        assert array.ndim == 2, "Array must be 2-dimensional"
        size = array.shape[0]
        return array[self.r_l : size - self.r_r, :]

    def __remove_ghost_cells_3D_ar_1D_sim(self, array):
        assert array.ndim == 3, "Array must be 3-dimensional"
        size = array.shape[0]
        return array[self.r_l : size - self.r_r, :, :]

    def __remove_ghost_cells_3D_ar_2D_sim(self, array):
        assert array.ndim == 3, "Array must be 3-dimensional"
        size_y = array.shape[0]
        size_x = array.shape[1]
        return array[self.t_l : size_y - self.t_r, self.r_l : size_x - self.r_r, :]

    def __remove_ghost_cells_4D_ar_2D_sim(self, array):
        assert array.ndim == 4, "Array must be 4-dimensional"
        size_y = array.shape[0]
        size_x = array.shape[1]
        return array[self.t_l : size_y - self.t_r, 
                     self.r_l : size_x - self.r_r,
                     :, :]

    def __remove_ghost_cells_4D_ar_3D_sim(self, array):
        assert array.ndim == 4, "Array must be 4-dimensional"
        size_z = array.shape[0]
        size_y = array.shape[1]
        size_x = array.shape[2]
        return array[self.p_l : size_z - self.p_r,
                     self.t_l : size_y - self.t_r, 
                     self.r_l : size_x - self.r_r,
                     :]
    
    def __remove_ghost_cells_5D_ar_3D_sim(self, array):
        assert array.ndim == 5, "Array must be 5-dimensional"
        size_z = array.shape[0]
        size_y = array.shape[1]
        size_x = array.shape[2]
        return array[self.p_l : size_z - self.p_r,
                     self.t_l : size_y - self.t_r, 
                     self.r_l : size_x - self.r_r,
                     :, :]
