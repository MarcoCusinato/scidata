import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),'../..'))
from scidata.math_functions.IDL_derivative import IDL_derivative
import numpy as np

def GW_strain(sim_dim, column_change, data, index):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_strain_1D(data)
    elif sim_dim == 2:
        return correct_zero(2, GW_strain_2D(data), index)
    else:
        GWs = GW_strain_3D(data)
        if column_change is not None:
            GWs[:column_change, 1] = GW_strain_2D(data)[:column_change, 1]
            GWs[:column_change,2:] = 0
        return correct_zero(3, GWs, index)

def GW_strain_1D(data):
    print("No GW for you :'(")
    return None

def GW_strain_2D(data):
    """
        GWs amplitudes calculated as the first partial time derivative
        of NE_220 and not with the second partial time derivative of ME_220.
        Moreover we only consider matter contribution to the amplitude, being
        the main contribution to it.
    """
    const = -0.125 *  np.sqrt(15/np.pi)
    return np.stack((data[:, 2], const * IDL_derivative(data[:,2], data[:,5])),
              axis = 1)

def GW_strain_3D(data):
    const = 1 #8 * np.sqrt(np.pi / 15)
    hD_pl_p = 2. * ( data[:,9] - data[:,13] )
    hD_pl_e = 2. * ( data[:,17] - data[:,13] )
    hD_cr_p = 2. * ( data[:,10] + data[:,12] )
    hD_cr_e = 2. * ( - data[:,14] - data[:,16] )

    hD_pl_p = const * IDL_derivative( data[:,2], hD_pl_p )
    hD_cr_p = const * IDL_derivative( data[:,2], hD_cr_p )
    hD_pl_e = const * IDL_derivative( data[:,2], hD_pl_e )
    hD_cr_e = const * IDL_derivative( data[:,2], hD_cr_e )

    return np.stack((data[:,2], hD_pl_e, hD_pl_p, hD_cr_e, hD_cr_p), axis = -1)

def correct_zero(sim_dim, GWs, index):
    if sim_dim == 1:
        pass
    else:
        if index is None:
            return GWs
        for i in range(1, GWs.shape[1]):
            GWs[:, i] -= GWs[:index, i].mean()
        return GWs
