import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),'../..'))
from scidata.math_functions.IDL_derivative import IDL_derivative
import numpy as np

def GW_strain(sim_dim, data):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_strain_1D(data)
    if sim_dim == 2:
        return GW_strain_2D(data)
    return GW_strain_3D(data)

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

