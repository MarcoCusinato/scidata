import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),'../..'))
import numpy as np
from scidata.math_functions.IDL_derivative import IDL_derivative
from scidata.units.units import units

u = units()



def calculate_AE220(simulation):
    """
    Calculates the AE220 from density and velocities for a
    2D simulation.
    ONLY 2D
    Returns
        time: arry of time step
        AE220: len(radius), len(time) array
    """
    radius = simulation.cell.radius(simulation.ghost)
    dV = -simulation.cell.dVolume_integration(simulation.ghost)
    costheta = np.cos(simulation.cell.theta(simulation.ghost))
    file_list = simulation.file_list_hdf()
    inner_radius, time, gh_inner = simulation.get_innercore_radius(innercore_radius=True, ret_time=True,
                                                    ghost_cells=True)
    convective_radius, gh_conv = simulation.get_convective_radius(convective_radius=True,
                                                            ghost_cells=True)
    NE220_full = np.zeros((len(radius), len(time)))
    NE220_inner = np.zeros((len(radius), len(time)))
    NE220_convection = np.zeros((len(radius), len(time)))
    NE220_outer = np.zeros((len(radius), len(time)))
    radius = radius[None, ...]
    costheta = costheta[..., None]
    for index in range(len(file_list)):
        # Get all the quantities needed for the calculation
        try:
            data_h5 = simulation.open_h5(file_list[index])
        except:
            NE220_full[:, index] = NE220_full[:, index - 1]
            NE220_inner[:, index] = NE220_inner[:, index - 1]
            NE220_convection[:, index] = NE220_convection[:, index - 1]
            NE220_outer[:, index] = NE220_outer[:, index - 1]
            continue
        rho = simulation.rho(data_h5)
        vR = simulation.radial_velocity(data_h5)
        vT = simulation.theta_velocity(data_h5)
        simulation.close_h5(data_h5)
        # Calculate the NE220
        # NE220 of the full space
        NE220 = dV * radius * rho *  ( vR * \
            ( 3 * costheta ** 2 - 1 ) - 3 * vT * costheta * np.sqrt( 1 - costheta ** 2 ) )
        # Create masks
        mask_inner = radius <= simulation.ghost.remove_ghost_cells_radii(convective_radius[..., index],
                                                                        simulation.dim, **gh_conv)[..., None]
        mask_convection = (radius <= simulation.ghost.remove_ghost_cells_radii(inner_radius[..., index],
                                                                            simulation.dim, **gh_inner)[..., None] + 2e6) & \
                            ( np.logical_not(mask_inner) )
        mask_outer = np.logical_not(mask_inner + mask_convection)
        # Calculate the NE220s
        NE220_full[:, index] = np. sum( NE220, axis = 0 )
        NE220_inner[:, index] = np. sum( NE220 * mask_inner, axis = 0 )
        NE220_convection[:, index] = np. sum( NE220 * mask_convection, axis = 0 )
        NE220_outer[:, index] = np. sum( NE220 * mask_outer, axis = 0 )
    
    const = ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )

    return time, IDL_derivative(time, NE220_full * const), calculate_strain(NE220_full, time), \
            calculate_strain(NE220_inner, time), calculate_strain(NE220_convection, time), \
            calculate_strain(NE220_outer, time)
            

def calculate_strain(NE220, time):
    """
    Calculates the strain from the NE220
    """
    # Put again the correct units
    NE220 *= ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )
    # Calculate the strain
    strain = IDL_derivative(time, NE220)
    return strain.sum(axis=0)
