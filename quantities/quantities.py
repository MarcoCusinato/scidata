import os
import sys
import warnings
from typing import Literal

import h5py
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import stft
from scipy.fft import fft, fftfreq

dirname = os.path.dirname(__file__)
sys.path.append(os.path.join(dirname,'../..'))
from scidata.units.units import units

from scidata.cell.cell import cell as cl
from scidata.cell.ghost import ghost as gh
from scidata.file_manipulation_functions.load_file import load_file, find_column_changing_line
from scidata.file_manipulation_functions.save_file import save_h5
from scidata.grid.grid import grid as gr
from scidata.legendre.legendre_polynomials import legendre_polynomials as LegP
from scidata.math_functions.average_profile import average_profile
from scidata.math_functions.IDL_derivative import IDL_derivative
from scidata.math_functions.streamlines import strfct2D
from scidata.par.read_parfile import get_indices_from_parfile as getPar
from scidata.paths.existence import check_path
from scidata.paths.path_managment import find_simulation
from scidata.platform.platform import local_storage_folder, pltf
from scidata.GWs_strain.GWs import GW_strain
from scidata.GWs_strain.spectrogram import GW_spectrogram

u = units()

class SimulationAnalysis:
    def __init__(self, simulation_name, dim = None, simulation_folder_path=None):
        platform = pltf()
        self.simulation_name = simulation_name
        self.path = find_simulation(simulation_name, platform, simulation_folder_path)
        self.log_path = os.path.join(self.path, 'log')
        self.hdf_path = os.path.join(self.path, 'outp-hdf')
        self.grid_path = os.path.join(self.path, 'grid')
        self.par_path = os.path.join(self.path, 'pars')
        self.integrated_nu = 'neu.dat'
        self.rho_max_file = 'rho.dat'
        self.grw = 'grw.dat'
        self.nu_e_grid = 'grid.e.dat'
        self.hydroTHD_index, self.ghost_cells = getPar('start.pars', self.par_path)
        self.cell = cl(self.path, dim)
        self.dim = self.cell.simulation_dimension()
        self.ghost = gh(self.ghost_cells)
        self.storage_path = local_storage_folder(platform, simulation_name, self.dim)
    
    #h5 methods
    def open_h5(self, file_name):
        file_path = os.path.join(self.hdf_path, file_name)
        return h5py.File(file_path)
    
    def close_h5(self, data_h5):
        data_h5.close()
    
    #hydro quatities
    def rho(self, data_h5):
        data = np.array(data_h5['hydro']['data'])[:,:,:,self.hydroTHD_index['hydro']['I_RH']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
   
    def energy_MHD(self, data_h5):
        data = np.array(data_h5['hydro']['data'])[:,:,:,self.hydroTHD_index['hydro']['I_EN']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def Ye(self, data_h5):
        data = np.array(data_h5['hydro']['data'])[:,:,:,self.hydroTHD_index['hydro']['I_YE']]
        data = self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
        rho = self.rho(data_h5)
        return data/rho
    
    #thd quatities
    def gas_pressure(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_PGAS']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def temperature(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_TMPR']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def nu_heat(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_HEAT']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def radial_velocity(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_VELX']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def theta_velocity(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_VELY']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def phi_velocity(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_VELZ']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def omega(self, data_h5):
        v_phi = self.phi_velocity(data_h5)
        radius = self.cell.radius(self.ghost)
        if self.dim == 1:
            return v_phi / radius[None, :]
        theta = np.sin(self.cell.theta(self.ghost))
        if self.dim == 2:
            return v_phi / (theta[:, None] * radius[None, :])
        if self.dim == 3:
            phi = np.cos(self.cell.phi(self.ghost))
            return v_phi / (phi[:, None, None] * theta[None, :, None] * radius[None, None, :])

    def speed_of_sound(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_CSND']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def entropy(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_ENTR']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def internal_energy(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_EINT']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def lorentz(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_LRTZ']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def chemical_potential_electrons(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_CPOT'][0]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def chemical_potential_neutrons(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_CPOT'][1]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def chemical_potential_protons(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_CPOT'][2]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def chemical_potential_neutrinos(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_CPOT'][3]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def neutron_fraction(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_COMP'][0]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def proton_fraction(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_COMP'][1]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def alpha_fraction(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_COMP'][2]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def heavy_nuclei_fraction(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_COMP'][3]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def heavy_nuclei_average_mass(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_COMP'][4]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def heavy_nuclei_average_proton_number(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_COMP'][5]]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def eta_electrons(self, data_h5):
        chem_pot = self.chemical_potential_electrons(data_h5)
        temp = self.temperature(data_h5)
        return chem_pot/temp
    
    def adiabatic_index(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_GAMM']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def eta_neutrons(self, data_h5):
        chem_pot = self.chemical_potential_neutrons(data_h5)
        temp = self.temperature(data_h5)
        return chem_pot/temp

    def eta_protons(self, data_h5):
        chem_pot = self.chemical_potential_protons(data_h5)
        temp = self.temperature(data_h5)
        return chem_pot/temp

    def eta_neutrinos(self, data_h5):
        chem_pot = self.chemical_potential_neutrinos(data_h5)
        temp = self.temperature(data_h5)
        return chem_pot/temp
    
    #error
    def errors(self, data_h5):
        data = np.array(data_h5['thd']['data'])[...,self.hydroTHD_index['thd']['I_EOSERR']]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    #neutrinos
    def neutrino_energy_grid(self):
        bin_energies = np.loadtxt(os.path.join(self.grid_path, self.nu_e_grid))
        return bin_energies[:,2]
    
    def neutrino_dE(self):
        bin_energies = np.loadtxt(os.path.join(self.grid_path, self.nu_e_grid))
        return bin_energies[:,3] - bin_energies[:,1]

    def neutrino_flux_tot(self, data_h5):
        """
        Last index is the neutrino flavour: nue, nua, nux.
        Second-last index is the energy bin index.
        """
        nu_out = np.array(data_h5['neutrino']['e'])[..., 1:]
        
        nu_out = nu_out.sum(axis=-1)
        nu_out[..., 2] /= 4

        return self.ghost.remove_ghost_cells(np.squeeze(nu_out), self.dim)

    def neutrino_energy_density(self, data_h5):
        """
        Returns the 0th moment (energy density) into neutrino flavour 
        """
        nu_out = np.array(data_h5['neutrino']['e'])[..., 0]
        nu_out[..., 2] /= 4
        return self.ghost.remove_ghost_cells(np.squeeze(nu_out), self.dim)
    
    def neutrino_energy_density_bin_integrated(self, data_h5):
        nu_out = self.neutrino_energy_density(data_h5)
        return nu_out.sum(axis=self.dim)

    def neutrino_fluxes(self, data_h5):
        """
        Returns the neutrino fluxes per direction and neutrino flavour
        Last index: direction, 0:r 1:theta 2:phi
        second to last neutrino flavour
        """
        nu_out = np.array(data_h5['neutrino']['e'])[..., 1:]
        nu_out[..., 2, :] /= 4
        return self.ghost.remove_ghost_cells(np.squeeze(nu_out), self.dim)

    def neutrino_flux_tot_bin_integrated(self, data_h5):
        """
        Last index is the neutrino flavour: nue, nua, nux
        """
        nu = self.neutrino_flux_tot(data_h5)
        return nu.sum(axis = self.dim)

    def neutrino_number_flux_tot(self, data_h5, ret_nu_flux = False):
        """
        Last index is the neutrino flavour: nue, nua, nux
        """
        nu = self.neutrino_flux_tot(data_h5)
        bin_energies = 1/u.convert_to_erg(self.neutrino_energy_grid())
        bin_energies = bin_energies[None, :, None]
        while bin_energies.ndim != nu.ndim:
            bin_energies = bin_energies[None, ...]
        if ret_nu_flux:
            return nu, nu*bin_energies
        return nu*bin_energies

    def neutrino_opacity_flux(self, data_h5):
        """
        Last index is the neutrino flavour: nue, nua, nux
        """
        data = np.array(data_h5['neutrino']['oe'])[...,1:]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def neutrino_luminosity_profile(self, data_h5):
        """
        Last index differentiate grey neutrino luminosities and number luminosities
        indices
        1: grey nue luminosity  4: grey nue number luminosity
        2: grey nua luminosity  5: grey nua number luminosity
        3: grey nux luminosity  6: grey nux number luminosity
        """
        dV = self.cell.dVolume_sum(self.ghost)
        neutrino_fl, neutrino_num_fl = self.neutrino_number_flux_tot(data_h5, True)
        neutrino_fl *= dV[..., None, None]
        neutrino_num_fl *= dV[..., None, None]
        return np.concatenate((neutrino_fl.sum(axis = self.dim),
                               neutrino_num_fl.sum(axis = self.dim)),
                               axis = self.dim)

    def neutrino_mean_energy_profile(self, data_h5):
        """
        Computation of neutrino mean energies
        indices
        1: nue mean energy
        2: nua mean energy
        3: nux mean energy
        """
        nu = self.neutrino_luminosity_profile(data_h5)
        return u.convert_to_MeV(nu[..., :3] / nu[..., 3:])

    def neutrino_averaged_radial_profile(self, data_h5):
        """
        index 1 to 3: max of nue, nua, nux radially
        index 4 to 6: average of nue nua nux radially
        """
        data = self.neutrino_luminosity_profile(data_h5)
        return np.stack((data[..., 0].reshape(-1, data[..., 0].shape[-1]).max(axis = 0),
                         data[..., 1].reshape(-1, data[..., 1].shape[-1]).max(axis = 0),
                         data[..., 2].reshape(-1, data[..., 2].shape[-1]).max(axis = 0),
                         data[..., 0].reshape(-1, data[..., 0].shape[-1]).mean(axis = 0),
                         data[..., 1].reshape(-1, data[..., 1].shape[-1]).mean(axis = 0),
                         data[..., 2].reshape(-1, data[..., 2].shape[-1]).mean(axis = 0)), axis=1)

    #Profiles
    def average_rho_profile(self, data_h5):
        return average_profile(self.dim, self.rho(data_h5))
    
    def average_radial_velocity_profile(self, data_h5):
        return average_profile(self.dim, self.radial_velocity(data_h5))

    def average_theta_velocity_profile(self, data_h5):
        return average_profile(self.dim, self.theta_velocity(data_h5))
    
    def average_phi_velocity_profile(self, data_h5):
        return average_profile(self.dim, self.phi_velocity(data_h5))

    #magnetic fields
    def magnetic_field_ct(self, data_h5):
        """
        Magnetic field at the cells border, to use only to calculate streamlines.
        If you want to plot the actual magnetic fields use the 'magnetic_field'
        method.
        """
        data = np.array(data_h5['mag_CT']['data'])
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)

    def magnetic_field(self, data_h5):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(data_h5['mag_vol']['data'])),
                                             self.dim) 

    def poloidal_magnetic_field(self, data_h5):
        data = self.magnetic_field(data_h5)
        return np.sqrt(data[...,0]**2+data[...,1]**2)

    def toroidal_magnetic_field(self, data_h5):
        return self.magnetic_field(data_h5)[...,2]  
    
    def magnetic_energy_per_volum_unit_components(self, data_h5):
        data = self.magnetic_field(data_h5)
        data[..., 0] = 0.5 * data[..., 0]**2 
        data[..., 1] = 0.5 * data[..., 1]**2 
        data[..., 2] = 0.5 * data[..., 2]**2 
        return data
    
    def magnetic_energy_per_volum_unit_poloidal(self, data_h5):
        data = self.magnetic_field(data_h5)
        return 0.5*(data[..., 0]**2 + data[..., 1]**2) 

    def magnetic_energy_per_volum_unit_toroidal(self, data_h5):
        data = self.magnetic_field(data_h5)
        return 0.5 * data[..., 2]**2
    
    def magnetic_energy_per_unit_volume(self, data_h5):
        data = self.magnetic_field(data_h5)
        return 0.5*(data[..., 0]**2 + data[..., 1]**2 + data[..., 2]**2)
    
    ##---------------------------------------------------------------------------------------
    ## Gravitational Wavess
    ##---------------------------------------------------------------------------------------

    def GW_Amplitudes(self, correct_for_tob=True, zero_correction=True, lower_refinement=True):
        """
        Params:
            zero_correction: shifts up or down the amplitude to make it centred with zero
            lower_refinement: on fine timestep refinements GWs are usually noisy, so we take 
            one every n points. n defined as the number to reach 0.1 ms
        Returns GWs amplitudes:
        1D:
            No GWs in spherical symmetry
        2D:
            Column 1: time
            Column 2: + ploarization
        3D:
            Column 1: time
            Column 2: + polarization equatorial plane
            Column 3: + polarization polar plane
            Column 4: x polarization equatorial plane
            Column 5: x polarization polar plane
        """
        data = load_file(self.log_path, self.grw)
        if lower_refinement:
            dt = data[1, 2] - data[0, 2]
            new_dt = dt
            n=1
            while new_dt < 5e-5:
                new_dt += dt
                n += 1
            data = data[::n, :]

        column_change = find_column_changing_line(self.log_path, self.grw)
        if zero_correction:
            index = np.argmax((data[:, 2] - self.time_of_bounce_rho())  >= -0.01)
        else:
            index = None
        if correct_for_tob and zero_correction:
            data[:,2] -= self.time_of_bounce_rho()
        
        return GW_strain(self.dim, column_change, data, index)
    
    def GW_dimensionless_strain(self, distance, angle = 0.5 * np.pi, correct_for_tob=True):
        """
        Returns GWs amplitudes:
        1D:
            No GWs in spherical symmetry
        2D:
            Column 1: time
            Column 2: + ploarization
        3D:
            Column 1: time
            Column 2: + polarization equatorial plane
            Column 3: + polarization polar plane
            Column 4: x polarization equatorial plane
            Column 5: x polarization polar plane
        """
        data = self.GW_Amplitudes(correct_for_tob)
        if self.dim == 1:
            return data
        const = np.sqrt(np.pi/15)*np.sin(angle)**2/distance
        data[:, 1:] *= const
        return data
    
    def AE220(self, correct_for_tob = True):
        """
        Calculates the AE220 from density and velocities for a
        2D simulation.
        ONLY 2D
        Returns
            radius
            time: array of time step
            AE220: len(radius), len(time) array
        """
        AE220_file = os.path.join(self.storage_path,'AE220.h5')
        if not os.path.exists(AE220_file):
            warnings.warn("AE220 file not found. Creating one. \nPlease wait...")
            time, AE220 = self.__AE220()
            AE220_hdf = h5py.File(AE220_file, 'w')
            AE220_hdf.create_dataset('time', data = time)
            AE220_hdf.create_dataset('AE220', data = AE220)
            AE220_hdf.close()
        data = self.open_h5(AE220_file)
        time = data["time"][...]
        AE220 = data["AE220"][...]
        self.close_h5(data)
        const = -0.125 *  np.sqrt(15/np.pi)
        if not correct_for_tob:
            time += self.time_of_bounce_rho()
        return self.cell.radius(self.ghost), time, AE220 * const


    def __AE220(self):
        """
        Calculates the AE220 from density and velocities for a
        2D simulation.
        ONLY 2D
        Returns
            time: arry of time step
            AE220: len(radius), len(time) array
        """
        radius = self.cell.radius(self.ghost)
        dV = -self.cell.dVolume_integration(self.ghost)
        costheta = np.cos(self.cell.theta(self.ghost))
        file_list = self.file_list_hdf()
        time = np.zeros(len(file_list))
        NE220 = np.zeros((len(radius), len(time)))
        radius = radius[None, ...]
        costheta = costheta[..., None]
        for index in range(len(file_list)):
            data_h5 = self.open_h5(file_list[index])
            rho = self.rho(data_h5)
            vR = self.radial_velocity(data_h5)
            vT = self.theta_velocity(data_h5)
            time[index] = self.time(data_h5)
            self.close_h5(data_h5)
            NE220[:, index] = np. sum( dV * radius * rho *  ( vR * \
                ( 3 * costheta ** 2 - 1 ) - 3 * vT * costheta * np.sqrt( 1 - costheta ** 2 ) ), 
                axis = 0 )
        NE220 *= ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )
            
        return time - self.time_of_bounce_rho(), IDL_derivative(time, NE220)
        
    
    def GW_spectrogram(self, window_size = 10, GW_norm = 1, correct_for_tob = True):
        """
        Parameters:
            window_size: value of the time window to use in ms
            GW_norm: factor that multiplies the GW strain
            correct_for_tob: if the returned timeseries has to be corrected for the tob
        Returns:
            time: timeseries in s
            frequency: aray of the frequencies in Hz
            Zxx: magnitude
        In 3D simulations:
            h_pl_e, h_pl_p, h_cr_e, h_cr_p
        """
        GW_strain = self.GW_Amplitudes(correct_for_tob=False)
        GW_strain[:, 1:] *= GW_norm
        window = 0
        while np.abs(GW_strain[window,0] - GW_strain[0,0]) < u.convert_to_s(window_size):
            window += 1
        
        time, frequency, Zxx = GW_spectrogram(self.dim, GW_strain, window)
        if correct_for_tob:
            time -= self.time_of_bounce_rho()
        return time, frequency, Zxx
    
    def Deltah(self, peak:Literal['bounce', 'highest'],
               interval = [None, None], min_time=1.75, max_time=2, use_derivative = False, 
               return_coordinates = False):
        """
        Returns the Delta h of the gravitational wave strain as defined in
        Richers et al. 2017 (https://arxiv.org/pdf/1701.02752.pdf).
        Basically the difference between the first maximum and minimum postbounce,
        the first peak that appears is not considered.
        In case the highest peak is selected the amplitude returned is the maxima
        between left and right.
        min and max time are windows in which to search the low and high peaks in the strain
        Return:
            amplitude in cm
            if return coordinates
                time, h of the highest and lowest peak
        """
        GWs = self.GW_Amplitudes()
        indices = self.__GWs_peak_indices(GWs, peak, interval, min_time, max_time, use_derivative)
        deltah = np.abs(GWs[indices[1], 1] - GWs[indices[2], 1])
        if return_coordinates:
            x = [GWs[indices[1], 0], GWs[indices[2], 0]]
            y = [GWs[indices[1], 1], GWs[indices[2], 1]]
            return deltah, np.array(x), np.array(y)
        else:
            return deltah


    def GWs_peak(self, peak:Literal['bounce', 'highest'] = 'bounce',
                 interval = [None, None], min_time=1.75, max_time=2, use_der=False,
                 return_time = False):
        """
        This method calculate the value of the GWs strain peak and its position in time
        Return:
            Peak value
            OR
            time, peak
        """
        GWs = self.GW_Amplitudes(True)
        indices = self.__GWs_peak_indices(GWs, peak, interval, min_time, max_time, use_der)
        if return_time:
            return GWs[indices[1], 0], GWs[indices[1], 1]
        else:
            return GWs[indices[1], 1]

    def GWs_fourier_transform(self, peak:Literal['bounce', 'highest'] = 'bounce',
                              min_time=1.75, max_time=2, use_der=False,
                              interval = [None, None]):
        """
        Method that calculates the fourier transform of a specific GWs peak.
        The rest of the array will be padded with zeros to increase the resolution
        Return:
            frequency array
            fourier transformed strain x sqrt{frequency}
        """
        GWs = self.GW_Amplitudes(True)
        indices = self.__GWs_peak_indices(GWs, peak, interval, min_time, max_time, use_der)
        return self.__GWs_fourier_transform(GWs, indices)
    
    def GWs_peak_frequencies(self, peak:Literal['bounce', 'highest'] = 'bounce',
                              min_time=1.75, max_time=2, use_der=False,
                              interval = [None, None], return_intensities = False, return_fourier = False):
        """
        Calculates the dominant and the second dominant frequency of a GWs peak
        Return:
            frequencies: dominat, second dominant
            if return intensities
                intensities: dominant, second dominant
            if return return fourier
                frequencies array
                htilde
        """
        frequency, htilde = self.GWs_fourier_transform(peak, min_time, max_time, use_der, interval)
        indices = self.__GWs_frequency_peak_indices(frequency, htilde)
        return_list = [frequency[indices]]
        if return_intensities:
            return_list.append(htilde[indices])
        if return_fourier:
            return_list.append([frequency, htilde])
        if len(return_list) == 1:
            return return_list[0]
        return return_list
            

    def __GWs_peak_indices(self, GWs, peak, interval, min_time, max_time, use_der):
        """
        Method that finds the coordinates of the minimum and maximum peak of the GWs strain as well as 
        the coordinates of the points before and after the oscillation. Namely the latter are the
        first intersection point with the x axis after the peak and the third one before it.
        Parameters:
            peak: which peak to find, can be the bounce peak, highest in an interval
            interval: interval (in ms) in which the peak has to be found, if only one value is provided, that
                    would be used as the right hand side
        Returns:
            list containing left index, peak index, right index
        """
        zeros = np.where(GWs[:-1, 1] * GWs[1:, 1] < 0 )[0] + 1
        if peak == 'bounce':
            ## FIND the bounce time
            bounce_index = np.argmax(GWs[:, 0] >= 0)
            ## FIND the min in the 1.5 ms after the bounce
            index_after_bounce = np.argmax(GWs[:,0] >= u.convert_to_s(min_time)) + 1
            x_min = np.argmin(GWs[bounce_index:index_after_bounce, 1]) + bounce_index
            ##FIND MAX AFTER the min in 1.5 ms
            index2_after_bounce = index_after_bounce = np.argmax(GWs[:,0] >= GWs[x_min, 0] + u.convert_to_s(max_time)) + 1
            
            if use_der:
                x_max = 0
                derivative = IDL_derivative(GWs[x_min:index2_after_bounce+1, 0], GWs[x_min:index2_after_bounce+1, 1])
                flex_points = np.where((derivative[1:] * derivative[:-1] < 0))[0] + x_min
            
                for i in range(flex_points.size):
                    
                    if GWs[flex_points[i], 1] > 0 and GWs[flex_points[i+1], 1] < 0:
                        x_max = flex_points[i]
                        break
            else:
                x_max = np.argmax(GWs[x_min:index2_after_bounce, 1]) + x_min
        elif peak == 'highest':
            ## CUT the GWS
            if interval[0] is not None:
                start_index = np.argmax(GWs[:, 0] >= u.convert_to_s(interval[0]))
                GWs = GWs[start_index:, :]
            else:
                start_index = None
            if interval[1] is not None:
                GWs = GWs[:np.argmax(GWs[:, 0] >= u.convert_to_s(interval[1])), :]
            ## FIND the peak
            x_max = np.argmax(GWs[:, 1])
            ## FIND the min
            min_index = np.argmax(GWs[:, 0] >= GWs[x_max, 0] - u.convert_to_s(max_time))
            x_min = np.argmin(GWs[min_index:np.argmax(GWs[:, 0] >= GWs[x_max, 0] + u.convert_to_s(max_time)), 1]) + min_index
            if start_index is not None:
                x_min += start_index
                x_max += start_index
        ## Find the beginning and end of the peak
        if x_max > x_min:
            zeros_end_index = np.argmax(zeros>x_max)
            end_index = zeros[zeros_end_index]
            start_index = zeros[np.argmax(zeros>x_min) - 4]
        elif x_max < x_min:
            zeros_end_index = np.argmax(zeros>x_min)
            end_index = zeros[zeros_end_index]
            start_index = zeros[np.argmax(zeros>x_max) - 4]
        return start_index, x_min, x_max, end_index

    def __GWs_fourier_transform(self, GWs, indices):
        """
        This method applies FFT to a small portion of the GW strain to find the domiunant frequency of a specific
        obsillation
        Returns
            positive frequency range
            $\tilde{h} * \sqrt{freq}$
        """
        dt = np.abs(GWs[1, 0] - GWs[0, 0])
        ## Cut the GWs signal
        strain = np.zeros(11000)
        ## Pad the strain with zeros to increase the resolution
        if (GWs[indices[0]:indices[-1], 1]).size < 11000:
            strain[11000 - (GWs[indices[0]:indices[-1], 1]).size:] = GWs[indices[0]:indices[-1], 1]
        else:
            strain = GWs[indices[0]:indices[-1], 1]
        #find the frequencies
        freq = fftfreq(strain.size, dt)
        #freq = np.where(freq < 0, 0, freq)
        dft = np.abs(fft(strain)) * np.sqrt(freq)
        dft = np.where(freq > 0, dft, 0)
        freq = np.where(freq < 0, 0, freq)
        return freq, dft

    def __GWs_frequency_peak_indices(self, frequency, htilde):
        """
        Method that find the indices of the first and second frequency peak on
        a fourier transformed GW strain
        Return:
            list containing: first peak index, second peak index 
        """
        dhtilde_df = IDL_derivative(frequency, htilde)
        sign = dhtilde_df[1:] * dhtilde_df[:-1]
        ## FInd the peak frequency of the strain
        peak_frequency_index = np.argmax(htilde)


        ## find derivative sign change
        indices = np.argwhere( sign < 0)[:, 0]
        second_peak_index = 0
        for index in indices:
            if ( htilde[index] > htilde[second_peak_index] and \
                htilde[index] < htilde[peak_frequency_index] and \
                np.abs( frequency[index] - frequency[peak_frequency_index] ) > 10 and \
                frequency[index] > 200 ):
                second_peak_index = index
        return [peak_frequency_index, second_peak_index]
    
    #misc from h5 file: grav pot, time
    def grav_pot(self, data_h5):
        data = np.array(data_h5['gravpot']['data'])
        if data.ndim == 4:
            data = data[..., 0]
        return self.ghost.remove_ghost_cells(np.squeeze(data), self.dim)
    
    def gravitational_energy(self, data_h5):
        return 0.5 * self.rho(data_h5) * self.grav_pot(data_h5)

    def time(self, data_h5, correct_for_tob=False):
        if correct_for_tob:
            return np.array(data_h5['Parameters']['t'])-self.time_of_bounce_rho()
        return np.array(data_h5['Parameters']['t'])

    def file_list_hdf(self):
        """
        Listo of all the 'timestep' files in the outp-hdf folder.
        """
        file_list = os.listdir(self.hdf_path)
        #remove x00 files
        file_list = [x for x in file_list if x.startswith('h')]
        #return the sorted list of files
        file_list.sort()
        return file_list
    
    def find_file_from_time(self, time_to_find, tob=False, time_in_ms=True, return_index=False):
        """
        time in ms
        """
        if time_in_ms:
            time_to_find = u.convert_to_s(time_to_find)
        if tob is not False:
            time_to_find += self.time_of_bounce_rho()
        file_list = self.file_list_hdf()
        for (file, file_index) in zip(file_list, range(len(file_list))):
            data_file = self.open_h5(file)
            time = self.time(data_file)
            if time>=time_to_find:
                self.close_h5(data_file)
                if return_index:
                    return file, file_index
                return file
            self.close_h5(data_file)
        return None

    #time points computation: time of bounce and of BH
    def time_of_bounce(self):
        warnings.warn("WARNING: Works fine in 1D, may give incorrect results in higher dimensions")
        """
        Time of bounce computed with the entropy profiles.
        BEWARE Does not work that well
        """
        file_list = self.file_list_hdf()
        mass_thresh = 1.01222
        entr_thresh = 3.0
        for file in file_list:
            data_h5 = self.open_h5(file)
            s = self.entropy(data_h5)
            M = self.mass_shells(data_h5)
            mass_index = np.argmax(M>=mass_thresh)
            entropy_index = np.unravel_index(np.argmax(s >= entr_thresh), s.shape)
            if entropy_index[-1] != 0 and entropy_index[-1] < mass_index:
                t = self.time(data_h5)
                self.close_h5(data_h5)
                return t
            self.close_h5(data_h5)
           
    def time_of_bounce_rho(self):
        """
        Empirical criterion: time of bounce defined as the time at which
        the central density (max desity) raises above 2.5e14 g/cm^3 before a rapid fall.
        If we do not reach that density we lower the threshold to 2
        """
        rho_data = self.rho_max(False)
        rho_index = np.argmax(rho_data[:,1] > 1.4e14)
        if rho_index == 0 or rho_data[rho_index, 0] >= 0.6:
            rho_index = np.argmax(rho_data[:,1]>2e14)
        return rho_data[rho_index, 0]

    def time_of_BH(self):
        """
        Attempt to find when and if a BH is formed.
        For now is when the rho max raises above a certain limit (6e15 g/cm3)
        """
        rho_data = self.rho_max()
        rho_BH = 6e15
        return rho_data[np.argmax(rho_data[:,1]>=rho_BH),0]

    #gravitational potential
    def central_grav_pot_at_time(self,data_h5):
        if self.dim == 1:
            self.grav_pot(data_h5)[0]
        if self.dim == 2:
            self.grav_pot(data_h5)[0, 0]
        return self.grav_pot(data_h5)[0, 0, 0]

    def get_central_grav_pot(self):
        """
        Evolution of the gravitational potential
        indices
        1: time
        2: central gravitational potential
        """
        file_list = self.file_list_hdf()
        data_out = np.zeros((len(file_list), 2))
        for (file, index) in zip(file_list, range(len(file_list))):
            data_file = self.open_h5(file)
            data_out[index, 0] = self.time(data_file)
            data_out[index, 1] = self.central_grav_pot_at_time(data_file)
            self.close_h5(data_file)
        return data_out

    #integrated quantities
    def integrated_neutrino_lum(self, corrected_by_tob = True):
        """
        indices
        1: time
        2: luminosity flux nue  5: number luminosity flux nue
        3: luminosity flux nua  6: number luminosity flux nua
        4: luminosity flux nux  7: number luminosity flux nux
        """
        nu_tmp = load_file(self.log_path, self.integrated_nu)
        if corrected_by_tob:
            nu_tmp[:, 2] -= self.time_of_bounce_rho()
        return np.stack((nu_tmp[:, 2], nu_tmp[:, 38], nu_tmp[:, 39],
                         0.25 * nu_tmp[:, 40], nu_tmp[:, 35],
                         nu_tmp[:, 36], 0.25 * nu_tmp[:, 37]), axis=1)
        
    def neutrino_mean_energies_grey(self):
        """
        indices
        1: time
        2: mean energy nue  
        3: mean energy nua  
        4: mean energy nux  
        """
        nu = self.integrated_neutrino_lum()
        ene_mean = np.zeros((nu.shape[0],4))
        ene_mean[:,0] = nu[:,0]
        ene_mean[:,1] = u.convert_to_MeV(nu[:,1]/nu[:,4])
        ene_mean[:,2] = u.convert_to_MeV(nu[:,2]/nu[:,5])
        ene_mean[:,3] = u.convert_to_MeV(nu[:,3]/nu[:,6])
        return np.stack((nu[:, 0], u.convert_to_MeV(nu[:, 1]/nu[:, 4]),
                         u.convert_to_MeV(nu[:, 2]/nu[:, 5]),
                         u.convert_to_MeV(nu[:, 3]/nu[:, 6])), axis = 1)

    def rho_max(self, correct_for_tob=True):
        """
        indices
        1: time
        2: rho max
        """
        rho = load_file(self.log_path, self.rho_max_file)
        if correct_for_tob:
            rho[:,2] -= self.time_of_bounce_rho()
        return np.stack((rho[:, 2], rho[:, 3]), axis = 1)
    
    def tot_mass(self, correct_for_tob=True):
        """
        indices
        1: time
        2: total mass
        """
        rho = load_file(self.log_path, self.rho_max_file)
        if correct_for_tob:
            rho[:,2] -= self.time_of_bounce_rho()
        return np.stack((rho[:, 2], u.convert_to_solar_masses(rho[:, 4])), axis = 1)

    def Yl_cent(self, correcetd_by_tob = True):
        """
        indices
        1: time
        2: Y_l cent
        """
        Y_l = load_file(self.log_path, self.rho_max_file)
        if correcetd_by_tob:
            Y_l[:,2] -= self.time_of_bounce_rho()
        return np.stack((Y_l[:,2], Y_l[:,9]), axis = 1)

    def Ye_max_min_cent(self, correcetd_by_tob = True):
        """
        indices
        1: time
        2: Ye max
        3: Ye min
        4: Ye cent
        """
        Ye = load_file(self.log_path, self.rho_max_file)
        if correcetd_by_tob:
            Ye[:,2] -= self.time_of_bounce_rho()
        return np.stack((Ye[:,2], Ye[:,6], Ye[:,7], Ye[:,8]), axis = 1)

    def rotational_energy_integrated(self, corrected_by_tob = True):
        """
        1: time
        2: total rotational energy
        """
        en = load_file(self.log_path, 'mag.dat')
        if corrected_by_tob:
            en[:,2] -= self.time_of_bounce_rho()
        return np.stack((en[:, 2], en[:, 3]), axis = 1)

    def BV_frequency(self, data_h5):
        """
        Returns the Brunt-Väisälä frequency for a specific timestep
        """
        
        omega_BV = ( 1 / self.speed_of_sound(data_h5) ** 2 * IDL_derivative( self.cell.radius(self.ghost), self.gas_pressure(data_h5) ) - \
            IDL_derivative( self.cell.radius(self.ghost), self.rho(data_h5) ) ) * IDL_derivative( self.cell.radius(self.ghost), self.grav_pot(data_h5) ) / self.rho(data_h5)
        return omega_BV
    
    def BV_frequency_profile(self, tob_corrected = True):
        """
        Calculates the Brunt-Väisälä frequency and returns three arrays
        time: array of time steps
        radius: array of radial distance from tne center
        omega_BV: array of len(radius) x len(time) containing the angular averaged BV frequency 
        """
        BV_file = os.path.join(self.storage_path, 'BV_frequency.h5')
        if not os.path.exists(BV_file):
            warnings.warn("BV frequency file not found. Creating one...\n" + \
                          "Please wait...")
            time, radius, omega_BV = self.__BV_frequency_profile()
            file_h5 = h5py.File(BV_file, 'w')
            file_h5.create_dataset("time", data=time)
            file_h5.create_dataset("radius", data=radius)
            file_h5.create_dataset("BV_frequency", data=omega_BV)
            file_h5.create_dataset("ToB_corrected", data=True)
            file_h5.close()
        data = self.open_h5(BV_file)
        time = data["time"][...]
        omega_BV = data["BV_frequency"][...]
        radius = data["radius"][...]
        if not tob_corrected:
            time += self.time_of_bounce_rho()
        return time, radius, omega_BV

    def __BV_frequency_profile(self):
        """
        Returns the Brunt-Väisälä frequency radial profile, using the angular averaged quantities
        """
        file_list_hdf = self.file_list_hdf()
        radius = self.cell.radius(self.ghost)
        indices = np.arange(len(file_list_hdf))
        omega_BV = np.zeros((len(radius), len(file_list_hdf)))
        time = np.zeros(len(file_list_hdf))
        for (file, i) in zip(file_list_hdf, indices):
            data_h5 = self.open_h5(file)
            time[i] = self.time(data_h5)
            if self.dim == 1:
               omega_BV[:, i] = self.BV_frequency(data_h5)
            else:
                cs2 = self.speed_of_sound(data_h5).mean(axis = tuple(range(self.dim - 1))) ** 2
                pgas =  self.gas_pressure(data_h5).mean(axis = tuple(range(self.dim - 1)))
                rho = self.rho(data_h5).mean(axis = tuple(range(self.dim - 1)))
                phi = self.rho(data_h5).mean(axis = tuple(range(self.dim - 1)))
                omega_BV[:, i] = ( 1 / cs2 * IDL_derivative( radius, pgas ) - IDL_derivative( radius, rho )) * IDL_derivative( radius, phi ) / rho
        return time - self.time_of_bounce_rho(), radius, omega_BV 

    #derivatives
    def velocity_radial_derivative(self, data_h5):
        radius = self.cell.radius(self.ghost)
        velocity = self.radial_velocity(data_h5)
        return IDL_derivative(radius, velocity) * radius / velocity

    def entropy_radial_derivative(self, data_h5):
        radius = self.cell.radius(self.ghost)
        entropy = self.entropy(data_h5)
        return IDL_derivative(radius, entropy) * radius / entropy
    
    def omega_radial_derivative(self, data_h5):
        omega = self.omega(data_h5)
        radius = self.cell.radius(self.ghost)
        return IDL_derivative(radius, omega)
    
    def innercore_radius_single(self, data_h5):
        """
        We define the inner core of a star as the region in sonic contact with the
        centre. This is the region where the velocity of the fluid is lower than the
        speed of sound.
        """
        cs2 = self.speed_of_sound(data_h5) ** 2
        velocity = self.radial_velocity(data_h5) ** 2
        radius = self.cell.radius(self.ghost)
        return radius[np.argmax(velocity >= cs2, axis = -1)]
    
    def innercore_radius(self, tob_correction, save_name = 'innercore_radius', **kwargs):
        self.__save_radii('innercore', tob_correction, save_name, **kwargs)

    def get_innercore_radius(self, innercore_radius = False, indices = False, min_max_average = False,
                                ret_time = False, ghost_cells = False, tob_corrected = True):
        innercore_file = os.path.join(self.storage_path,'innercore_radius.h5')
        if not os.path.exists(innercore_file):
            warnings.warn("Innercore radius file not found. Creating one with default settings.\n" + \
                          "If different settings are needed please refer to the " + \
                          "\"Innercore_radius(...)\" method")
            g = self.ghost_cells - 1
            self.innercore_radius(tob_corrected, t_l = g, t_r = g, p_l = g, p_r = g)
        return self.__get_radii(innercore_file, innercore_radius, indices, min_max_average,
                                ret_time, ghost_cells, tob_corrected)
       
    def PNS_radius_single(self, data_h5):
        """
        A bit empirical but should work fine. We consider everything with
        a density >= 10^11 g/cm^3 as PNS.
        """
        radius = self.cell.radius(self.ghost)
        rho_thres = 1e11
        rho = self.rho(data_h5)
        indices = np.argmax(rho <= rho_thres, axis = -1)
        return radius[indices]

    def PNS_radius(self, tob_correction, save_name = 'PNS_radius', **kwargs):
        self.__save_radii('PNS', tob_correction, save_name, **kwargs)
        
    def get_PNS_radius(self, PNS_radius = False, indices = False, min_max_average = False,
                       ret_time = False, ghost_cells = False, tob_corrected = True):
        PNS_file = os.path.join(self.storage_path,'PNS_radius.h5')
        if not os.path.exists(PNS_file):
            warnings.warn("PNS radius file not found. Creating one with default settings.\n" + \
                          "If different settings are needed please refer to the " + \
                          "\"PNS_radius(...)\" method")
            g = self.ghost_cells - 1
            self.PNS_radius(tob_corrected, t_l = g, t_r = g, p_l = g, p_r = g)
        return self.__get_radii(PNS_file, PNS_radius, indices, min_max_average,
                                ret_time, ghost_cells, tob_corrected)
 
    def gain_radius_single(self, data_h5, PNS_radius):
        """
        We find the gain radius as the region outside the PNS where neutrino
        heating dominates over neutrino cooling.
        """
        radius = self.cell.radius(self.ghost)
        nu_heat = self.nu_heat(data_h5)
        r = radius.copy()
        while r.ndim != nu_heat.ndim:
            r = r[None, ...]
        nu_heat = np.where(r >= PNS_radius[..., None], nu_heat,0)
        gain_index = np.argmax(nu_heat > 0, axis = -1)
        return radius[gain_index]
    
    def gain_radius(self, tob_correction, save_name = 'gain_radius', **kwargs):
        self.__save_radii('gain', tob_correction, save_name, **kwargs)

    def get_gain_radius(self, gain_radius = False, indices = False, min_max_average = False,
                        ret_time = False, ghost_cells = False, tob_corrected = True):
        gain_file = os.path.join(self.storage_path,'gain_radius.h5')
        if not os.path.exists(gain_file):
            warnings.warn("Gain radius file not found. Creating one with PNS radius settings.\n" + \
                          "If different settings are needed please refer to the " + \
                          "\"gain_radius(...)\" method")
            ghost_par = self.get_PNS_radius(ghost_cells = True)
            self.gain_radius(tob_corrected, **ghost_par)
        return self.__get_radii(gain_file, gain_radius, indices, min_max_average,
                                ret_time, ghost_cells, tob_corrected)

    def neutrino_sphere_radius_single(self, data_h5):
        """
        To calculate neutrino sphere we consider:
                            \sum_{ebin} (L_\nu(r, \theta, \phi, \nu_e) * \kappa(r, \theta, \phi, \nu_e)
        k(r, theta, phi) =  ---------------------------------------------------------------------------
                                            \sum_{ebin} (L_\nu(r, \theta, \phi, \nu_e) 
        where \nu is the neutrino flavour. Then the neutrino sphere is the R extrema of the integral 
        tau = int^R_\infty dr k(r, \theta, \phi)
        """
        tau = 1.0
        neutrino_fluxes = self.neutrino_fluxes(data_h5)
        neutrino_opacities = self.neutrino_opacity_flux(data_h5)
        dr = self.cell.dr(self.ghost)[..., None]
        if self.dim > 1:
            neutrino_data = np.sum(np.sum(neutrino_fluxes*neutrino_opacities, axis = self.dim), axis=-1) / \
                            np.sum(np.sum(neutrino_fluxes, axis = self.dim), axis=-1)
        else:
            neutrino_data = np.sum(neutrino_fluxes*neutrino_opacities, axis = self.dim) / \
                            np.sum(neutrino_fluxes, axis = self.dim)
        np.nan_to_num(neutrino_data, False, 0)
        while (neutrino_data.ndim) != dr.ndim:
            dr = dr[None, ...]
        neutrino_data = np.flip(neutrino_data * dr, axis = -2)
        ind = np.argmax(np.cumsum(neutrino_data, axis = -2) >= tau, axis = -2)
        radius = np.flip(self.cell.radius(self.ghost))
        if type(ind) == int:
            ind = np.zeros()
        #deallocate arrays
        del neutrino_fluxes
        del neutrino_opacities
        del dr
        del neutrino_data
        return radius[ind]
    
    def neutrino_sphere_radius(self, tob_correction, save_name = 'neutrino_spheres_radius', **kwargs):
        self.__save_radii('neutrino', tob_correction, save_name, **kwargs)

    def get_neutrino_sphere_radius(self, radius_nue = False, radius_nua = False, radius_nux = False,
                        indices = False, min_max_average_nue = False, min_max_average_nua = False,
                        min_max_average_nux = False, ret_time = False, ghost_cells = False, 
                        tob_corrected = True):
        nu_file = os.path.join(self.storage_path,'neutrino_spheres_radius.h5')
        if not os.path.exists(nu_file):
            warnings.warn("Neutrino sphere radius file not found. Creating one with default settings.\n" + \
                          "If different settings are needed please refer to the " + \
                          "\"neutrino_spere_radius(...)\" method")
            g = self.ghost_cells - 1
            self.neutrino_sphere_radius(tob_corrected, t_l = g, t_r = g, p_l = g, p_r = g)
        data = self.open_h5(nu_file)
        time = np.array(data['time'])
        if (not tob_corrected) and np.array(data['tob_correction']):
            time += self.time_of_bounce_rho()
        elif tob_corrected and not np.array(data['tob_correction']):
            time -= self.time_of_bounce_rho()
        quantities_to_return = []
        if radius_nue:
            quantities_to_return.append(np.array(data['radii']['nue']))
        if radius_nua:
            quantities_to_return.append(np.array(data['radii']['nua']))
        if radius_nux:
            quantities_to_return.append(np.array(data['radii']['nux']))
        if indices:
            quantities_to_return.append(np.array(data['indices']))
        if min_max_average_nue:
            quantities_to_return.append(np.stack((time, np.array(data['min']['nue']),
                                        np.array(data['max']['nue']), np.array(data['average']['nue'])),
                                        axis = -1))
        if min_max_average_nua:
            quantities_to_return.append(np.stack((time, np.array(data['min']['nua']),
                                        np.array(data['max']['nua']), np.array(data['average']['nua'])),
                                        axis = -1))
        if min_max_average_nux:
            quantities_to_return.append(np.stack((time, np.array(data['min']['nux']),
                                        np.array(data['max']['nux']), np.array(data['average']['nux'])),
                                        axis = -1))
        if ret_time:
            quantities_to_return.append(time)
        if ghost_cells:
            g_cells = {'p_l': list(data['ghost_cells']['phi'])[0],
                       'p_r':  list(data['ghost_cells']['phi'])[1],
                       't_l':  list(data['ghost_cells']['theta'])[0],
                       't_r':  list(data['ghost_cells']['theta'])[1],
                       'r_l':  list(data['ghost_cells']['radius'])[0],
                       'r_r':  list(data['ghost_cells']['radius'])[1]}
            quantities_to_return.append(g_cells)
        self.close_h5(data)
        if len(quantities_to_return) == 1:
            return quantities_to_return[0]
        if len(quantities_to_return) == 0:
            print("Nothing to return :(")
            return None
        return quantities_to_return

    def shock_radius_single(self, data_h5, zero = False):
        """
        This method to calculate the shock propagation is still under development, however
        for now it works acceptably.
        This method consists of three steps
        1 - Find where the entropy is above 6 kb/br and the radial velocity above -3e9 cm/s, while
            their radial derivatives are below -0.3 and -1.5. The search is conducted from the infinity
            to the simulation centre. We choose entropy because is a good indicator of matter that is moving
            and radial velocity because shock mostrly propagates radially.
        2 - Repeat search in point one neglecting derivatives thresholds and performing the search from
            centre to infinity, with the inequalities direction changed.
            If search in point 1 fails at some angles we replace those values with point 2.
        3 - If any point is above or below <R_sh> \pm 1.5\sigma then we replace it with a value got from
            cubic interpolation of the other values (only in 2D, 3D coming soon). This should get rid of
            some spikes that are not compatible with visual comparison.
        As you can see this is not perfect. If someone has better ideas, please let me know.
        """
        velocity = np.flip(self.radial_velocity(data_h5), axis = -1)
        if zero:
            return np.zeros(velocity.shape[:-1])
        r = self.cell.radius(self.ghost)
        velocity_derivative = np.flip(self.velocity_radial_derivative(data_h5), axis = -1)
        entropy = np.flip(self.entropy(data_h5), axis = -1)
        entropy_derivative = np.flip(self.entropy_radial_derivative(data_h5), axis = -1)
        #thresholds
        vel_der_tresh = -1.5
        vel_tresh = -3e9
        entr_tresh = 6
        entr_der_tresh = -0.3
        if self.dim == 1:
            ind = velocity.shape[-1] - 1 - np.argmax((np.flip(velocity, axis = -1) <= vel_tresh) & \
                                        (np.flip(entropy, axis = -1) <= entr_tresh), axis = -1)
            shock = np.flip(r)[ind]
        else:
            ind = np.argmax((velocity_derivative <= vel_der_tresh) & (entropy_derivative <= entr_der_tresh) & \
                        (velocity >= vel_tresh) & (entropy >= entr_tresh), axis = -1)
            r = self.cell.radius(self.ghost)
            ind2 = velocity.shape[-1] - 1 - np.argmax((np.flip(velocity, axis = -1) <= vel_tresh) & \
                                            (np.flip(entropy, axis = -1) <= entr_tresh), axis = -1)
            oob = ind2 ==  velocity.shape[-1] - 1
            if oob.any():
                ind2 = np.where(oob, ind, ind2)
            oob = ind == 0
            if oob.any():
                ind = np.where(oob, ind2, ind)
            shock = np.flip(r)[ind]
            mean_ind = shock != r[-1]
            oob = shock >= (shock.mean(where = mean_ind) + 1.5 * np.std(shock, where = mean_ind))
            if oob.any():
                shock2 = np.flip(r)[ind2]
                shock = np.where(oob, shock2, shock)
                mean_ind = shock != r[-1]
            sh = shock.mean(where = mean_ind)
            sh_std = np.std(shock, where = mean_ind)
            k = 1.5
            out_value = sh - k * sh_std
            while (out_value <= 0 and k > 0):
                out_value = sh - k * sh_std
                k -= 0.5
            oob = shock < (out_value)
            if oob.any() and self.dim == 2:
                theta = self.cell.theta(self.ghost)
                ind_to_interp = np.invert(oob)
                f = interp1d(theta[ind_to_interp], shock[ind_to_interp], kind = 'linear', fill_value="extrapolate")
                shock2 = f(theta)
                shock = np.where((oob) & (shock2 > 0), shock2, shock)
            if (shock < 0).any():
                raise TypeError("Shock radius is negative!!!!!!")
        return shock

    def shock_radius(self, tob_correction, save_name = 'shock_radius', **kwargs):
        self.__save_radii('shock', tob_correction, save_name, **kwargs)

    def get_shock_radius(self, shock_radius = False, indices = False, min_max_average = False,
                        ret_time = False, ghost_cells = False, tob_corrected = True):
        shock_file = os.path.join(self.storage_path,'shock_radius.h5')
        if not os.path.exists(shock_file):
            warnings.warn("Shock radius file not found. Creating one with default settings.\n" + \
                          "If different settings are needed please refer to the " + \
                          "\"shock_radius(...)\" method")
            g = self.ghost_cells - 1
            self.shock_radius(tob_corrected, t_l = g, t_r = g, p_l = g, p_r = g)
        return self.__get_radii(shock_file, shock_radius, indices, min_max_average,
                                ret_time, ghost_cells, tob_corrected)

    def __save_radii(self, radius_to_calculate: Literal['PNS', 'gain', 'shock', 'neutrino', 'innercore'], tob_correction, 
                save_name, **kwargs):
        methods = {'PNS': self.PNS_radius_single,
                   'gain': self.gain_radius_single,
                   'neutrino': self.neutrino_sphere_radius_single, 
                   'shock': self.shock_radius_single,
                   'innercore': self.innercore_radius_single}
        save_path = check_path(self.storage_path, save_name)
        self.ghost.update_ghost_cells(**kwargs)
        file_list = self.file_list_hdf()
        if radius_to_calculate == 'neutrino':
            data = np.zeros((len(file_list), 10))
        else:
            data = np.zeros((len(file_list), 4))
        indices = np.arange(len(file_list))
        calculate_radius = methods[radius_to_calculate]
        tob = self.time_of_bounce_rho() 
        if radius_to_calculate == 'gain':
            PNS_radius = self.get_PNS_radius(PNS_radius = True)
        for (f, i) in zip(file_list, indices):
            no_neutrino_error = False #in some simulations neutrino output is not saved 
                                      #for all timesteps to save storage space
            data_h5 = self.open_h5(f)
            data[i, 0] = self.time(data_h5)
            if radius_to_calculate == 'gain':
                radius_single = calculate_radius(data_h5, PNS_radius[..., i])
            elif radius_to_calculate == 'shock' and data[i, 0] <= tob:
                radius_single = calculate_radius(data_h5, True)
            elif radius_to_calculate == 'neutrino':
                try:
                    radius_single = calculate_radius(data_h5)
                except:
                    no_neutrino_error = True
            else:
                radius_single = calculate_radius(data_h5)
            
            if i == 0:
                rad_out = radius_single[..., None]
            else:
                if no_neutrino_error:
                    if i == 0:
                        raise ValueError("No neutrino saved in this simulation.")
                    radius_single = rad_out[..., i-1]
                rad_out = np.concatenate((rad_out, radius_single[..., None]),
                                        axis = -1)
            if radius_to_calculate == 'neutrino':
                if no_neutrino_error:
                    data[i, 1:] = data[i-1, 1:]
                else:
                    for nu in range(3):
                        rr = radius_single[..., nu]
                        data[i, 3 * nu + 1] = self.ghost.remove_ghost_cells_radii(rr,
                                                                                self.dim).min()
                        data[i, 3 * nu + 2] = self.ghost.remove_ghost_cells_radii(rr,
                                                                                self.dim).max()
                        data[i, 3 * nu + 3] = self.ghost.remove_ghost_cells_radii(rr,
                                                                                self.dim).mean()
            else:
                data[i, 1] = self.ghost.remove_ghost_cells_radii(radius_single, self.dim).min()
                data[i, 2] = self.ghost.remove_ghost_cells_radii(radius_single, self.dim).max()
                data[i, 3] = self.ghost.remove_ghost_cells_radii(radius_single, self.dim).mean()
            self.close_h5(data_h5)
        if tob_correction:
            data[:, 0] -= tob
        save_h5(save_path, rad_out, data,
                indices, self.ghost.return_ghost_dictionary(), tob_correction)
        self.ghost.restore_default()
    
    def __get_radii(self, file, radius = False, indices = False, min_max_average = False,
                    ret_time = False,  ghost_cells = False, tob_corrected = False):
        data = self.open_h5(file)
        time = np.array(data['time'])
        if (not tob_corrected) and np.array(data['tob_correction']):
            time += self.time_of_bounce_rho()
        elif tob_corrected and not np.array(data['tob_correction']):
            time -= self.time_of_bounce_rho()
        quantities_to_return = []
        if radius:
            quantities_to_return.append(np.array(data['radii']))
        if indices:
            quantities_to_return.append(np.array(data['indices']))
        if min_max_average:
            quantities_to_return.append(np.stack((time, np.array(data['min']), 
                                                  np.array(data['max']), np.array(data['average'])),
                                        axis = -1))
        if ret_time:
            quantities_to_return.append(time)
        if ghost_cells:
            g_cells = {'p_l': list(data['ghost_cells']['phi'])[0],
                       'p_r':  list(data['ghost_cells']['phi'])[1],
                       't_l':  list(data['ghost_cells']['theta'])[0],
                       't_r':  list(data['ghost_cells']['theta'])[1],
                       'r_l':  list(data['ghost_cells']['radius'])[0],
                       'r_r':  list(data['ghost_cells']['radius'])[1]}
            quantities_to_return.append(g_cells)
        self.close_h5(data)
        if len(quantities_to_return) == 1:
            return quantities_to_return[0]
        if len(quantities_to_return) == 0:
            print("Nothing to return :(")
            return None
        return quantities_to_return

    #Masses
    def mass_shells(self, data_h5):
        """
        It returns the mass contained inside each radius.
        """
        rho = self.rho(data_h5) * self.cell.dVolume_integration(self.ghost)
        rho = rho.sum(axis=tuple(range(rho.ndim - 1)))    
        return u.convert_to_solar_masses(np.cumsum(rho))
    
    def get_masses_and_energies(self, time = False, PNS_mass = False, PNS_mag_ene = False,
                        PNS_rot_ene = False, PNS_pol_kin_ene = False, PNS_conv_ene= False, 
                        ejected_mass = False, explosion_energy = False, gain_mass = False, 
                        gain_heating = False, mass_accretion_500km = False, tob_corrected = True, 
                        name = 'mass_energy'):
        path = os.path.join(self.storage_path, name + '.h5')
        if not os.path.exists(path):
            warnings.warn(name + " file does not exists. Creating one now.")
            self.__save_energies_and_masses(path, tob_corrected)
        data = self.open_h5(path)
        quantities_list = []
        if time:
            if (tob_corrected and data['tob_correction']) or \
                not (tob_corrected and data['tob_correction']):
                quantities_list.append(np.array(data['time']))
            elif tob_corrected and not data['tob_correction']:
                quantities_list.append(np.array(data['time']) - self.time_of_bounce_rho())
            else:
                quantities_list.append(np.array(data['time']) + self.time_of_bounce_rho())
        if PNS_mass:
            quantities_list.append(np.array(data['PNS_mass']))
        if PNS_mag_ene:
            quantities_list.append(np.array(data['PNS_mag_ene']))
        if PNS_rot_ene:
            quantities_list.append(np.array(data['PNS_rot_ene']))
        if PNS_pol_kin_ene:
            quantities_list.append(np.array(data['PNS_pol_kin_ene']))
        if PNS_conv_ene:
            quantities_list.append(np.array(data['PNS_convec_ene']))
        if ejected_mass:
            quantities_list.append(np.array(data['ejected_mass']))
        if explosion_energy:
            quantities_list.append(np.array(data['explosion_ene']))
        if gain_mass:
            quantities_list.append(np.array(data['gain_mass']))
        if gain_heating:
            quantities_list.append(np.array(data['gain_heating']))
        if mass_accretion_500km:
            quantities_list.append(np.array(data['mass_accr_500km'])) 
        self.close_h5(data)
        if len(quantities_list) == 1:
            return quantities_list[0]
        if len(quantities_list) == 0:
            print("Nothing to return :(")
            return None
        return quantities_list
    
    def get_innercore_mass_enregies(self, time = False, mass = False, kin_rot_ene = False,
                                    kin_tot_ene = False, mag_ene = False, grav_ene = False,
                                    T_W = False, tob_corrected = True):
        innercore_file = os.path.join(self.storage_path,'innercore_mass_energies.h5')
        if not os.path.exists(innercore_file):
            warnings.warn("Innercore quantities file not found. Creating one with default settings.\n" + \
                          "If different settings are needed please refer to the " + \
                          "\"shock_radius(...)\" method")
            self.__save_innercore_mass_energies()
        data = h5py.File(innercore_file, 'r')
        quantities_list = []
        if time:
            if (tob_corrected and data['tob_correction']) or \
                not (tob_corrected and data['tob_correction']):
                quantities_list.append(np.array(data['time']))
            elif tob_corrected and not data['tob_correction']:
                quantities_list.append(np.array(data['time']) - self.time_of_bounce_rho())
            else:
                quantities_list.append(np.array(data['time']) + self.time_of_bounce_rho())
        if mass:
            quantities_list.append(np.array(data['mass']))
        if kin_rot_ene:
            quantities_list.append(np.array(data['kin_rot_ene']))
        if kin_tot_ene:
            quantities_list.append(np.array(data['kin_tot_ene']))
        if mag_ene:
            quantities_list.append(np.array(data['mag_ene']))
        if grav_ene:
            quantities_list.append(np.array(data['grav_ene']))
        if T_W:
            quantities_list.append(np.array(data['T_W']))
        self.close_h5(data)
        if len(quantities_list) == 1:
            return quantities_list[0]
        if len(quantities_list) == 0:
            print("Nothing to return :(")
            return None
        return quantities_list
        
    
    def __save_innercore_mass_energies(self):
        rad_in, indices, time, g_cells = self.get_innercore_radius(innercore_radius = True,
                                        indices = True, min_max_average = False, ret_time = True,
                                        ghost_cells = True, tob_corrected = True)
        data_out = np.zeros((time.size, 6))
        data_out[:, 0] = time
        file_list = self.file_list_hdf()
        dV = self.cell.dVolume_integration(self.ghost)
        radius = self.cell.radius(self.ghost)
        for file, index in zip(file_list, indices):
            ## Quantities calculation
            data_h5 = self.open_h5(file)
            rho = self.rho(data_h5)
            rot_ene = 0.5 * rho * self.phi_velocity(data_h5) ** 2
            kin_ene = 0.5 * (self.theta_velocity(data_h5) + self.radial_velocity(data_h5)) ** 2 +\
                    rot_ene
            ene_mag = self.magnetic_energy_per_unit_volume(data_h5)
            ene_grav = 0.5 * self.grav_pot(data_h5) * rho
            self.close_h5(data_h5)
            dVolume = np.where(radius <= self.ghost.remove_ghost_cells_radii(rad_in[..., index],
                                        self.dim, **g_cells)[..., None], dV, 0)
            data_out[index, 1] = u.convert_to_solar_masses(np.sum(rho * dVolume))
            data_out[index, 2] = np.sum(rot_ene * dVolume)
            data_out[index, 3] = np.sum(kin_ene * dVolume)
            data_out[index, 4] = np.sum(ene_mag * dVolume)
            data_out[index, 5] = np.sum(ene_grav * dVolume)
        file_out = h5py.File(os.path.join(self.storage_path, 'innercore_mass_energies.h5'), 'w')
        file_out.create_dataset('time', data = data_out[:, 0])
        file_out.create_dataset('mass', data = data_out[:, 1])
        file_out.create_dataset('kin_rot_ene', data = data_out[:, 2])
        file_out.create_dataset('kin_tot_ene', data = data_out[:, 3])
        file_out.create_dataset('mag_ene', data = data_out[:, 4])
        file_out.create_dataset('grav_ene', data = data_out[:, 5])
        file_out.create_dataset('T_W', data = np.abs(data_out[:, 2] / data_out[:, 5]))
        file_out.create_dataset('tob_correction', data = True)
        file_out.close()

    def __save_energies_and_masses(self, save_path, tob_corrected):
        PNS = self.__PNS_energies_and_mass(tob_corrected)
        mass_accretion = self.__mass_accretion(500, True)
        unbound = self.__unbound_mass_energy()
        gain = self.__gain_mass_and_heating()
        mass_trajectories = self.__fluid_trajectories()
        file_out = h5py.File(save_path, 'w')
        file_out.create_dataset('time', data = PNS[:, 0])
        file_out.create_dataset('tob_correction', data = tob_corrected)
        file_out.create_dataset('PNS_mass', data = PNS[:, 1])
        file_out.create_dataset('PNS_pol_kin_ene', data = PNS[:, 2])
        file_out.create_dataset('PNS_rot_ene', data = PNS[:, 3])
        file_out.create_dataset('PNS_mag_ene', data = PNS[:, 4])
        file_out.create_dataset('PNS_convec_ene', data = PNS[:, 5])
        file_out.create_dataset('mass_accr_500km', data = mass_accretion)
        file_out.create_dataset('ejected_mass', data = unbound[:, 1])
        file_out.create_dataset('explosion_ene', data = unbound[:, 0])
        file_out.create_dataset('gain_mass', data = gain[:, 0])
        file_out.create_dataset('gain_heating', data = gain[:, 1])
        file_out.create_dataset('fluid_trajectories', data = mass_trajectories)
        self.close_h5(file_out)

    #masses and energies 
    def __gain_mass_and_heating(self):
        dV = self.cell.dVolume_integration(self.ghost)
        file_list = self.file_list_hdf()
        shock_radius, indices, ghost_cells_shock = self.get_shock_radius(shock_radius = True, 
                                        indices = True, ghost_cells = True)
        gain_radius, ghost_cells_gain = self.get_gain_radius(gain_radius = True, 
                                        ghost_cells = True)
        r = self.cell.radius(self.ghost)
        while r.ndim != dV.ndim:
            r = r[None, ...]
        data_out = np.zeros((len(file_list), 2))
        for (file, index) in zip(file_list, indices):
            s_r = self.ghost.remove_ghost_cells_radii(shock_radius[..., index],
                                                      self.dim, **ghost_cells_shock)
            g_r = self.ghost.remove_ghost_cells_radii(gain_radius[..., index],
                                                      self.dim, **ghost_cells_gain)
            up_lim = np.maximum(s_r, g_r)
            low_lim = np.minimum(s_r, g_r)
            data_h5 = self.open_h5(file)
            rho = self.rho(data_h5) * dV
            nu_heat = self.nu_heat(data_h5) * dV
            summ_index = np.where((r >= low_lim[..., None]) & (r <= up_lim[..., None]), 1, 0)
            data_out[index, 0] = (rho * summ_index).sum()
            data_out[index, 1] = (nu_heat * summ_index).sum()
            self.close_h5(data_h5)
        data_out[:, 0] = u.convert_to_solar_masses(data_out[:, 0])
        return data_out

    def __mass_accretion(self, distance, distance_in_km=False):
        if distance_in_km:
            distance = u.convert_to_cm(distance)
        dV = self.cell.dVolume_integration(self.ghost) / self.cell.dr(self.ghost)
        r_index = np.argmax(self.cell.radius(self.ghost) >= distance)
        file_list = self.file_list_hdf()
        data_out = np.zeros((len(file_list)))
        for (f, index) in zip(file_list, range(len(file_list))):
            data_h5 = self.open_h5(f)
            data_out[index] = np.sum(dV[..., r_index] * self.rho(data_h5)[..., r_index] * \
                                        self.radial_velocity(data_h5)[..., r_index])
            self.close_h5(data_h5)
        return -u.convert_to_solar_masses(data_out)
      
    def __PNS_energies_and_mass(self, tob_corrected):
        dV = self.cell.dVolume_integration(self.ghost)
        file_list = self.file_list_hdf()
        data_out = np.zeros((len(file_list), 6))
        PNS_radius, indices, data_out[:, 0], ghost_cells = self.get_PNS_radius(PNS_radius = True,
                                           indices = True, ret_time = True, ghost_cells = True,
                                           tob_corrected = tob_corrected)
        r = self.cell.radius(self.ghost)
        while r.ndim != dV.ndim:
            r = r[None, ...]
        for (index, file) in zip(indices, file_list):
            data_h5 = self.open_h5(file)
            PNS = np.where(r <= self.ghost.remove_ghost_cells_radii(PNS_radius[..., index], 
                           self.dim, **ghost_cells)[..., None], 1, 0)
            rho = self.rho(data_h5) * PNS * dV
            if self.dim > 1:
                rot_energy = rho * self.phi_velocity(data_h5) ** 2
                conv_energy = rho * self.theta_velocity(data_h5) ** 2
                kin_energy = rho * (self.radial_velocity(data_h5) ** 2 \
                            + self.theta_velocity(data_h5) ** 2)
                mag_energy = self.magnetic_energy_per_unit_volume(data_h5) * PNS * dV
            else:
                rot_energy = np.zeros(rho.shape)
                kin_energy = np.zeros(rho.shape)
                mag_energy = np.zeros(rho.shape)
            data_out[index, 1] = rho.sum()
            data_out[index, 2] = kin_energy.sum()
            data_out[index, 3] = rot_energy.sum()
            data_out[index, 4] = mag_energy.sum()
            data_out[index, 5] = conv_energy.sum()
            self.close_h5(data_h5)
        data_out[:, 1] = u.convert_to_solar_masses(data_out[:, 1])
        return data_out

    def __unbound_mass_energy(self):
        dV = self.cell.dVolume_integration(self.ghost)
        file_list = self.file_list_hdf()
        data_out = np.zeros((len(file_list), 2))
        r = self.cell.radius(self.ghost)
        while r.ndim != dV.ndim:
            r = r[None, ...]
        for (index, file) in zip(range(len(file_list)), file_list):
            data_h5 = self.open_h5(file)
            rho = self.rho(data_h5)
            en = dV * (self.energy_MHD(data_h5) + rho * self.grav_pot(data_h5))
            rho *= dV
            summ_index = np.where((r <= 1e10) & (en > 0), 1, 0)
            data_out[index, 0] = (en * summ_index).sum()
            data_out[index, 1] = (rho * summ_index).sum()
        data_out[:, 1] = u.convert_to_solar_masses(data_out[:, 1])
        return data_out

    def stream_function(self, data_h5, plane):
        return strfct2D(self.magnetic_field_ct(data_h5),self.cell, self.ghost, plane)

    def get_growth(self):
        path = os.path.join(self.storage_path, 'magnetic_growth.h5')
        if not os.path.exists(path):
            g = self.__growth()
            file_out = h5py.File(path, 'w')
            file_out.create_dataset('time', data = g[:, 0])
            file_out.create_dataset('growth', data = g[:, 1])
            self.close_h5(file_out)
            return g
        g = h5py.File(path)
        d_out = np.stack((np.array(g['time']), np.array(g['growth'])), axis = -1)
        return d_out

    def __B_growth(self, data_h5):
        return self.magnetic_field(data_h5)[..., 0] * \
               self.omega_radial_derivative(data_h5) * \
               self.cell.radius(self.ghost)

    def __growth(self):
        if self.dim != 2:
            raise FileNotFoundError("Not implemented yet :\"(")
        file_list = self.file_list_hdf()
        data_out = np.zeros((len(file_list), 2))
        indices, PNS_radius, data_out[:, 0] = self.get_PNS_radius(indices = True,
                                              min_max_average = True, ret_time = True,
                                              tob_corrected = True)
        radius = self.cell.radius(self.ghost)
        dr = self.cell.dr(self.ghost)
        rad_fs = np.argmax(radius >= 5e5)
        DeltaR = radius[rad_fs]
        theta = self.cell.theta(self.ghost)
        theta_fs, theta_ls = np.argmax(theta >= (np.pi / 3)), np.argmax(theta > (2 / 3 * np.pi))
        for (i, f) in zip(indices, file_list):
            data_h5 = self.open_h5(f)
            rad_ls = np.argmax(radius >= PNS_radius[i, 3])
            g = self.__B_growth(data_h5)[theta_fs : theta_ls, rad_fs : rad_ls] * \
                        dr[rad_fs : rad_ls]
            data_out[i, 1] = np.mean(np.sum(g, axis = -1)) / (radius[rad_ls] - DeltaR)
            self.close_h5(data_h5)
        return data_out
#Modes


#rho decomposition
    def pns_rho_decomposition(self):
        PNS_radius, indices, \
        time, ghost_cells =  self.get_PNS_radius(PNS_radius = True,indices = True, 
                                                 ret_time = True, ghost_cells = True,)
        data_out = np.zeros((time.shape[0], 3))
        data_out[:, 0] = time
        dtheta = self.cell.dtheta_integration(self.ghost)
        legp = LegP(np.cos(self.cell.theta(self.ghost)))
        polynomial_0 = getattr(legp, "P_" + str(0))
        polynomial_1 = getattr(legp, "P_" + str(1))
        polynomial_2 = getattr(legp, "P_" + str(2))
        dr = self.cell.dr_integration(self.ghost)
        r = np.ones((dtheta.shape[0], dr.shape[0])) * self.cell.radius(self.ghost)
        index = np.argmax(r >= u.convert_to_cm(30))
        file_list = self.file_list_hdf()
        for (i, file) in zip(indices, file_list):
            data_h5 = self.open_h5(file)
            rho = self.rho(data_h5)
            data_out[i, 1] = (rho[..., index] * dtheta * polynomial_1()).sum() / (rho[..., index] * dtheta * polynomial_0()).sum()
            data_out[i, 2] = (rho[..., index] * dtheta * polynomial_2()).sum() / (rho[..., index] * dtheta * polynomial_0()).sum()
        return data_out

    def frequency(self):
        data = self.pns_rho_decomposition()
        return np.stack((data[:, 0], stft(data[:, 1], nperseg= 15), stft(data[:, 2], nperseg= 15)), axis=-1)

    def get_rho_modes(self, order, tob_corrected = True):
        assert type(order) == int, "Order of the polynomial must be an integer"
        assert order < 7, "Order of the polynomials mus be at most 6"
        path = os.path.join(self.storage_path, "rho_modes.h5")
        if not os.path.exists(path):
            warnings.warn("Creating file...\nPlease wait.")
            self.__save_modes_rho(path, tob_corrected)
    
    def __save_modes_rho(self, path, tob_corrected):
        file_out = h5py.File(path, 'w')
        file_out.create_dataset("tob_correction", data = tob_corrected)
        modes_rho_group = file_out.create_group("rho_modes")
        modes_rho_v_group = file_out.create_group("rho_vel_modes")
        for o in range(7):
            mode = self.__rho_modes(o, tob_corrected)
            modes_rho_group.create_dataset("in_PNS_order_" + str(o), data = mode[:, 1])
            modes_rho_group.create_dataset("out_PNS_order_" + str(o), data = mode[:, 2])
            modes_rho_v_group.create_dataset("in_PNS_order_" + str(o), data = mode[:, 3])
            modes_rho_v_group.create_dataset("out_PNS_order_" + str(o), data = mode[:, 4])
        file_out.create_dataset("time", data = mode[:, 0])
        self.close_h5(file_out)
    
    def __rho_modes(self, order, tob_corrected):
        shock_radius, indices, \
        time, ghost_cells_shock = self.get_shock_radius(shock_radius = True,
                            indices = True, ret_time = True,ghost_cells = True,
                            tob_corrected = tob_corrected)
        PNS_radius, ghost_cells_pns =  self.get_PNS_radius(PNS_radius = True,
                                                           ghost_cells = True)
        data_out = np.zeros((time.shape[0], 5))
        data_out[:, 0] = time
        dtheta = self.cell.dtheta_integration(self.ghost)
        legp = LegP(np.cos(self.cell.theta(self.ghost)))
        polynomial = getattr(legp, "P_" + str(order))
        dr = self.cell.dr_integration(self.ghost)
        r = np.ones((dtheta.shape[0], dr.shape[0])) * self.cell.radius(self.ghost)
        file_list = self.file_list_hdf()
        for (i, file) in zip(indices, file_list):
            data_h5 = self.open_h5(file)
            rho = self.rho(data_h5) * dr
            v_rad = self.radial_velocity(data_h5)
            self.close_h5(data_h5)
            mask_in_PNS = r <= self.ghost.remove_ghost_cells_radii(PNS_radius[..., i], 
                                                                    self.dim, **ghost_cells_pns)[:, None]
            mask_out_PNS = (r > self.ghost.remove_ghost_cells_radii(PNS_radius[..., i], 
                                                                    self.dim, **ghost_cells_pns)[:, None] )& \
                            (r <= self.ghost.remove_ghost_cells_radii(shock_radius[..., i], 
                                                                    self.dim, **ghost_cells_shock)[:, None])
            data_out[i, 1] = 0.5 * (rho[mask_in_PNS].sum(axis = -1) * dtheta * polynomial()).sum()
            data_out[i, 2] = 0.5 * (rho[mask_out_PNS].sum(axis = -1) * dtheta * polynomial()).sum()
            data_out[i, 3] = 0.5 * ((rho[mask_in_PNS] * v_rad[mask_in_PNS]).sum(axis = -1) \
                                    * dtheta * polynomial()).sum()
            data_out[i, 4] = 0.5 * ((rho[mask_out_PNS] * v_rad[mask_out_PNS]).sum(axis = -1) \
                                    * dtheta * polynomial()).sum()
        return data_out
           
# dipole
    def get_modes(self, order, normalized = False, tob_corrected = True):
        assert type(order) == int, "Order of the polynomial must be an integer"
        assert order < 7, "Order of the polynomials mus be at most 6"
        path = os.path.join(self.storage_path, "GW_modes.h5")
        if not os.path.exists(path):
            warnings.warn("Creating file...\nPlease wait.")
            self.__save_modes_GW(path, tob_corrected)
        data = h5py.File(path)
        time = np.array(data['dipole_modes']['time'])
        mode = np.array(data['dipole_modes']['order_' + str(order)])
        if normalized:
            mode /= np.array(data['dipole_modes']['order_0'])
        if not tob_corrected and data['tob_correction']:
            time += self.time_of_bounce_rho()
        elif tob_corrected and not data['tob_correction']:
            time -= self.time_of_bounce_rho()
        return np.stack((time, mode), axis = -1)

    def __save_modes_GW(self, path, tob_corrected):
        file_out = h5py.File(path, 'w')
        GW = self.GW_Amplitudes(tob_corrected)
        file_out.create_dataset("tob_correction", data = tob_corrected)
        GW_group = file_out.create_group("GWs")
        if self.dim == 1:
            GW_group.create_dataset("time", data = 0)
            GW_group.create_dataset("GWs", data = 0)
        else:
            GW_group.create_dataset("time", data = GW[:, 0])
            GW_group.create_dataset("GWs", data = GW[:, 1])
        modes_group = file_out.create_group("dipole_modes")
        for o in range(7):
            mode = self.__dipole(o, tob_corrected)
            modes_group.create_dataset("order_" + str(o), data = mode[:, 1])
        modes_group.create_dataset("time", data = mode[:, 0])
        self.close_h5(file_out)
        
    def __dipole(self, order, tob_corrected):
        """
        The method used to derive the dipole modes is implemented following
        A. Summa, F. Hanke, H-T Janka, T Melson, A. Marek, and B. Müller
        DOI 10.3847/0004-637X/825/1/6
        https://iopscience.iop.org/article/10.3847/0004-637X/825/1/6
        """
        shock_radius, indices, \
        time, ghost_cells = self.get_shock_radius(shock_radius = True,
                            indices = True, ret_time = True,ghost_cells = True,
                            tob_corrected = tob_corrected)
        data_out = np.zeros((time.shape[0], 2))
        data_out[:, 0] = time
        dOmega = self.cell.dOmega(self.ghost)
        legp = LegP(np.cos(self.cell.theta(self.ghost)))
        polynomial = getattr(legp, "P_" + str(order))()
        for i in indices:
            data_out[i, 1] = 0.5 * (self.ghost.remove_ghost_cells_radii(shock_radius[..., i], 
                                                                    self.dim, **ghost_cells) * \
                             polynomial * dOmega).sum()
        return data_out

    def __fluid_trajectories(self):
        file_list = self.file_list_hdf()
        mass_step = 0.02
        mass_max = 2.5
        mass_min = mass_step
        radius = self.cell.radius(self.ghost)
        data_out = np.zeros((len(file_list), int(mass_max / mass_step)))
        for (i, f) in zip(range(data_out.shape[0]), file_list):
            data_h5 = self.open_h5(f)
            mass = self.mass_shells(data_h5)
            tot_mass = mass_min
            self.close_h5(data_h5)
            for m_index in range(data_out.shape[1]):
                data_out[i, m_index] = radius[np.argmax(mass >= tot_mass)]
                tot_mass += mass_step
        return data_out
    
    def bounce_compactness(self, mass = None):
        """
        Compactness parameter as defined in O'Connor 2011  
        `10.1088/0004-637X/730/2/72`
        It returns the compacteness at bounce (maximum compactness of
        the core) for the model. 
        If a mass at which the compactness has to be provided is specified
        it will return just a value
        """
        bounce_file = self.find_file_from_time(0, True)
        data_h5 = self.open_h5(bounce_file)
        rho = self.rho(data_h5)
        radius = u.convert_to_km(self.cell.radius(self.ghost)) / 1000
        dV = self.cell.dVolume_integration(self.ghost)
        self.close_h5(data_h5)
        if self.dim > 1:
            enclosed_mass = np.sum( rho * dV, axis = tuple(range( self.dim - 1) ) )
        enclosed_mass = u.convert_to_solar_masses( np.cumsum(enclosed_mass) )
        compactness = enclosed_mass / radius
        if mass is not None:
            return compactness[ np.argmax( enclosed_mass >= mass ) ]
        else:
            return compactness
