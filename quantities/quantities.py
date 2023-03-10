import os
import sys
import warnings
from typing import Literal

import h5py
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import stft

dirname = os.path.dirname(__file__)
sys.path.append(os.path.join(dirname,'../..'))
from scidata.units.units import units

from scidata.cell.cell import cell as cl
from scidata.cell.ghost import ghost as gh
from scidata.file_manipulation_functions.load_file import load_file
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

    def entropy(self, data_h5):
        data = np.array(data_h5['thd']['data'])[:,:,:,self.hydroTHD_index['thd']['I_ENTR']]
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
    
    def eta_electrons(self, data_h5):
        chem_pot = self.chemical_potential_electrons(data_h5)
        temp = self.temperature(data_h5)
        return chem_pot/temp
    
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

    #neutrinos
    def neutrino_energy_grid(self):
        bin_energies = np.loadtxt(os.path.join(self.grid_path, self.nu_e_grid))
        return bin_energies[:,2]
    
    def neutrino_dE(self):
        bin_energies = np.loadtxt(os.path.join(self.grid_path, self.nu_e_grid))
        return bin_energies[:,3] - bin_energies[:,1]

    def neutrino_flux(self, data_h5):
        """
        Last index is the neutrino flavour: nue, nua, nux.
        Second-last index is the energy bin index.
        """
        nu_out = np.array(data_h5['neutrino']['e'])[:,:,:,:,:,1]
        nu_out[:, :, :, :, 2] /= 4
        return self.ghost.remove_ghost_cells(np.squeeze(nu_out), self.dim)

    def neutrino_flux_bin_integrated(self, data_h5):
        """
        Last index is the neutrino flavour: nue, nua, nux
        """
        nu = self.neutrino_flux(data_h5)
        return nu.sum(axis = self.dim)

    def neutrino_number_flux(self, data_h5, ret_nu_flux = False):
        """
        Last index is the neutrino flavour: nue, nua, nux
        """
        nu = self.neutrino_flux(data_h5)
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
        data = np.array(data_h5['neutrino']['oe'])[:,:,:,:,:,1]
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
        neutrino_fl, neutrino_num_fl = self.neutrino_number_flux(data_h5, True)
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
    
    #GWs
    def GW_Amplitudes(self, correct_for_tob=True):
        """
        GWs amplitudes calculated as the first partial time derivative
        of NE_220 and not with the second partial time derivative of ME_220.
        Moreover we only consider matter contribution to the amplitude, being
        the main contribution to it.
        """
        data = load_file(self.log_path, self.grw)
        const = -np.sqrt(15/np.pi)
        if correct_for_tob:
            data[:,2] -= self.time_of_bounce_rho()
        return np.stack((data[:, 2],
                         const * IDL_derivative(data[:,2], data[:,5])),
                         axis = 1)
    
    def GW_dimensionless_strain(self, distance, angle, correct_for_tob=True):
        data = self.GW_Amplitudes(correct_for_tob)
        const = np.sqrt(np.pi/15)*np.sin(angle)**2/distance
        data[:,1] *= const
        return data
    
    #misc from h5 file: grav pot, time
    def grav_pot(self, data_h5):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(data_h5['gravpot']['data'])), self.dim)

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
        rho_index = np.argmax(rho_data[:,1]>2.5e14)
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
        neutrino_fluxes = self.neutrino_flux(data_h5)
        neutrino_opacities = self.neutrino_opacity_flux(data_h5)
        dr = self.cell.dr(self.ghost)[..., None]
        while (neutrino_fluxes.ndim - 1) != dr.ndim:
            dr = dr[None, ...]
        neutrino_data = np.sum(neutrino_fluxes*neutrino_opacities, axis = self.dim) / \
                        np.sum(neutrino_fluxes, axis = self.dim) * dr
        np.nan_to_num(neutrino_data, False, 0)
        neutrino_data = np.flip(neutrino_data, axis = -2)
        ind = np.argmax(np.cumsum(neutrino_data, axis = -2) >= tau, axis = -2)
        radius = np.flip(self.cell.radius(self.ghost))
        if type(ind) == int:
            ind = np.zeros()
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
        
        ind = np.argmax((velocity_derivative <= vel_der_tresh) & (entropy_derivative <= entr_der_tresh) & \
                        (velocity >= vel_tresh) & (entropy >= entr_tresh), axis = -1)
        r = self.cell.radius(self.ghost)
        ind2 = velocity.shape[-1] - 1 - np.argmax((np.flip(velocity, axis = -1) <= vel_tresh) & \
                                        (np.flip(entropy, axis = -1) <= entr_tresh), axis = -1)
        oob = ind2 ==  velocity.shape[-1] - 1
        if any(oob):
            ind2 = np.where(oob, ind, ind2)
        oob = ind == 0
        if any(oob) and self.dim > 1:
            ind = np.where(oob, ind2, ind)
        shock = np.flip(r)[ind]
        mean_ind = shock != r[-1]
        oob = shock >= (shock.mean(where = mean_ind) + 1.5 * np.std(shock, where = mean_ind))
        if any(oob) and self.dim > 1:
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
        if any(oob) and self.dim == 2:
            theta = self.cell.theta(self.ghost)
            ind_to_interp = np.invert(oob)
            f = interp1d(theta[ind_to_interp], shock[ind_to_interp], kind = 'linear', fill_value="extrapolate")
            shock2 = f(theta)
            shock = np.where((oob) & (shock2 > 0), shock2, shock)
        if any(shock < 0):
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

    def __save_radii(self, radius_to_calculate: Literal['PNS', 'gain', 'shock', 'neutrino'], tob_correction, 
                save_name, **kwargs):
        methods = {'PNS': self.PNS_radius_single,
                   'gain': self.gain_radius_single,
                   'neutrino': self.neutrino_sphere_radius_single, 
                   'shock': self.shock_radius_single}
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
            data_h5 = self.open_h5(f)
            data[i, 0] = self.time(data_h5)
            if radius_to_calculate == 'gain':
                radius_single = calculate_radius(data_h5, PNS_radius[..., i])
            elif radius_to_calculate == 'shock' and data[i, 0] <= tob:
                radius_single = calculate_radius(data_h5, True)
            else:
                radius_single = calculate_radius(data_h5)
            if i == 0:
                rad_out = radius_single[..., None]
            else:
                rad_out = np.concatenate((rad_out, radius_single[..., None]),
                                          axis = -1)
            if radius_to_calculate == 'neutrino':
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
        #self.ghost.restore_default()
    
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
                        PNS_rot_ene = False, PNS_pol_kin_ene = False, ejected_mass = False, 
                        explosion_energy = False, gain_mass = False, gain_heating = False,
                        mass_accretion_500km = False, tob_corrected = True, name = 'mass_energy'):
        path = os.path.join(self.storage_path, name + '.h5')
        if not os.path.exists(path):
            warnings.warn("Selected file does not exists. Creating one now.")
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
        data_out = np.zeros((len(file_list), 5))
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
            rot_energy = rho * self.phi_velocity(data_h5) ** 2
            kin_energy = rho * (self.radial_velocity(data_h5) ** 2 \
                         + self.theta_velocity(data_h5) ** 2)
            mag_energy = self.magnetic_energy_per_unit_volume(data_h5) * PNS * dV
            data_out[index, 1] = rho.sum()
            data_out[index, 2] = kin_energy.sum()
            data_out[index, 3] = rot_energy.sum()
            data_out[index, 4] = mag_energy.sum()
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
            warnings.warn("Creatiing file...\nPlease wait.")
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
        A. Summa, F. Hanke, H-T Janka, T Melson, A. Marek, and B. M??ller
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
        polynomial = getattr(legp, "P_" + str(order))
        for i in indices:
            data_out[i, 1] = 0.5 * (self.ghost.remove_ghost_cells_radii(shock_radius[..., i], 
                                                                    self.dim, **ghost_cells) * \
                             polynomial() * dOmega).sum()
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
