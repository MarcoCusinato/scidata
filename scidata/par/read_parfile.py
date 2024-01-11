import os
import f90nml

def get_indices_from_parfile(file_name, path_folder):
    """
    Reads the simulation parfile and return a dictionary with the indices
    of the hydro and thd quantities.
    keywords: 
        first dictionary:
            'hydro' : hydrodynamics quantities
            'thd' : thermodynamics quantities
        second dictionary:

                                      HYDRO
            |----------------------------------------------------------|
            |'I_RH' : density           |   'I_EN' : entropy           |
            |'I_VX' : radial velocity   |   'I_YE' : electron fraction |
            |'I_VY' : polar velocity    |   'I_YZ' :                   |
            |'I_VZ' : azimutal velocity |   'I_LRTZ' : lorentz factor  |
            ------------------------------------------------------------
                                       THD
    """
    namelist = f90nml.read(os.path.join(path_folder, file_name))
    indices = {'hydro':{}, 'thd':{}}
    if 'I_RH'in namelist["IINDICES"]:
        indices['hydro']['I_RH'] = namelist["IINDICES"]["I_RH"] - 1
    if 'I_EN'in namelist["IINDICES"]:
        indices['hydro']['I_EN'] = namelist["IINDICES"]["I_EN"] - 1
    if 'I_VX'in namelist["IINDICES"]:
        indices['hydro']['I_VX'] = namelist["IINDICES"]["I_VX"] - 1
    if 'I_VY'in namelist["IINDICES"]:
        indices['hydro']['I_VY'] = namelist["IINDICES"]["I_VY"] - 1
    if 'I_VZ'in namelist["IINDICES"]:
        indices['hydro']['I_VZ'] = namelist["IINDICES"]["I_VZ"] - 1
    if 'I_YE'in namelist["IINDICES"]:
        indices['hydro']['I_YE'] = namelist["IINDICES"]["I_YE"] - 1
    if 'I_YN'in namelist["IINDICES"]:
        indices['hydro']['I_YN'] = namelist["IINDICES"]["I_YN"] - 1
    if 'I_EOSERR'in namelist["IINDICES"]:
        indices['thd']['I_EOSERR'] = namelist["IINDICES"]["I_EOSERR"] - 1
    if 'I_LRTZ'in namelist["IINDICES"]:
        indices['thd']['I_LRTZ'] = namelist["IINDICES"]["I_LRTZ"] - 1
    if 'I_DENS'in namelist["IINDICES"]:
        indices['thd']['I_DENS'] = namelist["IINDICES"]["I_DENS"] - 1    
    if 'I_EINT'in namelist["IINDICES"]:
        indices['thd']['I_EINT'] = namelist["IINDICES"]["I_EINT"] - 1
    if 'I_ENTH'in namelist["IINDICES"]:
        indices['thd']['I_ENTH'] = namelist["IINDICES"]["I_ENTH"] - 1
    if 'I_PELE'in namelist["IINDICES"]:
        indices['thd']['I_PELE'] = namelist["IINDICES"]["I_PELE"] - 1    
    if 'I_TELE'in namelist["IINDICES"]:
        indices['thd']['I_TELE'] = namelist["IINDICES"]["I_TELE"] - 1
    if 'I_NELE'in namelist["IINDICES"]:
        indices['thd']['I_NELE'] = namelist["IINDICES"]["I_NELE"] - 1
    if 'I_PION'in namelist["IINDICES"]:
        indices['thd']['I_PION'] = namelist["IINDICES"]["I_PION"] - 1
    if 'I_TION'in namelist["IINDICES"]:
        indices['thd']['I_TION'] = namelist["IINDICES"]["I_TION"] - 1    
    if 'I_NION'in namelist["IINDICES"]:
        indices['thd']['I_NION'] = namelist["IINDICES"]["I_NION"] - 1
    if 'I_VELX'in namelist["IINDICES"]:
        indices['thd']['I_VELX'] = namelist["IINDICES"]["I_VELX"] - 1
    if 'I_VELY'in namelist["IINDICES"]:
        indices['thd']['I_VELY'] = namelist["IINDICES"]["I_VELY"] - 1    
    if 'I_VELZ'in namelist["IINDICES"]:
        indices['thd']['I_VELZ'] = namelist["IINDICES"]["I_VELZ"] - 1
    if 'I_TMPR'in namelist["IINDICES"]:
        indices['thd']['I_TMPR'] = namelist["IINDICES"]["I_TMPR"] - 1
    if 'I_ENTR'in namelist["IINDICES"]:
        indices['thd']['I_ENTR'] = namelist["IINDICES"]["I_ENTR"] - 1
    if 'I_GAMM'in namelist["IINDICES"]:
        indices['thd']['I_GAMM'] = namelist["IINDICES"]["I_GAMM"] - 1    
    if 'I_HEAT'in namelist["IINDICES"]:
        indices['thd']['I_HEAT'] = namelist["IINDICES"]["I_HEAT"] - 1
    if 'I_DELP'in namelist["IINDICES"]:
        indices['thd']['I_DELP'] = namelist["IINDICES"]["I_DELP"] - 1
    if 'I_SMOMX'in namelist["IINDICES"]:
        indices['thd']['I_SMOMX'] = namelist["IINDICES"]["I_SMOMX"] - 1    
    if 'I_SMOMY'in namelist["IINDICES"]:
        indices['thd']['I_SMOMY'] = namelist["IINDICES"]["I_SMOMY"] - 1
    if 'I_SMOMZ'in namelist["IINDICES"]:
        indices['thd']['I_SMOMZ'] = namelist["IINDICES"]["I_SMOMZ"] - 1
    if 'I_PGAS'in namelist["IINDICES"]:
        indices['thd']['I_PGAS'] = namelist["IINDICES"]["I_PGAS"] - 1
    if 'I_CSND'in namelist["IINDICES"]:
        indices['thd']['I_CSND'] = namelist["IINDICES"]["I_CSND"] - 1    
    if 'I_COMP'in namelist["IINDICES"]:
        list_comp_ind = namelist["IINDICES"]["I_COMP"]
        list_comp_ind = [x -1 for x in list_comp_ind]
        indices['thd']['I_COMP'] = list_comp_ind
    if 'I_CPOT'in namelist["IINDICES"]:
        list_cpot_ind = namelist["IINDICES"]["I_CPOT"]
        list_cpot_ind = [x -1 for x in list_cpot_ind]
        indices['thd']['I_CPOT'] = list_cpot_ind
    if 'I_BHEX'in namelist["IINDICES"]:
        indices['thd']['I_BHEX'] = namelist["IINDICES"]["I_BHEX"] - 1    
    if 'STENCIL'in namelist["GRIDPARS"]:
        STENCIL =  namelist["GRIDPARS"]["STENCIL"]
    for  k1 in indices['hydro']:
        if type(indices['hydro'][k1]) == list:
            for i in range(0,len(indices['hydro'][k1])):
                if indices['hydro'][k1][i] < 0:
                    indices['hydro'][k1][i] = None
        elif indices['hydro'][k1] < 0:
            indices['hydro'][k1] = None
    for  k1 in indices['thd']:
        if type(indices['thd'][k1]) == list:
            for i in range(0, len(indices['thd'][k1])):
                if indices['thd'][k1][i] < 0:
                    indices['thd'][k1][i] = None
        elif indices['thd'][k1] < 0:
            indices['thd'][k1] = None
    return indices, STENCIL

def get_initial_parameters(path_folder):
    initial_parameters = {"omgadd": 0.0,
                          "omgmult": 0.0,
                          "bfact": 0.0,
                          "btorfact": 0.0,
                          "b0": 0.0,
                          "bt": 0.0
                          }
    start_pars = f90nml.read(os.path.join(path_folder, 'start.pars'))
    if start_pars["PHYSSYST"]["RELATIVISTIC"]:
        initial_parameters['gravity'] = 'Pseudo-relativistic'
    else:
        initial_parameters['gravity'] = 'Newtonian'
    initial_parameters["gravitational_potential"] = start_pars["GRAVPARS"]["MDPOT"]
    initial_parameters["lapse_function"] = start_pars["GRAVPARS"]["LAPSE_FORM"]
    initial_parameters["NS_EOS"] = (os.path.basename(start_pars["SHENEOSPARS"]["SHEN_TBFILE"]))
    initial_parameters["NS_EOS"] = initial_parameters["NS_EOS"].split('.')[:-1]
    if type(initial_parameters["NS_EOS"]) == list:
        if len(initial_parameters["NS_EOS"]) == 1:
            initial_parameters["NS_EOS"] = initial_parameters["NS_EOS"][0]
        else:
            tmp = initial_parameters["NS_EOS"][0]
            for i in initial_parameters["NS_EOS"][1:]:
                tmp += '.' + i
            initial_parameters["NS_EOS"] = tmp
    
    files = os.listdir(path_folder)
    files = [x for x in files if not x == 'start.pars' and not x == '.run']
    files.sort()
    Found_heger = False
    for file in files:
        namelist = f90nml.read(os.path.join(path_folder, file))
        if 'Hegerpars' in namelist:
            Found_heger = True
            initial_parameters["Heger_model"] = namelist["Hegerpars"]["Heger_model"]
            if "omgadd" in namelist["Hegerpars"]:
                initial_parameters["omgadd"] = namelist["Hegerpars"]["omgadd"]
            if "omgmult" in namelist["Hegerpars"]:
                initial_parameters["omgmult"] = namelist["Hegerpars"]["omgmult"]
            if "bfield" in namelist["Hegerpars"]:
                if namelist["Hegerpars"]["bfield"]:
                    if "bfact" in namelist["Hegerpars"]:
                        initial_parameters["bfact"] = namelist["Hegerpars"]["bfact"]
                    if "btorfact" in namelist["Hegerpars"]:
                        initial_parameters["btorfact"] = namelist["Hegerpars"]["btorfact"]
        if "axivecpotpars" in namelist:
            if "b0" in namelist["axivecpotpars"]:
                initial_parameters["b0"] = namelist["axivecpotpars"]["b0"]
            if "bt" in namelist["axivecpotpars"]:
                initial_parameters["bt"] = namelist["axivecpotpars"]["bt"]
        if Found_heger:
            break
    if not Found_heger:
        run_path = os.path.join(path_folder, '.run')
        run_files = os.listdir(run_path).sort()
        for file in run_files:
            namelist = f90nml.read(os.path.join(run_path, file))
            if 'Hegerpars' in namelist and namelist["Hegerpars"]["Heger_model"]:
                Found_heger = True
                initial_parameters["Heger_model"] = namelist["Hegerpars"]["Heger_model"]
                if "omgadd" in namelist["Hegerpars"]:
                    initial_parameters["omgadd"] = namelist["Hegerpars"]["omgadd"]
                if "omgmult" in namelist["Hegerpars"]:
                    initial_parameters["omgmult"] = namelist["Hegerpars"]["omgmult"]
                if "bfield" in namelist["Hegerpars"]:
                    if namelist["Hegerpars"]["bfield"]:
                        if "bfact" in namelist["Hegerpars"]:
                            initial_parameters["bfact"] = namelist["Hegerpars"]["bfact"]
                        if "btorfact" in namelist["Hegerpars"]:
                            initial_parameters["btorfact"] = namelist["Hegerpars"]["btorfact"]
            if "axivecpotpars" in namelist:
                if "b0" in namelist["axivecpotpars"]:
                    initial_parameters["b0"] = namelist["axivecpotpars"]["b0"]
                if "bt" in namelist["axivecpotpars"]:
                    initial_parameters["bt"] = namelist["axivecpotpars"]["bt"]
            if Found_heger:
                break
    if not Found_heger:
        initial_parameters["omgadd"] = None
        initial_parameters["omgmult"] = None
        initial_parameters["Heger_model"] = None
        initial_parameters["bfact"] = None
        initial_parameters["btorfact"] = None
        initial_parameters["b0"] = None
        initial_parameters["bt"] = None
    return initial_parameters
