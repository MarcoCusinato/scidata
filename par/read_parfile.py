import os
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
    f = open(os.path.join(path_folder, file_name), 'r')
    Lines = [line for line in f.readlines() if line.strip()]
    f.close()
    indices = {'hydro':{}, 'thd':{}}
    for li in Lines:
        li.replace(",", "")
        words = li.split()
        if words[0] == 'I_RH':
            indices['hydro']['I_RH'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_EN':
            indices['hydro']['I_EN'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_VX':
            indices['hydro']['I_VX'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_VY':
            indices['hydro']['I_VY'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_VZ':
            indices['hydro']['I_VZ'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_YE':
            indices['hydro']['I_YE'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_YZ':
            indices['hydro']['I_YZ'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_LRTZ':
            indices['thd']['I_LRTZ'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_DENS':
            indices['thd']['I_DENS'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_EINT':
            indices['thd']['I_EINT'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_ENTH':
            indices['thd']['I_ENTH'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_PELE':
            indices['thd']['I_PELE'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_TELE':
            indices['thd']['I_TELE'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_NELE':
            indices['thd']['I_NELE'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_PION':
            indices['thd']['I_PION'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_TION':
            indices['thd']['I_TION'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_NION':
            indices['thd']['I_NION'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_VELX':
            indices['thd']['I_VELX'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_VELY':
            indices['thd']['I_VELY'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_VELZ':
            indices['thd']['I_VELZ'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_TMPR':
            indices['thd']['I_TMPR'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_ENTR':
            indices['thd']['I_ENTR'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_GAMM':
            indices['thd']['I_GAMM'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_HEAT':
            indices['thd']['I_HEAT'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_DELP':
            indices['thd']['I_DELP'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_SMOMX':
            indices['thd']['I_SMOMX'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_SMOMY':
            indices['thd']['I_SMOMY'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_SMOMZ':
            indices['thd']['I_SMOMZ'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_PGAS':
            indices['thd']['I_PGAS'] = int(words[2].replace(",", "")) - 1
        if words[0] == 'I_CSND':
            indices['thd']['I_CSND'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'I_COMP':
            list_comp_ind = []
            for i in range(2,len(words)):
                list_comp_ind.append(int(words[i].replace(",", "")) - 1)
            indices['thd']['I_COMP'] = list_comp_ind
        if words[0] == 'I_CPOT':
            list_cpot_ind = []
            for i in range(2,len(words)):
                list_cpot_ind.append(int(words[i].replace(",", "")) - 1)
            indices['thd']['I_CPOT'] = list_cpot_ind
        if words[0] == 'I_BHEX':
            indices['thd']['I_BHEX'] = int(words[2].replace(",", "")) - 1    
        if words[0] == 'STENCIL':
           STENCIL =  int(words[2].replace(",", ""))
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

