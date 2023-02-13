import platform
import os

def pltf():
    platf = platform.system()
    """
    The only platforms "supported" are Linux and Windows.
    """
    if not platf in ['Windows', 'Linux']:
        raise TypeError("Platform not suppoerted, sorry :)")
    return platf

def local_storage_folder(platf, folder_name, dim):
    """
    Some calculated quantities need are computationally expensive, therefore
    they need to be stored in a convenient place (home folder in Linux or 
    Desktop in Windows).
    """
    if platf == 'Windows':
        storage_path = os.path.expanduser('~/Desktop/Simulations_analysis_quantities')
    else:
        storage_path = os.path.expanduser('~/Simulations_analysis_quantities')
    if dim == 1:
        dim_folder = '1D'
    elif dim == 2:
        dim_folder = '2D'
    else:
        dim_folder = '3D'
    storage_path = os.path.join(storage_path, dim_folder)
    if not os.path.exists(storage_path):
        os.mkdir(storage_path)
    storage_path = os.path.join(storage_path, folder_name)
    if not os.path.exists(storage_path):
        os.mkdir(storage_path)
    return storage_path