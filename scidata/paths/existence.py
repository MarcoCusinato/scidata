import os

def check_path(storage_path, save_name):
    path = os.path.join(storage_path, save_name + '.h5')
    if os.path.exists(path):
        raise FileExistsError("This file already exists. Please change name.")
    return path
