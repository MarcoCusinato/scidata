first = True

if first:
    from scidata.platform.platform import local_storage_folder, pltf
    from scidata.paths.path_managment import get_paths_dictionary, save_paths_dictionary
    import os, scidata
    print("First import of scidata detected.")
    ## Create folder to store simulations paths
    print("Creating folder and dictionary for storing simulations paths.")
    if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(scidata.__file__)),
                          'paths/paths')):
        os.mkdir(os.path.join(os.path.dirname(os.path.abspath(scidata.__file__)),
                    'paths/paths'))
    ## Create dictionary to store simulations paths
    if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(scidata.__file__)),
                          'paths/paths/simulations_dictionary.py')):
        save_paths_dictionary(get_paths_dictionary(creating=True))
    ## Create folder to store quantities that need heavy computation
    print("Creating folder for storing quantities that need heavy computation.")
    local_storage_folder(pltf(), None, None)
    ## Creating symlink to scidata
    print("Creating symlink to scidata.")
    if pltf() == 'Windows':
        if not os.path.exists(os.path.expanduser('~/Desktop/scidata')):
            os.symlink(os.path.dirname(os.path.abspath(scidata.__file__)),
                       os.path.expanduser('~/Desktop/scidata'))
        else:
            print("Folder already exists.")
    elif pltf() == 'Linux':
        if not os.path.exists(os.path.expanduser('~/scidata')):
            os.symlink(os.path.dirname(os.path.abspath(scidata.__file__)),
                       os.path.expanduser('~/scidata'))
        else:
            print("Folder already exists.")
    else:
        raise TypeError("Platform not supported.")    
    with open(scidata.__file__) as f:
        lines = f.readlines()
    for (li, il) in zip(lines, range(len(lines))):
        if 'first' in li:
            li_index = il
            break
    lines[li_index] = 'first = False\n'
    with open(scidata.__file__, 'w') as f:
        f.writelines(lines)
