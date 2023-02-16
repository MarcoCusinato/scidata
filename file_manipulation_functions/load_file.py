import os
import numpy as np
import pandas as pd

def load_file(path_folder, file_name):
    """
    Load data files, three attempts are made
    - Attempt 1: use of loadtxt method of NumPy
    - Attempt 2: use of loadtxt method of NumPy, reading only the first n columns,
                 with n specified on the 2nd row of the .txt header file
    - Attempt 3: use of Pandas to generate missing data 
    """
    path = os.path.join(path_folder, file_name)
    assert os.path.exists(path), "Selected file does not exists"
    try:
        data = np.loadtxt(path)
    except:
        try:
            col_number = int(np.loadtxt(path.replace('.dat', '.txt'), skiprows=1) + 2)
            data = np.loadtxt(path, usecols=tuple(range(col_number)))
        except:
            head = list(np.genfromtxt(path.replace('.dat', '.txt'), dtype=str,
                                    delimiter=',', skip_footer=1))
            data_str = pd.read_table(path, dtype=str, sep='\s+', names=head)
            data_str = data_str.fillna('0')
            data_str = data_str.to_numpy()
            index_list = []
            for i in range(data_str.shape[0]):
                try:
                    data_str[i,:].astype('float')
                except:
                    index_list.append(i)
            data_str = np.delete(data_str, index_list,0)
            data = data_str.astype('float')
    return data
