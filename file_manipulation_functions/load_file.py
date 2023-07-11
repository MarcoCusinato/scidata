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
    - Attempt 4: read file line by line filling a matrix with dimension NxM, where N
                 is the number of rows in the file and M the maximum number of colums.
                 It assumes the last line as the one with the more colums.
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
            try:
                head = list(np.genfromtxt(path.replace('.dat', '.txt'), dtype=str,
                                        delimiter=',', skip_footer=1))
                data_str = pd.read_table(path, dtype=str, sep='\s+', names=head,
                                        usecols=range(col_number))
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
            except:
                with open(path, 'r') as f:
                    lines = f.readlines()
                    data = np.zeros( ( len( lines ), len( lines[-1].split() ) ) )
                    for (i, line) in zip( range( len( lines ) ), lines ):
                        line_data = np.array( line.split() ).astype(float)
                        data[i, :len( line_data )] = line_data
    return data

def find_column_changing_line(path_folder, file_name):
    """
    Loads a data file and returns the number of the line at which
    there is a numbero of columns change.
    """
    path = os.path.join(path_folder, file_name)
    assert os.path.exists(path), "Selected file does not exists"
    number_of_colums = None
    line_number = 1
    with open(path, 'r') as f:
        for line in f:
            columns = len(line.strip().split())
            if number_of_colums is None:
                number_of_colums = columns
            if number_of_colums != columns:
                line_number += 1
                break
            line_number += 1
    if line_number < 3:
        line_number = None
    return line_number
