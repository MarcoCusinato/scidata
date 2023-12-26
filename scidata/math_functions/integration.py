import numpy as np

def backward_integration_index(data_to_integrate, threshold):
    """
    Calclulates the index at which the integrated quantity goes
    above a certain threshold.
    """
    np.nan_to_num(data_to_integrate, False,0)
    data_to_integrate = np.cumsum(np.flip(data_to_integrate, axis = -1), 
                                  axis = -1)
    return np.argmax(data_to_integrate >= threshold, axis = -1)
