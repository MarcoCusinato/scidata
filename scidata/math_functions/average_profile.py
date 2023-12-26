import numpy as np

def average_profile(dim, quantity):
    """
    Returs the average quantity profile wrt the angle(s).
    """
    if dim == 1:
        return average_1D(quantity)
    if dim == 2:
        return average_2D(quantity)
    return average_3D(quantity)

def average_1D(quantity):
    return quantity

def average_2D(quantity):
    return np.average(quantity, axis=0)

def average_3D(quantity):
    return np.average(quantity, axis=(0,1))
