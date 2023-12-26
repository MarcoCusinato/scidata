import numpy as np
from typing import Literal

def IDL_derivative(x, y, xvariable: Literal['radius', 'theta', 'phi'] = 'radius'):
    """
    Derivatie performed using three point Lagrangian interpolation, as in:
    `https://www.l3harrisgeospatial.com/docs/deriv.html` 
    """
    if not xvariable == 'radius':
        raise TypeError("Not implemented yet, sorry :)")
    assert x.shape == y.shape or x.shape[0] == y[..., :].shape[-1], \
                      "Arrays must have equal last dimension"
    while x.ndim != y.ndim:
        x = x[None, ...]
    assert x.shape[-1] >= 3, "To calculate this derivative you need AT LEAST three points."
    #first point
    x01 = x[..., 0] - x[..., 1]
    x02 = x[..., 0] - x[..., 2]
    x12 = x[..., 1] - x[..., 2]
    derivative = y[..., 0] * (x01 + x02) / \
        (x01 * x02) - y[..., 1] * x02 / (x01 * x12) \
        + y[..., 2] * x01 / (x02 * x12)
    #mid points
    x01 = x[..., : -2] - x[..., 1 : -1]
    x02 = x[..., : -2] - x[..., 2 :]
    x12 = x[..., 1 : -1] - x[..., 2 :]
    derivative = np.concatenate((derivative[..., None], 
                                 y[..., : -2] * x12 / (x01 * x02) + \
                                 y[..., 1 : -1] * (1. / x12 - 1. / x01) - \
                                 y[..., 2 :] * x01 / (x02 * x12)), axis = -1)
    #last point
    x01 = x[..., -3] - x[..., -2]
    x02 = x[..., -3] - x[..., -1]
    x12 = x[..., -2] - x[..., -1]
    
    return np.concatenate((derivative, (-y[...,-3] * x12 / (x01 * x02) +\
                                        y[..., -2] * x02 / (x01 * x12) - \
                                        y[..., -1] * (x02 + x12) / (x02 * x12))[..., None]),
                                        axis = -1)