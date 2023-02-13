#Based on an Martin Obergaulinger's IDL script
import numpy as np

"""
Stramlines calculation by integrating tangents to the vector field.
In this particular case we are interested in calculating the streamlines
of the magnetic field.
"""
def strfunction2D(b1, b2, ax, ay, az, lx, ly, lz, plane):
    if plane == 'yz':
        dF1 = 0
        dF2 = 0
        dl1 = 1
        dl2 = 2
    elif plane == 'xz':
        dF1 =   b2 * az
        dF2 = - b1 * ax
        dl1 = ly
        dl2 = ly
    elif plane == 'xy':
        dF1 =   b2 * ay
        dF2 = - b1 * ax
        dl1 = lz
        dl2 = lz
    
    sb = b1.shape
    n1 = sb[0]+1 #theta
    n2 = sb[1]+1 #radius
    
    #stream fct
    F = np.zeros((n1, n2))
    dl1inv = 1./dl1
    dl1inv = np.nan_to_num(dl1inv, nan=0)
    dl2inv = 1./dl2
    dl2inv = np.nan_to_num(dl2inv, nan=0)

    #integrate the streamfunction
    for j in range(1, n1):
        jj = min(j, n1-2)
        F[j, 0] = (F[j, 0] * dl1[j-1, 0] + dF2[j-1, 0]) * dl2inv[jj, 0]
    for i in range(1, n2):
        ii = min(i, n2-2)
        F[0, i] = (F[0, i-1] * dl2[0, i-1] + dF1[0, i-1]) * dl2inv[0, ii]
        for j in range(1, n1):
            jj = min(j, n1-2)
            F[j, i] = (F[j-1, i] * dl1[j-1, ii] + dF2[j-1,ii]) * dl1inv[jj, ii]

    F0 = F[0:n1-1, 0:n2-1] + F[1:n1, 0:n2-1] + F [0:n1-1, 1:n2] + F[1:n1, 1:n2]

    if plane == 'yz':
        F0[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]
        F[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]*lx
    elif plane=='xz':
        F0[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]
        F[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]*ly
    elif plane == 'xy':
        F0[0:n1-1, 0:n2-1] = F[0:n1-1, 0:n2-1]
        F[0:n1-1, 0:n2-1 ] = F[0:n1-1, 0:n2-1]*lz
    

    FF = F[0:n1-1, 0:n2-1] + F[1:n1, 0:n2-1] + F [0:n1-1, 1:n2] + F[1:n1, 1:n2]
    FF =  FF / 4.
    F0 = F0 / 4.
    return FF

def strfct2D(b, cell, ghost, plane):
    """
    plane: 'yz', 'xz', 'xy'
    cell: object cell
    b magnetic field
    """
    assert plane in ['yz', 'xz', 'xy'], "plane must be one of: \'yz\', \'xz\', \'xy\'"
    ax = cell.ax(ghost)
    ay = cell.ay(ghost)
    az = cell.az(ghost)
    lx = cell.lx(ghost)
    ly = cell.ly(ghost)
    lz = cell.lz(ghost)

    if plane == 'yz':
        b1 = b[:,:,1] #theta b
        b2 = b[:,:,2] #phi b
    elif plane == 'xz':
        b1 = b[:,:,0] #rad b
        b2 = b[:,:,2] #phi b
    elif plane == 'xy':
        b1 = b[:,:,0] #rad b
        b2 = b[:,:,1] #theta b

    F =  strfunction2D(b1, b2, ax, ay, az, lx, ly, lz, plane)  

    return F
