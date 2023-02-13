from scidata.legendre.associated_legendre_polynomials import associated_legendre_polynomials as legP
import numpy as np
from scipy.special import factorial

class spherical_harmonics:
    """
    Class to calculate spherical hermonics decomposition of a 2D or 3D
    quantity.
    Conventions follow Bugli et al.
    `https://arxiv.org/pdf/2210.05012.pdf`
    in parameters:
        l: polynomial degree
        m: polynomial order
        theta: polar angle
        phi: azimutal angle
    """
    def __init__(self, l, m = 0, theta = None, phi = None):
        if not phi and m!=0:
            raise ValueError("m MUST be 0 in a 2D decomposition.")
        assert theta is not None, "polar angle must be an array"
        self.theta = theta
        self.phi = phi
        self.l = l
        self.m = m

    def Nm_l(self):
        """
        returns N^m_l coefficient
        """
        first_term = (2 * self.l +1) / (4. * np.pi)
        if self.m == 0:
            return np.sqrt(first_term)
        else:
            second_term = factorial(self.l - self.m) / \
                          factorial(self.l + self.m)
            return np.sqrt(first_term * second_term)

    def Ym_l(self):
        """
        returns the spherical harmonics Y^m_l
        """
        N = self.Nm_l()
        P = legP(np.cos(self.theta), self.l, self.m).polynomial()
        if self.m == 0:
            if self.phi:
                return N * P[..., None] * np.ones(self.phi.shape)
            else:
                return N * P
        else:
            return N * P[..., None] * np.exp(1j * self.m * self.phi)
    
    def Y_ml(self):
        """
        returns Y_{ml} spherical harmonics
        """
        Ym_l = self.Ym_l()
        if self.m == 0:
            return Ym_l
        elif self.m < 0:
            return np.sqrt(2) * (-1) ** self.m * Ym_l.imag
        else:
            return np.sqrt(2) * (-1) ** self.m * Ym_l.real