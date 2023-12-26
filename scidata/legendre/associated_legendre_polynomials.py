import numpy as np
from scipy.special import factorial, binom
from scidata.legendre.legendre_polynomials import legendre_polynomials as legP

class associated_legendre_polynomials:
    """
    Class that allows the calculation of the legendre associated polynomials
    Parameters:
        x: (npumpy array) quantity to use as polinomial generator
    Methods:
        polynomial: 
    Legendre associated polynomials are calculated following:
    `https://en.wikipedia.org/wiki/Associated_Legendre_polynomials`
    It calculates the associated polynomials 
    ```
        P^m_l
    ```
    to change the indices positioning use the "change_indices" method.
    """
    def __init__(self, x, l, m = 0):
        assert l >= m, "Must insert l >= m."
        self.variable = x
        self.l = l
        self.m = m
    
    def polynomial(self):
        if self.m == 0 and self.l <= 6:
            pol = legP(self.variable)
            return getattr(pol, "P_" + str(self.l))()
        elif self.l == self.m and self.l <= 4:
            return getattr(self, "_P" + str(self.m) + "_" + str(self.l))()
        else:
            return self._norm() * self._sum()

    def change_indices(self):
        return (-1) ** self.m * self.polynomial()
    
    def _sum(self):
        k = np.arange(self.m, self.l + 1)
        sum = factorial(k) / factorial(k - self.m)
        sum *= binom(self.l, k) * binom((self.l + k - 1) / 2, self.l) 
        sum = self.variable[..., None] **(k - self.m) * sum
        return sum.sum(axis = -1) 
    

    def _norm(self):
        if self.m == 0:
            return 2 ** self.l
        return (-1) ** self.m * 2 ** self.l * (1 - self.variable ** 2) ** (self.m / 2)
    
    def _P1_1(self):
        return -(1 - self.variable ** 2) ** 0.5
    
    def _P2_2(self):
        return 3. * (1 - self.variable **2)
    
    def _P3_3(self):
        return -15 * (1 - self.variable ** 2) ** 1.5
    
    def _P4_4(self):
        return 105. * (1 - self.variable ** 2) ** 2
    