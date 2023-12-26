import numpy as np

class legendre_polynomials:
    """
    Class that allows the calculation of the legendre polynomials
    Parameters:
        x: (npumpy array) quantity to use as polinomial generator
    Methods:
        P_(from 0 to 4): legendre polynomials from order 0 to 4
            result: (numpy array) polynomial
    Legendre polynomials are calculated following:
    `https://en.wikipedia.org/wiki/Legendre_polynomials`
    """
    def __init__(self, x):
        self.variable = x
    
    def P_0(self):
        return np.ones(self.variable.shape)
    
    def P_1(self):
        return self.variable
    
    def P_2(self):
        return 0.5 * (3 * self.variable**2 - 1)
    
    def P_3(self):
        return 0.5 * (5 * self.variable**3 - 3 * self.variable)

    def P_4(self):
        return 0.125 * (35 * self.variable**4 - 30 * self.variable**2 + 3)
    
    def P_5(self):
        return 0.125 * (63 * self.variable**5 - 70 * self.variable**3 + 
               15 * self.variable)

    def P_6(self):
        return 0.0625 * (231 * self.variable**6 - 315 * self.variable**4 + 
               105 * self.variable**2 - 5)