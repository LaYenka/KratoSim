# LBM settings

from numpy import *

class LBM_Init:

    def get_inivel(self,d):
     v = fromfunction(self.inivel,(d,self.nx,self.ny))
     return v
     
    def inivel(self,d, x, y):
       return (1-d) * self.uLB * (1 + 1e-4*sin(y/self.ly*2*pi))

