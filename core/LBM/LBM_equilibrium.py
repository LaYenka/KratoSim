# compute equilibrium of the ditribution function

from numpy import *

class LBM_equilibrium:

   def get_equilibrium(self,rho, vel):
    feq = self.equilibrium(rho, vel)
    return feq

   def equilibrium(self,rho, u):              # Equilibrium distribution function.
    usqr = 3/2 * (u[0]**2 + u[1]**2)
    feq = zeros((9,self.nx,self.ny))
    for i in range(9):
        cu = 3 * (self.v[i,0]*u[0,:,:] + self.v[i,1]*u[1,:,:])
        feq[i,:,:] = rho*self.t[i] * (1 + cu + 0.5*cu**2 - usqr)
    return feq
