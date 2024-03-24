# compute macroscopic values

from numpy import *


class LBM_macro:

  ###### Function Definitions ####################################################
  def macroscopic(self,fin):
    rho = sum(fin, axis=0)
    u = zeros((2, self.nx, self.ny))
    for i in range(9):
        u[0,:,:] += self.v[i,0] * fin[i,:,:]
        u[1,:,:] += self.v[i,1] * fin[i,:,:]
    u /= rho
    return rho, u
