# Lattice-Boltzmann class

from numpy import *
from .LBM_Init import *
from .LBM_equilibrium import *
from .macroscopic import *

class LBM(LBM_Init,LBM_equilibrium,LBM_macro):

   def __init__(self):
             self.maxIter = 200                        # Total number of time iterations.
             self.Re = 10000                           # Reynolds number.
             self.nx, self.ny = 420, 180               # Numer of lattice nodes.
             self.ly = self.ny-1                       # Height of the domain in lattice units.
             self.cx, self.cy, self.r = self.nx//4, self.ny//2, self.ny//9 # Coordinates of the cylinder.
             self.uLB     = 0.04                                 # Velocity in lattice units.
             self.nulb    = self.uLB*self.r/self.Re;             # Viscoscity in lattice units.
             self.omega = 1 / (3*self.nulb+0.5);                 # Relaxation parameter.
             self.v,self.t,self.col1,self.col2,self.col3 = [],[],[],[],[] # lattice characterstics
     

   def obstacle_fun(self,x, y):
              return (x-self.cx)**2+(y-self.cy)**2<self.r**2    
                         
     
   def run(self):
   
             # Imported methods
             from .latticeconstants import D2Q9
             from .latticesettings import SimSettings
             from .savedata import plot_field
             from .collision import collision
             from .streaming import streaming
   
             # initialize lattice
             self.v,self.t,self.col1,self.col2,self.col3 = D2Q9(self)
             try:
              print('----------------------------------------------------')
              print('lattice configuration')
              print('velocity',self.v.T)
              print('Gauss Weighta',self.t)
              print('----------------------------------------------------')
             except:
               print("An exception occurred - lattice velocity or Gauss weights not specified")

             
             # initialise perturbation
             vel = self.get_inivel(2)
     
             # Initialization of the populations at equilibrium with the given velocity.
             fin = self.get_equilibrium(1, vel)
                   
             ###### Setup: cylindrical obstacle and velocity inlet with perturbation ########
             # Creation of a mask with 1/0 values, defining the shape of the obstacle.
             obstacle = fromfunction(self.obstacle_fun, (self.nx,self.ny),dtype=int)
             
             ##### Run Simulation:
             for time in range(0,self.maxIter):
             
                # Right wall: outflow condition.
                fin[self.col3,-1,:] = fin[self.col3,-2,:] 

                # Compute macroscopic variables, density and velocity.
                rho, u = self.macroscopic(fin)

                # Left wall: inflow condition.
                u[:,0,:] = vel[:,0,:]
                rho[0,:] = 1/(1-u[0,0,:]) * (sum(fin[self.col2,0,:], axis=0) +
                                  2*sum(fin[self.col3,0,:], axis=0) )
                # Compute equilibrium.
                feq = self.get_equilibrium(rho, u)
                fin[[0,1,2],0,:] = feq[[0,1,2],0,:] + fin[[8,7,6],0,:] - feq[[8,7,6],0,:]
             
                # Collision step.
                fout = collision(self,self.omega,feq,fin,obstacle)

                # Streaming step.
                fin = streaming(self,fout,fin,self.v)
                
                # save data
                plot_field(self,time,u)

     
     
