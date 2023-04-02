# flow around cylinder


from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm

from core.LBM.LBM import LBM


if __name__== "__main__" :
###### Flow definition #########################################################
 maxIter = 200  # Total number of time iterations.
 Re = 40         # Reynolds number.
 print("Re",Re)
 nx, ny = 420, 180 # Numer of lattice nodes.
 ly = ny-1         # Height of the domain in lattice units.
 cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
 uLB     = 0.04                  # Velocity in lattice units.
 nulb    = uLB*r/Re;             # Viscoscity in lattice units.
 omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
 
 
 Sim = LBM()
 
 Sim.maxIter = 200000                # Total number of time iterations.
 Sim.Re = 50                         # Reynolds number.
 Sim.nx, Sim.ny = 420, 180          # Numer of lattice nodes.
 Sim.ly = ny-1                       # Height of the domain in lattice units.
 Sim.cx, Sim.cy, Sim.r = Sim.nx//4, Sim.ny//2, Sim.ny//9 # Coordinates of the cylinder.
 Sim.uLB     = 0.04                  # Velocity in lattice units.
 Sim.nulb    = Sim.uLB*r/Sim.Re;             # Viscoscity in lattice units.
 Sim.omega = 1 / (3*Sim.nulb+0.5);  # Relaxation parameter.
 
 Sim.run()
