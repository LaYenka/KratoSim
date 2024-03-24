# Simulation settings

def SimSettings(self):
 maxIter = 200000  # Total number of time iterations.
 Re = Ree         # Reynolds number.
 print("Re",Re)
 nx, ny = 420, 180 # Numer of lattice nodes.
 ly = ny-1         # Height of the domain in lattice units.
 cx, cy, r = nx//4, ny//2, ny//9 # Coordinates of the cylinder.
 uLB     = 0.04                  # Velocity in lattice units.
 nulb    = uLB*r/Re;             # Viscoscity in lattice units.
 omega = 1 / (3*nulb+0.5);    # Relaxation parameter.
