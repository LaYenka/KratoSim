import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from src.FEM.FE_model import FE_model

if __name__=="__main__":
    # parameters
    input_file = "Dogbone_Tension-1D.input"
    
    # setup shape function
    # shape_fun = "linear_shape_2D"

    # Define material properties
    ###############################
	# Plane-strain material tangent (see Bathe p. 194)
	# C is 3x3
    E = 100.0
    v = 0.3
    C = E/(1.0+v)/(1.0-2.0*v) * np.array([[1.0-v, v, 0.0], [v, 1.0-v, 0.0], [0.0, 0.0, 0.5-v]])
    C_beam = 0

    # instatiate class object
    Sim = FE_model(input_file)
    # including __call__ objects
    # Sim()

    # setup material matrix
    Sim.C = C
    Sim.C_beam = C_beam
    
    # solve FE model
    Sim.FE_Solve()

    # postprocess data
    Sim.plot_type = 'e11'
    Sim.stress_strain()
