import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from src.FEM.FE_model import FE_model

if __name__=="__main__":
    # parameters
    input_file = "Beam_2D.input" #"Dogbone_Tension-1D.input"
    
    # setup shape function
    # shape_fun = "linear_shape_2D"

    # Define material properties
    ###############################
	# Plane-strain material tangent (see Bathe p. 194)
	# C is 3x3
    E = 68.0*10**9  #Pa 
    v = 0.3
    C = E/(1.0+v)/(1.0-2.0*v) * np.array([[1.0-v, v, 0.0], [v, 1.0-v, 0.0], [0.0, 0.0, 0.5-v]])
    I_beam = (0.08*0.08**3)/3.0
    C_beam = E*I_beam


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
    Sim.plot_type = 'u1'
    Sim.stress_strain()
