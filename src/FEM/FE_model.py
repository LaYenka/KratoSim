#-------------------------------------------------------
#  Generate FE Model
# 
#-------------------------------------------------------

import numpy as np
import math
import sys

from .shapefun import *
from .ReadMesh import *
from .Stiffness import *
from .FE_forcing import *
from .FE_PostProcessing import *

class FE_model(read_mesh,shapefun,FE_forcing,Stiffness,FE_PostProcessing):

  def __init__(self,input_file):
    self.nodes = []
    self.num_nodes = []
    self.conn = []
    self.boundary = []
    self.input_file = input_file
    self.name_shapef, self.shape_f = None, None
    self.mesh_dict = {'Elem_type':[], 'Elem_num':[],'shape_fun':[]}

    # Gauss points
    self.q4 = np.array([[-1,-1],[1,-1],[-1,1],[1,1]]) / math.sqrt(3.0)
    self.q2 = np.array([[-1.0],[1.0]]) / math.sqrt(3.0)         # 0.5 factor scaling the interval [-1,1] to [0,1] reflected in the integral

    
    # overall variables
    self. u = []           # displacement
    self.node_strain = []  # strain
    self.node_stress = []  # stress
    # material properties - Young Modulus E, Poisson Ratio v, Tensor Matrix C
    self.E=[]
    self.v=[]
    self.C = []
    self.C_beam = []

    # post-processing
    self.plot_type = []

	
  def read_mesh(self): 
    self.read_inp_file(self.input_file,self.nodes,self.conn,self.boundary)
    self.nodes = np.array(self.nodes)
    self.num_nodes = len(self.nodes)
    print('   number of nodes:', len(self.nodes))
    print('   number of elements:', len(self.conn))
    print('   number of displacement boundary conditions:', len(self.boundary))

  def FE_Solve(self):
    """Build up matrices, frocing terms and solve"""
         
     # read model mesh
    self.read_mesh()

     # assemble stiffness matrix
    K = self.build_K_mat()
         
	  # build r.h.s
    f, K = self.build_forcing(self.num_nodes,self.boundary,K)

    # solve system
    ###############################
    print('\n** Solve linear system: Ku = f')	# [K] = 2N x 2N, [f] = 2N x 1, [u] = 2N x 1
    self.u = np.linalg.solve(K, f)
	  ###############################
    print('\n** finish solving the system')

