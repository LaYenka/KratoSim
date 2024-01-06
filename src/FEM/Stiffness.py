#------------------------------------------------------------------------
# Compute Stiffnness Matrix
#
# -----------------------------------------------------------------------
import math
import numpy as np
from .shapefun import shapefun


class Stiffness():
	"""Class building up the stiffness matrix for 
		   FE Solution
	"""

	def __init__(self):
		self.name_shapef, self.shape_f = None, None
		# 1st Order quadrature - 2D
		# 2x2 Gauss Quadrature (4 Gauss points)
	    # q4 is 4x2
		self.q4 = np.array([[-1,-1],[1,-1],[-1,1],[1,1]]) / math.sqrt(3.0)
		# 1st Order quadrature - 1D
		# 1x2 Gauss Quadrature (2 Gauss points)
	    # q2 is 1x2
		self.q2 = np.array([[1.0,-1.0]]) / math.sqrt(3.0)
		#
		self.C = []
		self.C_beam = []
 

	"""Build Stiffnness matrix for FE solver
	"""
	def build_K_mat(self):
	# Make stiffness matrix
	# if N is the number of DOF, then K is NxN
		K = np.zeros((2*self.num_nodes, 2*self.num_nodes))    # square zero matrix
	#
		print('\n** Assemble stiffness matrix')
	# strain in an element: [strain] = B    U
	#                        3x1     = 3x8  8x1
	#
	# strain11 = B11 U1 + B12 U2 + B13 U3 + B14 U4 + B15 U5 + B16 U6 + B17 U7 + B18 U8
	#          = B11 u1          + B13 u1          + B15 u1          + B17 u1
	#          = dN1/dx u1       + dN2/dx u1       + dN3/dx u1       + dN4/dx u1
	# conn[0] is node numbers of the element
		for n_el in self.mesh_dict['Elem_num']:     # loop through each element
		# coordinates of each node in the element
		# e.g. shape = 4x2
		# for example:
		#    nodePts = [[0.0,   0.0],
		#               [0.033, 0.0],
		#               [0.033, 0.066],
		#               [0.0,   0.066]]
			c = self.conn[n_el-1]        # connectivtiy
			node_pts = self.nodes[c,:]    # extract node coordinates
			self.shape_f=getattr(shapefun,self.mesh_dict['shape_fun'][n_el-1])

			if len(node_pts) == 4:
				B = np.zeros((3,8))     # 
				Ke = np.zeros((8,8))    # element stiffness matrix is 8x8
				for q in self.q4:		# for each Gauss point
					# q is 1x2, N(xi,eta)
					# dN = self.gradshapefun(q)    # partial derivative of N wrt (xi,eta): 2x4
					N,dN = self.shape_f(q)         # shape function N and partial derivatives dN
					J  = np.dot(dN, node_pts).T     # Jacobian - J is 3x8
					dN = np.dot(np.linalg.inv(J), dN)    # partial derivative of N wrt (x,y): 2x4
					# assemble B matrix  [3x8]
					B[0,0::2] = dN[0,:]
					B[1,1::2] = dN[1,:]
					B[2,0::2] = dN[1,:]
					B[2,1::2] = dN[0,:]
					# element stiffness matrix (8X8)
					Ke += np.dot(np.dot(B.T,self.C),B) * np.linalg.det(J)

			if len(node_pts) == 2:
				B = np.zeros((1,4))     # 
				Rot = np.identity(4)
				Ke = np.zeros((4,4))    # element stiffness matrix is 4x4
				# print('linear element to be implemented')
				for q in self.q2:		# for each Gauss point
					# q is 1x2, N(xi,eta)
					# dN = self.gradshapefun(q)    # partial derivative of N wrt (xi): 1x4
					N,dN = self.shape_f(q)         # N and partial derivatives dN
					J  = np.dot(dN[0::2].T, node_pts).T     # Jacobian - J is 4x4
					L = np.linalg.norm(node_pts[1,:]-node_pts[0,:])
					# print(L,q)
					# integral factor scaling from [-1,1] -> Gauss formula to [0,1] reference element
					fact_int_gauss = 0.5
					# assemble B matrix  [1x4] 
					# ! check Sign of the matrix -> self.C_beam = 1
					B[0,0] = 6*q/L**2
					B[0,1] = (3*q-1)/L
					B[0,2] = -6*q/L**2
					B[0,3] = (3*q+1)/L
					# compute rotation matrxi to orientate element
					## cos = x/L 
					lox=(node_pts[1,0]-node_pts[0,0])/L
					## sin = y/L
					mox=(node_pts[1,1]-node_pts[0,1])/L
					# element stiffness matrix - bending (4X4)
					Ke += fact_int_gauss*np.dot(np.dot(B.T,self.C_beam),B)
					# print(Ke)
					# print(B,J)
					# rotation matrix: align element to the global coordinate system
					Rot[0,0] = lox
					Rot[0,1] = 0
					Rot[1,0] = 0
					Rot[1,1] = 1
					Rot[2,2] = lox
					Rot[3,3] = 1
					#  Rotate Global matrix
					Ke = np.dot(np.dot(Rot.T,Ke),Rot)



		# Scatter operation = assign to global matrix	
			for i,I in enumerate(c):
				for j,J in enumerate(c):
					K[2*I,2*J]     += Ke[2*i,2*j]
					K[2*I+1,2*J]   += Ke[2*i+1,2*j]
					K[2*I+1,2*J+1] += Ke[2*i+1,2*j+1]
					K[2*I,2*J+1]   += Ke[2*i,2*j+1]
		return K