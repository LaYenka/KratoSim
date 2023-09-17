
import numpy as np

class FE_forcing():

	def build_forcing(self,num_nodes,boundary,K):
	###############################
	# Assign nodal forces and boundary conditions
	#    if N is the number of nodes, then f is 2xN
		f = np.zeros((2*num_nodes))          # initialize to 0 forces
	# How about displacement boundary conditions:
	#    [k11 k12 k13] [u1] = [f1]
	#    [k21 k22 k23] [u2]   [f2]
	#    [k31 k32 k33] [u3]   [f3]
	#
	#    if u3=x then
	#       [k11 k12 k13] [u1] = [f1]
	#       [k21 k22 k23] [u2]   [f2]
	#       [k31 k32 k33] [ x]   [f3]
	#   =>
	#       [k11 k12 k13] [u1] = [f1]
	#       [k21 k22 k23] [u2]   [f2]
	#       [  0   0   1] [u3]   [ x]
	#   the reaction force is
	#       f3 = [k31 k32 k33] * [u1 u2 u3]
		for i in range(len(boundary)):  # apply all boundary displacements
			nn  = boundary[i][0]
			dof = boundary[i][1]
			val = boundary[i][2]
			j = 2*nn
			if dof == 2: j = j + 1
			K[j,:] = 0.0
			K[j,j] = 1.0
			f[j] = val
	# return forcing terms f and modify stiffness matrix K accordingly
		return f, K
	###############################