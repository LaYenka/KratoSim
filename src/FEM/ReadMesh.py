#-------------------------------------------------------------
#
# Read Mesh
#
#-------------------------------------------------------------

import numpy as np
import math
import sys

class read_mesh:

	def __init__(self) -> None:
		self.mesh_dict = {'Elem_type':[], 'Elem_num':[],'shape_fun':[]}

	#def __call__(self, *a, **k):
      ##
    #if self.shape_f is None:
    #  modn, funcn = self.name_shapef.rsplit('.', 1)
      #if modn not in sys.modules:
      #  __import__(modn)
      #self.f = getattr(sys.modules[modn],
      #                 funcn)
    # definition of the shape function
    #	self.shape_f = getattr(shapefun,self.name_shapef)
    #print(self.shape_f)
    # self.shape_f(*a, **k) 


	def read_inp_file(self,inpFileName, nodes, conn, boundary):
		print('\n** Read input file')
		inpFile = open(inpFileName, 'r')
		lines = inpFile.readlines()
		inpFile.close()
		state = 0
		for line in lines:
			line = line.strip()
			if len(line) <= 0: continue
			if line[0] == '*':
				state = 0
			if line.lower() == "*node":
				state = 1
				continue
			if line.lower() == "*element":
				state = 2
				continue
			if line.lower() == "*boundary":
				state = 3
				continue
			if state == 0:
				continue
			if state == 1:
				# read nodes
				values = line.split(",")
				if len(values) != 3:
					local_error("A node definition needs 3 values")
				nodeNr = int(values[0]) - 1  # zero indexed
				xx = float(values[1])
				yy = float(values[2])
				nodes.append([xx,yy])   # assume the nodes are ordered 1, 2, 3...
				continue
			if state == 2:
				# read elements
				values = line.split(",")
				if (len(values) != 5) and (len(values) != 3) :
					local_error("An element definition needs 5 or 3 values")
				elemNr = int(values[0])

				# populate dictionary
				self.mesh_dict["Elem_type"].append(len(values)) # populate dictionary: Element_type
				self.mesh_dict["Elem_num"].append(elemNr)       # populate dictionary: Element_number

				# read 2D element
				if len(values) == 5:
					
					self.mesh_dict["shape_fun"].append('linear_2D')  # populate dictionary: shape function name

					n1 = int(values[1]) - 1  # zero indexed
					n2 = int(values[2]) - 1
					n3 = int(values[3]) - 1
					n4 = int(values[4]) - 1
					#conn.append([n1, n2, n3, n4]) # assume elements ordered 1, 2, 3
					conn.append([n1, n4, n3, n2]) # assume elements ordered 1, 2, 3
					continue

				# read 1D elements
				if len(values) == 3:

					self.mesh_dict["shape_fun"].append("beam_1D")  # populate dictionary: shape function name

					n1 = int(values[1]) - 1  # zero indexed
					n2 = int(values[2]) - 1
					conn.append([n1, n2]) # assume elements ordered 1, 2, 3
					continue
			if state == 3:
				# read displacement boundary conditions
				values = line.split(",")
				if len(values) != 4:
					local_error("A displacement boundary condition needs 4 values")
				nodeNr = int(values[0]) - 1  # zero indexed
				dof1 = int(values[1])
				dof2 = int(values[2])
				val = float(values[3])
				if dof1 == 1:
					boundary.append([nodeNr,1,val])
				if dof2 == 2:
					boundary.append([nodeNr,2,val])
				continue

