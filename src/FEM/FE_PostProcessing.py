import numpy as np
from matplotlib import pyplot as plt
from .shapefun import shapefun



class FE_PostProcessing:

	def __init__(self,name):
		self.node_strain = []
		self.node_stress = []
		self.nodes = []
		self.u =[]
		self.conn = []
		self.name_shapef, self.shape_f = name, None
		self.plot_type = 'e11'


	def stress_strain(self):
		print('\n** Post process the data')
	# (pre-allocate space for nodal stress and strain)
		for ni in range(len(self.nodes)):
			self.node_strain.append([0.0, 0.0, 0.0])
			self.node_stress.append([0.0, 0.0, 0.0])
		node_strain = np.array(self.node_strain)
		node_stress = np.array(self.node_stress)

		print(f'   min displacements: u1={min(self.u[0::2]):.4g}, u2={min(self.u[1::2]):.4g}')
		print(f'   max displacements: u1={max(self.u[0::2]):.4g}, u2={max(self.u[1::2]):.4g}')
		emin = np.array([ 9.0e9,  9.0e9,  9.0e9])
		emax = np.array([-9.0e9, -9.0e9, -9.0e9])
		smin = np.array([ 9.0e9,  9.0e9,  9.0e9])
		smax = np.array([-9.0e9, -9.0e9, -9.0e9])

		for n_el in self.mesh_dict['Elem_num']:     # loop through each element
			                                        # for each element (conn is Nx4)
										 # c is like [2,5,22,53]			
			c = self.conn[n_el-1]        # connectivtiy						
			self.shape_f=getattr(shapefun,self.mesh_dict['shape_fun'][n_el-1])
			nodePts = self.nodes[c,:]			# 4x2, eg: [[1.1,0.2], [1.2,0.3], [1.3,0.4], [1.4, 0.5]]
		
			if len(nodePts) == 4:
				B = np.zeros((3,8))     # 
				for q in self.q4:					# for each integration pt, eg: [-0.7,-0.7]
					N,dN = self.shape_f(q)              # 2x4
					J  = np.dot(dN, nodePts).T			# 2x2
					dN = np.dot(np.linalg.inv(J), dN)	# 2x4
					B[0,0::2] = dN[0,:]					# 3x8
					B[1,1::2] = dN[1,:]
					B[2,0::2] = dN[1,:]
					B[2,1::2] = dN[0,:]

					UU = np.zeros((8,1))				# 8x1
					UU[0] = self.u[2*c[0]]
					UU[1] = self.u[2*c[0] + 1]
					UU[2] = self.u[2*c[1]]
					UU[3] = self.u[2*c[1] + 1]
					UU[4] = self.u[2*c[2]]
					UU[5] = self.u[2*c[2] + 1]
					UU[6] = self.u[2*c[3]]
					UU[7] = self.u[2*c[3] + 1]
					# get the strain and stress at the integration point
					strain = B @ UU		# (B is 3x8) (UU is 8x1) 		=> (strain is 3x1)
					stress = self.C @ strain	# (C is 3x3) (strain is 3x1) 	=> (stress is 3x1)
					emin[0] = min(emin[0], strain[0][0])
					emin[1] = min(emin[1], strain[1][0])
					emin[2] = min(emin[2], strain[2][0])
					emax[0] = max(emax[0], strain[0][0])
					emax[1] = max(emax[1], strain[1][0])
					emax[2] = max(emax[2], strain[2][0])

					node_strain[c[0]][:] = strain.T[0]
					node_strain[c[1]][:] = strain.T[0]
					node_strain[c[2]][:] = strain.T[0]
					node_strain[c[3]][:] = strain.T[0]
					node_stress[c[0]][:] = stress.T[0]
					node_stress[c[1]][:] = stress.T[0]
					node_stress[c[2]][:] = stress.T[0]
					node_stress[c[3]][:] = stress.T[0]
					smax[0] = max(smax[0], stress[0][0])
					smax[1] = max(smax[1], stress[1][0])
					smax[2] = max(smax[2], stress[2][0])
					smin[0] = min(smin[0], stress[0][0])
					smin[1] = min(smin[1], stress[1][0])
					smin[2] = min(smin[2], stress[2][0])

			if len(nodePts) == 2:
					B = np.zeros((1,4))     #
					for q in self.q2:		# for each Gauss point
						# q is 1x2, N(xi,eta)
						# dN = self.gradshapefun(q)    # partial derivative of N wrt (xi): 1x4
						N,dN = self.shape_f(q)         # N and partial derivatives dN
						J  = np.dot(dN[0::2].T, nodePts).T     # Jacobian - J is 1
						L = np.linalg.norm(nodePts[0,:]-nodePts[1,:])
						# assemble B matrix  [1x4]
						B[0,0] = -6*q/L**2
						B[0,1] = (3*q-1)/L
						B[0,2] = -6*q/L**2
						B[0,3] = (3*q+1)/L 

						UU = np.zeros((4,1))				# 4x1
						UU[0] = self.u[2*c[0]]
						UU[1] = self.u[2*c[0] + 1]
						UU[2] = self.u[2*c[1]]
						UU[3] = self.u[2*c[1] + 1]

											# get the strain and stress at the integration point
						strain = B @ UU		# (B is 1x4) (UU is 4x1) 		=> (strain is 1x1)
						stress = self.C_beam * strain	# (C is 1x1) (strain is 1x1) 	=> (stress is 1x1)
						print(strain)
						emin[0] = min(emin[0], strain[0])
						emax[0] = max(emax[0], strain[0])

						node_strain[c[0]][:] = strain.T[0]
						node_strain[c[1]][:] = strain.T[0]
						smax[0] = max(smax[0], stress[0][0])

		print(f'   min strains: e11={emin[0]:.4g}, e22={emin[1]:.4g}, e12={emin[2]:.4g}')
		print(f'   max strains: e11={emax[0]:.4g}, e22={emax[1]:.4g}, e12={emax[2]:.4g}')
		print(f'   min stress:  s11={smin[0]:.4g}, s22={smin[1]:.4g}, s12={smin[2]:.4g}')
		print(f'   max stress:  s11={smax[0]:.4g}, s22={smax[1]:.4g}, s12={smax[2]:.4g}')
		
		
		###############################
		print('\n** Plot displacement')
		xvec = []
		yvec = []
		res  = []
		for ni,pt in enumerate(self.nodes):
			xvec.append(pt[0] + self.u[2*ni])
			yvec.append(pt[1] + self.u[2*ni+1])
			if self.plot_type=='u1':  res.append(self.u[2*ni])			# x-disp
			if self.plot_type=='u2':  res.append(self.u[2*ni+1])		# y-disp
			if self.plot_type=='s11': res.append(node_stress[ni])		# s11
			if self.plot_type=='s22': res.append(node_stress[ni])		# s22
			if self.plot_type=='s12': res.append(node_stress[ni])		# s12
			if self.plot_type=='e11': res.append(node_strain[ni])		# e11
			if self.plot_type=='e22': res.append(node_strain[ni])		# e22
			if self.plot_type=='e12': res.append(node_strain[ni])		# e12
		tri = []
		if len(nodePts) == 4:
			for c in self.conn:
				if len(c) == 4:
					tri.append( [c[0], c[1], c[2]] )
					tri.append( [c[0], c[2], c[3]] )
			t = plt.tricontourf(xvec, yvec, res, triangles=tri, levels=14, cmap=plt.cm.jet)
		 	#plt.scatter(xvec, yvec, marker='o', c='b', s=0.5) # (plot the nodes)
			plt.grid()
			plt.colorbar(t)
			plt.title(self.plot_type)
			plt.axis('equal')
			plt.show()
			print('Done.')
		bi = []
		if len(nodePts) == 2:
			for c in self.conn:
				if len(c) == 2:
					bi.append( [c[0], c[1]])
			print('plotting beam')
			fig, ax = plt.subplots()
			ax.plot(xvec, res)
		 	#plt.scatter(xvec, yvec, marker='o', c='b', s=0.5) # (plot the nodes)
			ax.grid()
			plt.title(self.plot_type)
			ax.set_ylim(min(res),max(res))
			plt.show()
			print('Done.')