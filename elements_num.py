import numpy as np 

def basis_interp_func_d2t3_1(xi,eta) : 
	return 1 - xi - eta 

def basis_interp_func_d2t3_2(xi,eta) : 
	return xi 

def basis_interp_func_d2t3_3(xi,eta) : 
	return eta 
	
def grad_basis_interp_func_d2t3_1(xi = None,eta = None) :
	return np.array([-1, -1])
	
def grad_basis_interp_func_d2t3_2(xi = None,eta = None) :
	return np.array([1, 0])

def grad_basis_interp_func_d2t3_3(xi = None,eta = None) :
	return np.array([0, 1])

class element() : 
	
	def __init__(self,type):
		self.type = type 
	
	def set_coordinates(self,coordinates_table):
		self.coordinates = coordinates_table
	
	def calc_hpp_rigidity_matrix(self,poisson = 0.,elastic_mod = 4.,thikness = 1.):
		# Get Material law 
		material_relation = self.calc_material_mat(elastic_mod,poisson,thikness)
		# Loop for gauss integration
		element_k = np.zeros([6,6])
		for i in range(len(self.gauss_weights)): 
			weight = self.gauss_weights[i]
			xi = self.gauss_points[i][0]
			eta = self.gauss_points[i][1]
			# Calc local sym part of gradient matrix 
			self.set_jacobian(xi,eta)
			B = self.calc_sym_grad_mat(xi,eta)
			# Gauss Integration 
			add_ = self.det_jac* np.dot( np.dot(np.transpose(B),material_relation), B)
			element_k = element_k + weight*add_ 
		#
		return element_k
		
	def calc_strain(self,nodal_displacements):
		# Calc B matrix / Sym Gradient matrix 
		xi, eta = self.middle_point
		B = self.calc_sym_grad_mat(xi, eta)
		strain = np.dot(B,nodal_displacements)
		return strain 
		
	def calc_stress(self,strain,young_modulus,poisson_modulus,thikness) : 
		# 
		material_matrix = self.calc_material_mat(young_modulus,poisson_modulus,thikness)
		stress = np.dot(material_matrix,strain)
		return stress 
		
class element_d2t3(element) : 
	
	def __init__(self):
		self.type = 'D2T3'
		self.dimension = 2
			
		self.interp1 = basis_interp_func_d2t3_1
		self.interp2 = basis_interp_func_d2t3_2
		self.interp3 = basis_interp_func_d2t3_3
		
		self.grad1 = grad_basis_interp_func_d2t3_1
		self.grad2 = grad_basis_interp_func_d2t3_2
		self.grad3 = grad_basis_interp_func_d2t3_3
		
		self.gauss_points = [[1/6, 1/6],[2/3, 1/6],[1/6, 2/3]] 
		self.gauss_weights = [1/6, 1/6, 1/6]
		
		self.middle_point = (1/3,1/3)
		
	def calc_correspondencies(self,element): 
		el = element
		correspondance = [[0, 2*el[0]],
		                  [1, 2*el[0]+1],
		                  [2, 2*el[1]],
		                  [3, 2*el[1]+1],
		                  [4, 2*el[2]],
		                  [5, 2*el[2]+1]]
		correspondance = np.asarray(correspondance)
		return correspondance
			
		
	def calc_material_mat(self,young_modulus,poisson_modulus,thikness):
		pois = poisson_modulus
		E = young_modulus
		thikness = thikness
		fac_tmp = (thikness*E)/((1+pois**2))
		material_relation = fac_tmp*np.array([[1, pois, 0],[pois, 1., 0],[0., 0., (1-pois)/2.]])
		return material_relation
		
	def calc_basis_interpolation(self,xi,eta): 
		interp_mat = [self.interp1(xi,eta),self.interp2(xi,eta),self.interp3(xi,eta)]
		interp_mat = np.asarray(interp_mat)#np.asarray([interp_mat,interp_mat])
		return interp_mat
	
	def calc_local_interpolation(self,xi,eta):
		basis_interp = self.calc_basis_interpolation(xi,eta)
		self.set
		
	def transform_basis_to_local(self,xi,eta):
		interp = self.calc_basis_interpolation(xi,eta)
		node_coord = self.coordinates
		local_coordinates = np.dot(interp,node_coord)
		return local_coordinates
	
	def calc_jacobian(self,xi,eta): 
		grad_n = self.calc_basis_gradient(xi,eta)
		jacobian = np.dot(grad_n,self.coordinates)
		return jacobian
		
	def set_jacobian(self,xi,eta):
		jacobian = self.calc_jacobian(xi,eta)
		self.jacobian = jacobian
		self.det_jac = np.linalg.det(self.jacobian)
	
	def calc_basis_gradient(self,xi,eta):
		grad_n = [self.grad1(xi,eta).tolist(),self.grad2(xi,eta).tolist(),self.grad3(xi,eta).tolist()]
		grad_n = np.asarray(np.transpose(grad_n))
		return grad_n
		
	def calc_local_gradient(self,xi,eta):
		jacobian = self.calc_jacobian(xi,eta)
		inv_jacobian = np.linalg.inv(jacobian)
		basis_gradient = self.calc_basis_gradient(xi,eta) 
		local_gradient = np.dot(inv_jacobian,basis_gradient)
		return local_gradient
		
	def calc_sym_grad_mat(self,xi,eta): 
		mat_tmp = self.calc_local_gradient(xi,eta)
		dn1dx = mat_tmp[0,0]
		dn1dy = mat_tmp[1,0]
		dn2dx = mat_tmp[0,1]
		dn2dy = mat_tmp[1,1]
		dn3dx = mat_tmp[0,2]
		dn3dy = mat_tmp[1,2]
		B = np.array([[dn1dx, 0, dn2dx, 0, dn3dx, 0],[0, dn1dy, 0, dn2dy, 0, dn3dy],[dn1dy, dn1dx, dn2dy, dn2dx, dn3dy, dn3dx]])
		return B 
			
	def calc_bc_load(self,force_density,thikness = 1.) : 
		# Loop for gauss integration
		element_load = np.zeros([6,1])
		flatten_load = force_density.flatten()
		for i in range(len(self.gauss_weights)): 
			weight = self.gauss_weights[i]
			xi = self.gauss_points[i][0]
			eta = self.gauss_points[i][1]
			# get interpolation mat 
			N1 = self.interp1(xi,eta)
			N2 = self.interp2(xi,eta)
			N3 = self.interp3(xi,eta)
			self.set_jacobian(xi,eta)
			interp_mat = np.array([ [N1, 0, N2, 0, N3, 0], 
			                       [0, N1, 0, N2, 0, N3] ])
			#print(flatten_load)
			#print(interp_mat)
			add_ = self.det_jac * np.dot(interp_mat,flatten_load)
			#print(add_)
			element_load = element_load + weight*add_ 
		return element_load
			
			
		
			
	#def calc_rigidity_matrix(self): 
	#	stiffness = np.zeros([3,3])
	#	for i in range(3):
	#		for j in range(3):
	#			gradi = self.gradient(i)
	#			gradj = self.gradient(j)
	#			jact_jac = np.transpose(self.jacobian) * self.jacobian 
	#			stiffness[i,j] = np.dot(gradi,np.dot(jact_jac,gradj))
	#	return self.det_jac*stiffness
	
	
		
		
		

coord = np.array([0,0])+np.array([[0, 0],
          [1, 0],
          [0, 1]])
   
#print(coord)       

element = element_d2t3()
element.set_coordinates(coord)


#element.calc_local_gradient(0,0)

xi = 0.1
eta = 0.3 
print('set basis coordinates : ',xi,' ',eta)
loc_coord = element.transform_basis_to_local(xi,eta)
print('local coordinates')
print(loc_coord)
loc_grad = element.calc_local_gradient(xi,eta)
print('local gradient given basis coordinates')
print(loc_grad)
basis_grad = element.calc_basis_gradient(xi,eta)
print('basis gradient given basis coordinates')
print(basis_grad)
kmat = element.calc_hpp_rigidity_matrix()
print('local stiffness matrix ')
print(kmat)



#print(element.calc_rigidity_matrix())
#print(element.calc_hpp_rigidity_matrix())

		
		 
