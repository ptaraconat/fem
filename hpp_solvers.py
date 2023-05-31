from mesh_gen import * 
from elements_num import * 


class MecaSolver():
	
	def __init__(self,dim,element_type,poisson=1,young=1):
		self.dim = dim 
		self.element_type = element_type
		self.poisson = poisson
		self.young = young
		if element_type == 'D2T3' : 
			self.element = element_d2t3()
			self.dim = 2 
			self.nelem_nodes = 3
			self.thikness = 1
		
	def set_grid(self,nodes, connectivity):
		self.nodes = nodes 
		self.connectivity = connectivity-1
		self.nnodes = np.size(nodes,0)
		self.nelem = np.size(connectivity,0)
		
	def calc_stifness_matrix(self) : 
		stiff_matrices = np.zeros((self.dim*self.nelem_nodes,self.dim*self.nelem_nodes,self.nelem))
		i = 0 
		for el in self.connectivity.tolist():
			nodes = self.nodes[el,:]
			self.element.set_coordinates(nodes)
			stiff_matrices[:,:,i] = self.element.calc_hpp_rigidity_matrix(poisson = self.poisson,elastic_mod = self.young, thikness = self.thikness)
			i = i + 1 
		print(stiff_matrices)
		print(np.shape(stiff_matrices))
			


####### Define a mesh ########
d1 = 1 
d2 = 1 
p = 4
m = 3
NL, EL = uniform_mesh(d1,d2,p,m,element_type = 'D2T3')
solver = MecaSolver(2,'D2T3')
solver.set_grid(NL,EL)
solver.calc_stifness_matrix()
