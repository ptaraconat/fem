from mesh_gen import * 
from elements_num import * 


class MecaSolver():
	
	def __init__(self,dim,element_type,poisson=1,young=1,thikness = 1):
		self.dim = dim 
		self.element_type = element_type
		self.poisson = poisson
		self.young = young
		if element_type == 'D2T3' : 
			self.element = element_d2t3()
			self.dim = 2 
			self.nelem_nodes = 3
			self.thikness = thikness
			
		
	def set_grid(self,nodes, connectivity):
		self.nodes = nodes 
		self.connectivity = connectivity-1
		self.nnodes = np.size(nodes,0)
		self.nelem = np.size(connectivity,0)
		self.force = np.zeros([self.dim*self.nnodes,1])
	
	def assemble_stiffness(self):
		assembled_matrix = np.zeros([self.dim*self.nnodes,self.dim*self.nnodes])
		for el in self.connectivity.tolist(): 
			# Calc element stiffness
			nodes = self.nodes[el,:]
			self.element.set_coordinates(nodes)
			element_mat = self.element.calc_hpp_rigidity_matrix(poisson = self.poisson,elastic_mod = self.young, thikness = self.thikness)
			# 
			correspondance = [[0, 2*el[0]],
			                  [1, 2*el[0]+1],
			                  [2, 2*el[1]],
			                  [3, 2*el[1]+1],
			                  [4, 2*el[2]],
			                  [5, 2*el[2]+1]]
			correspondance = np.asarray(correspondance)
			#print(correspondance)
			#
			for i in range(len(correspondance)): 
				for j in range(len(correspondance)):
					loci = correspondance[i,0]
					locj = correspondance[j,0]
					ass_i = correspondance[i,1]
					ass_j = correspondance[j,1]
					#print(loci,locj,ass_i,ass_j)
					assembled_matrix[ass_i,ass_j] = assembled_matrix[ass_i,ass_j] + element_mat[loci,locj]
		return assembled_matrix
	
	def set_stiffness(self): 
		stiffness = self.assemble_stiffness()
		self.stiffness = stiffness 
	
	def add_nodal_force(self,node_id,force) : 
		#print('Hey im gona add some forces')
		self.force[2*node_id,0] += force[0]
		self.force[2*node_id+1,0] += force[1]
		
	def add_boundary_constant_force(self,bnd,force):
		for node_id in bnd : 
			self.add_nodal_force(node_id,force)
	
	def set_nodal_displacement(self,node_id,dispx = 0., dispy =0., dispz = None):
		if dispx != None : 
			self.stiffness[2*node_id,:] = 0
			self.stiffness[2*node_id,2*node_id] = 1.
			self.force[2*node_id,0] = dispx
		if dispy != None : 
			self.stiffness[2*node_id+1,:] = 0
			self.stiffness[2*node_id+1,2*node_id+1] = 1.
			self.force[2*node_id+1,0] = dispy
	
	def set_boundary_nodal_displacement(self,bnd,dispx = 0., dispy =0., dispz = None):
		for node_id in bnd : 
			self.set_nodal_displacement(node_id,dispx=dispx,dispy=dispy,dispz=dispz)
	
	def solve(self): 
		#return np.linalg.solve(self.stiffness,self.force)
		return np.dot(np.linalg.pinv(self.stiffness),self.force)
	
	def display_displacement(self,disp,amplify = 1.) : 
		new_loc = self.nodes
		disp_tmp = disp.reshape((self.nnodes,self.dim))
		new_loc = new_loc + amplify*disp_tmp
		plt.plot(self.nodes[:,0],self.nodes[:,1],'r+')
		plt.plot(new_loc[:,0],new_loc[:,1],'ko')
		plt.show()
		

def set_mesh_bounds(nodes_table,xmin=0.,xmax=1.,ymin=0.,ymax = 1.) : 
	bound1 = []
	bound2 = []
	bound3 = []
	bound4 = []
	for i in range(np.size(nodes_table,0)): 
		coords = nodes_table[i,:]
		if coords[0] == xmin : 
			bound1.append(i)
		if coords[0] == xmax : 
			bound2.append(i)
		if coords[1] == ymin : 
			bound3.append(i)
		if coords[1] == ymax : 
			bound4.append(i)
	dico = {} 
	dico['bnd1'] = bound1
	dico['bnd2'] = bound2
	dico['bnd3'] = bound3
	dico['bnd4'] = bound4
	return dico 
		

####### Define a mesh ########
d1 = 10 
d2 = 5
p = 20
m = 20
NL, EL = uniform_mesh(d1,d2,p,m,element_type = 'D2T3')
plt.close()
dico_bnd = set_mesh_bounds(NL,xmax = d1, ymax = d2)
#######
f = 10
solver = MecaSolver(2,'D2T3',poisson=0.8,young=40000,thikness = 2)
solver.set_grid(NL,EL)
solver.set_stiffness()
solver.add_boundary_constant_force(dico_bnd['bnd2'],[0,f])
solver.set_boundary_nodal_displacement(dico_bnd['bnd1'],dispy = 0. )
displacement = solver.solve()
solver.display_displacement(displacement,amplify = 1)





