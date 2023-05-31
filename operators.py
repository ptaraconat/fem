import numpy as np 
import math 
from mesh_gen import * 

def stiffness(x,GPE):
	NPE = np.size(x,0)
	PD = np.size(x,1)
	
	K = np.zeros([NPE*PD,NPE*PD])
	
	coord = np.transpose(x)
	for i in range(1,NPE+1): 
		for j in range(1, NPE+1):
			k = np.zeros([PD,PD])
			# loop over Gauss points 
			for gp in range(1,GPE):
l				# init Jacobian
				J = np.zeros([PD,PD])
				# init gradient of interpolation functions
				grad = np.zeros([PD,NPE])

def GaussPoint(NPE, GPE, gp):
	

np.array([[0,0],
         [1,0],
         [1,2],
         [0.2]])
# Gauss point 
GP2 = 1. 


d1 = 1 
d2 = 1 
p = 4
m = 3
NL, EL = uniform_mesh(d1,d2,p,m,element_type = 'D2T3')
