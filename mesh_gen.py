import numpy as np
import matplotlib.pyplot as plt 

def uniform_mesh(d1,d2,p,m,element_type = 'D2T3'):
	PD = 2 # problem dimension
	NPE = 4 # number of npde êr element (true for D2Q4)
	# Define the domain 
	q = np.array([[0,0],[d1,0],[0,d2],[d1,d2]])
	# calculate number of nodes 
	NoN = (p+1)*(m+1)
	# calculate number of elements 
	NoE = p*m
	
	### Nodes table ### 
	NL = np.zeros([NoN,PD])
	a = (q[1,0]-q[0,0])/p#(d1/p) # increment in horizontal axis (delta x )
	b = (q[2,1]-q[0,1])/m# d2/m  # increment in vertical axis (delta y)
	# Nodes Loop 
	n = 0 #allow to go through rows in NL 
	for i in range(1 , m +2):
		for j in range(1 , p+2) :
			NL[n,0] = q[0,0] + (j-1)*a
			NL[n,1] = q[0,0] + (i-1)*b
			n += 1
			
	### Element table ### 
	EL = np.zeros([NoE,NPE],dtype = int)
	# Elements loop 
	for i in range(1,m+1):
		for j in range(1,p+1):
			if j == 1 :
				EL[(i-1)*p+j-1,0] = (i-1)*(p+1) + j 
				EL[(i-1)*p+j-1,1] = EL[(i-1)*p+j-1,0] + 1 
				EL[(i-1)*p+j-1,3] = EL[(i-1)*p+j-1,0] + (p+1)
				EL[(i-1)*p+j-1,2] = EL[(i-1)*p+j-1,3] + 1 
			else : 
				EL[(i-1)*p+j-1,0] = EL[(i-1)*p+j-2,1]
				EL[(i-1)*p+j-1,3] = EL[(i-1)*p+j-2,2]
				EL[(i-1)*p+j-1,1] = EL[(i-1)*p+j-1,0] + 1
				EL[(i-1)*p+j-1,2] = EL[(i-1)*p+j-1,3] + 1
	
	if element_type == 'D2T3': 
		NPE_new = 3 
		NoE_new = NoE*2
		EL_new = np.zeros([NoE_new,NPE_new],dtype = int)
		for i in range(1,NoE+1):
			# First triangle within the rectangle
			EL_new[2*(i-1), 0] = EL[i-1,0]
			EL_new[2*(i-1), 1] = EL[i-1,1]
			EL_new[2*(i-1), 2] = EL[i-1,2]
			# Second triangle within the rectangle 
			EL_new[2*(i-1)+1, 0] = EL[i-1,0]
			EL_new[2*(i-1)+1, 1] = EL[i-1,2]
			EL_new[2*(i-1)+1, 2] = EL[i-1,3]
		EL = EL_new
		
	display = True 
	if display : 
		print(EL)
		NoN = np.size(NL,0)
		NoE = np.size(EL,0)
		plt.figure(1)
		count = 1 
		for i in range(0,NoN):
			plt.annotate(count, xy =(NL[i,0],NL[i,1]))
			count +=1
		if element_type == 'D2Q4' : 
			count2 = 1 
			for j in range(0,NoE):
				plt.annotate(count2, xy=( (NL[EL[j,0]-1,0]+NL[EL[j,1]-1,0]+NL[EL[j,2]-1,0]+NL[EL[j,3]-1,0])/4,(NL[EL[j,0]-1,1]+NL[EL[j,1]-1,1]+NL[EL[j,2]-1,1]+NL[EL[j,3]-1,1])/4) )
				count2 +=1
			x0,y0 = NL[EL[:,0]-1,0] , NL[EL[:,0]-1,1]
			x1,y1 = NL[EL[:,1]-1,0] , NL[EL[:,1]-1,1]
			x2,y2 = NL[EL[:,2]-1,0] , NL[EL[:,2]-1,1]
			x3,y3 = NL[EL[:,3]-1,0] , NL[EL[:,3]-1,1]
			
			plt.plot(np.array([x0,x1]),np.array([y0,y1]),'red',linewidth = 3)
			plt.plot(np.array([x1,x2]),np.array([y1,y2]),'red',linewidth = 3)
			plt.plot(np.array([x2,x3]),np.array([y2,y3]),'red',linewidth = 3)
			plt.plot(np.array([x3,x0]),np.array([y3,y0]),'red',linewidth = 3)
			
		if element_type == 'D2T3' : 
			count2 = 1 
			for j in range(0,NoE):
				plt.annotate(count2, xy=( (NL[EL[j,0]-1,0]+NL[EL[j,1]-1,0]+NL[EL[j,2]-1,0])/3,(NL[EL[j,0]-1,1]+NL[EL[j,1]-1,1]+NL[EL[j,2]-1,1])/3) )
				count2 +=1
			x0,y0 = NL[EL[:,0]-1,0] , NL[EL[:,0]-1,1]
			x1,y1 = NL[EL[:,1]-1,0] , NL[EL[:,1]-1,1]
			x2,y2 = NL[EL[:,2]-1,0] , NL[EL[:,2]-1,1]

			
			plt.plot(np.array([x0,x1]),np.array([y0,y1]),'red',linewidth = 3)
			plt.plot(np.array([x1,x2]),np.array([y1,y2]),'red',linewidth = 3)
			plt.plot(np.array([x2,x0]),np.array([y2,y0]),'red',linewidth = 3)
		plt.show()
					
	return NL, EL 
				

