#functions used in SCF loop of Hartree-Fock method

import numpy as np

#########################
class SCF:
	
#########################
	def getTransform(self, S):
		#builds matrix to transform operators from AO to MO basis
		#based off of canonical orthogonalization presented in "Modern Quantom Chemistry," on page 159, equation 3.170

		#number of basis functions used
		nbf = len(S)
	
		#init transformation matrix X
		X = np.zeros([nbf,nbf])

		#init unitary matrix U
		U = np.identity(nbf)
		
		#get overlap matrix eigen values and vectors
		eVal, eVec = np.linalg.eig(S)

		#build transformation matrix
		for b1 in range(nbf):
			for b2 in range(nbf):
				X[b1][b2] = U[b1][b2] / eVal[b2]

		return X 



		
