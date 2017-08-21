#functions used in SCF loop of Hartree-Fock method

import numpy as np
from scipy import linalg
import math

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
                 

                #X = eVec * pow(eVal, -0.5) * eVec.transpose()
            
                test = np.linalg.qr(eVec)

                length = len(test[1])
                S = np.zeros([length,length])
                for b1 in range(len(test[1])):
                    for b2 in range(len(test[1])):

                        val = test[1][b1][b2]
                        if(b1 == b2):
                            if(val > 0):
                                S[b1][b2] = 1.0
                            else: 
                                S[b1][b2] = -1.0

                        else:
                            S[b1][b2] = 0.0


                
                Z = test[0] * S

		#build transformation matrix
		for b1 in range(nbf):
                    for b2 in range(nbf):
		            X[b1][b2] = Z[b1][b2] / math.sqrt(eVal[b1])
                 
		return X



		
