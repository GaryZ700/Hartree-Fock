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
                #uses equations present in step 4 from the following website as reference: http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3
		
                #number of basis functions used
		nbf = len(S)
	
		#init transformation matrix X
		X = np.asmatrix(np.zeros([nbf,nbf]))

		#get overlap matrix eigen values and vectors
		eVal, eVec = np.linalg.eigh(S)

                #get inverse square root of eigen values
                for basis in range(nbf):
                    eVal[basis] = eVal[basis]**-.5
                
                #sort eigen values and vectors
                print(eVal.argsort())
                print(eVec.argsort())

                #create diagonalized eigen value matrix
	        eValMatrix = self.eValMatrix(eVal)

                #build transformation matrix
                X = eVec * eValMatrix * eVec.transpose()

		return X

#########################
        def eValMatrix(self, eVal):
            #builds diagonolized eigen value matrix
           
            #get number of basis functions used
            nbf = len(eVal)

            #build zero matrix 
            matrix = np.asmatrix(np.zeros([nbf,nbf]))
            
            #construct eigen value matrix
            for basis in range(nbf):
                matrix[basis,basis] = eVal[basis]

            return matrix

#########################
        def zero(self, matrix):
            #set very small decimal values to 0
            #meant to fix numerical errors in python
            
            #cutoff value
            cutoff = 2.0 * pow(10.0, -15)
        
            #get number of basis functions
            nbf = len(matrix)
        
            #loop through all matrix values,
            #if value less than cutoff, set to zero
            for b1 in range(nbf):
                for b2 in range(nbf):
                    if(abs(matrix[b1,b2]) <= cutoff):
                        matrix[b1,b2] = 0.0

            return matrix        

#########################
        def buildDensity(self, eVec, N):
            #builds density from fock matrix AO eigen vectors and number of electrons
            
            #number of basis functions
            nbf = len(eVec)

            #init zero density matrix 
            D = np.asmatrix(np.zeros([nbf,nbf]))

            for b1 in range(nbf):
                for b2 in range(nbf):
                    for e in range(N/2):
                        D[b1,b2] += eVec[b1,e] * eVec[b2,e]

            return D
#########################
        def energy(self, Hcore, F, D):
            #calculates energy of system

            #get number of basis functions
            nbf = len(F)
            
            #init energy
            E = 0.0

            for b1 in range(nbf):
                for b2 in range(nbf):
                    E +=  D[b1,b2] * (F[b1,b2] + Hcore[b1,b2])

            return E
