#basis system used for Python implementation of Hartree-Fock method

import math

#########################
class Basis:

#########################
    def buildBasis(self, system):
        #builds basis dictionary for main program
        #uses notation as established in page 161, step 1 of Szabo

        basisSet = []

	nbf = -1	

	#loop through all atoms in system
	for atom in range(len(system["Z"])):
	    
	    #load all coeff data needed for atom basis
	    alphas, contractions  = self.loadBasis(system["Z"][atom])
	    center = system["R"][atom]	    

	    #for each basis used to represent one atom
            for basis in range(len(alphas)):
	        
                basisSet.append([])
		nbf += 1
	    
                for primative in range(len(alphas[basis])):
    	
      	            primativeData = {
    
    	               "CC" : contractions[basis][primative],
    		       "A" : alphas[basis][primative],
    		       "N" : (2.0 * alphas[basis][primative] / math.pi) ** 0.75,
    		       "R" : center 
                    }
    		
                    basisSet[nbf].append(primativeData)

        return basisSet

#########################
    def loadBasis(self, atom):
	#load basis given atom number
	#returns two lists
	#first of alpha values for the guassian
	#second of contraction coeffs for the guassians
	#all basis set from the basis set exchange
        	
	#hydrogen atom
	if(atom == 1):

	    #return [[3.42525091,  0.62391373, 0.16885540]], [[0.15432897,0.53532814 ,0.44463454]]
            return [ [ 18.7311370, 2.8253937, 0.6401217 ], [ 0.1612778 ] ], [[ 0.0334945, 0.2347269,  0.81375733], [ 1.0 ] ]

        #if(atom == 2):
           # return [ [6.36242139, 1.15892300, 0.31364979 ], [  ]

