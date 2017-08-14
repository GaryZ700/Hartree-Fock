#main Hartree-Fock script

import numpy as np
import math

#########################
#init main variables

#scf energy limit
E_limit = 1.0 * pow(10.0, -6)

#########################
#main code goes here

#init integral operators

#hamiltonian
H = np.zeros()

#start of scf loop

#init energy
E = [-float("inf")]

#set guess fock matrix equal to the Hamiltonian
F = H

#scf loop
while(0 == 0):
	
	#transform Fock matrix from AO to MO basis
 
	
	#diagnolize MO fock matrix
	eVal, eVec = np.linalg.eig(Fp)



