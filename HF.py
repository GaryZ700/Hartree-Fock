#main Hartree-Fock script
#numbers refer to steps listed in Szabo page 161

import numpy as np
import math
from scipy.sparse.linalg import eigs as eig

from scf import SCF
from basis import Basis
from tests.testValues import testSCF
from integrals import Integrals

scf = SCF() 
basis = Basis()
testValues = testSCF()
integrals = Integrals()

#########################
#init main variables

#1
#system contains information for physical atoms/electrons
#uses notation est. in Szabo page 161, step 1
system = {
        
        #atomic coordinates
        "R":[[0,0,0.5], [3.0,2.0,5.0], [0.0,4.0,6.0]],
        #atomic numbers
        "Z":[1,1,1],
        #number of electrons
        "N":10
        
        }

#build basis set for system
basisSet = basis.buildBasis(system)

print(basisSet)

#scf energy limit
E_limit = 1.0 * pow(10.0, -6)

#########################
#main code goes here

#init integral operators
print("\n \n")
print(integrals.overlap(basisSet))
#print("###########")
#init test operators to check program is working
S, Vext, T = testValues.testerParse()

#hamiltonian is equal to kinetic energy plus external potential
#pg.176 Equ. 3.233
Hcore = T + Vext

#########################
#start of scf loop

#get transformation matrix X
X = scf.getTransform(S)
X = scf.zero(X)

#print("\n X:")
#print(X)
#print("\n")

#init energy
E = [-float("inf")]

#4
#set guess fock matrix equal to the Hamiltonian
F = Hcore

#scf loop
while(0 == 0):
	
        #7
	#transform Fock matrix from AO to MO basis
        FMO = X.conjugate().transpose() * F * X
        FMO = scf.zero(FMO)

#        print("FMO")
#        print(FMO)
#        print("\n")

        #8
        #get eigen values and vectors of FMO
        #and get eigen value matrix
        eVal, eVec = np.linalg.eigh(FMO)
    
        print(eVec)
        print(eVal)

        #9
        #transform eigen vectors to AO basis
        eVecAO = np.dot(X, eVec)
        
#        print("eVec AO")
#        print(eVecAO)
#        print("\n")
            
        #calculate density matrix
        D = scf.buildDensity(eVecAO, system["N"])
        
#        print("Density")
#        print(D)
#        print("\n")
    
        #calculate electronic energy
        E = scf.energy(Hcore, F, D)
        
        print(E)

        break
