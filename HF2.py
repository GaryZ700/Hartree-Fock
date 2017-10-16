#2nd iteration of Hartree Fock Python Implementation, with use of Pyquante as an integral library
#Each numbered section refers to coressponding step in Szabo QM Textbook on page 161

from PyQuante import Molecule, Ints
from PyQuante import LA2 as linalg
from PyQuante.NumWrap import eigh, matrixmultiply
from PyQuante import hartree_fock as HF

#Global Variables############################
convergenceLimit = 1.0 * pow(10, -6)
maxCycle = 50

#Section 1############################
#specify a molecule
molecule = Molecule(
        
        "H2", 
        [
            (1, (0,0,0)),
            (1, (0,0,1)),
            (8, (-1,0,0))
            
            ]
        
        )

basisSet = Ints.getbasis(molecule, "sto-3g")

#Section 2############################
#Overlap Matrix
S = Ints.getS(basisSet)

#Follwing Two matrices compose the core Hamiltonian
#KE Matrix
KE = Ints.getT(basisSet)

#External Potential, Nuclear - Electron Attraction
Vext = Ints.getV(basisSet, molecule)

#Form Hcore
Hcore = KE + Vext

#calculate two electron integrals
elecRepulsion = Ints.get2ints(basisSet)

#Section 3############################
#Calculate Transformation Matrix X from Overlap matrix
X = linalg.SymOrth(S)

#Section 4############################
#guess density matrix from Core Hamiltonian

#diagnolize Hamiltonian and transform hamiltonian to get guess orbitals
eVal, eVec = eigh(linalg.simx(Hcore, X))

#transform orbitals from AO to MO basis
orbitals = matrixmultiply(Hcore, X)

#guess Density from guess orbitals and number of closed electron shells
D = HF.mkdens(orbitals, 0, molecule.get_closedopen()[0])

#Start of SCF Procedures############################
#init variables for SCF
cycle = 0
energy = [float("-inf")]

while(cycle < maxCycle):
    #while max. number of cycles not reached
    #continues with SCF iterations

#Start of SCF Procedures############################

    #Section 5############################
    #calculate hartree exchange energy, G
    #from two electron integrals, and density matrix
    G = HF.get2JmK(elecRepulsion, D)

    #Section 6############################
    #create Fock matrix from hartree exchange and core hamiltonian
    F = Hcore + G

    #Section 7############################
    #Transform Fock matrix to MO basis
    Fmo = linalg.simx(F, X)

    #Section 8############################
    #get eigen values and vectors of transformed Fock matrix
    eVal, eVec = eigh(Fmo)

    #Section 9############################
    #transform eigen vectors back to AO basis to get C
    C = matrixmultiply(X, eVec)

    #Section 10############################
    #calculate a new density matrix from C
    D = HF.mkdens(C, 0, molecule.get_closedopen()[0])

    #Section 11############################
    #check convergence critera, 
    #which for this particular program is checking 
    #the delta energy of the system
    E = HF.get_energy(Hcore, F, D, molecule.get_enuke())

    #append to energy list
    energy.append(E)

    #increment cycle counter
    cycle += 1

    #if the change in energy from one iteration to the next
    #is less than the specfied convergence difference, then end SCF procedure
    if(abs(energy[cycle] - energy[cycle-1]) <= convergenceLimit):
        break


#End of SCF procedure,
#print important information
print("Emergy: " + str(energy[cycle]) + " Hartrees")

