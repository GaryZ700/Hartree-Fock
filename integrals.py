#python implementation of Hartree Fock integrals as described in Szabo

import numpy as np

########################
class Integrals:

########################
    def GPT(self, p1, p2):
	#performs guassian product theorem
	#as specified in Szabo pg. 411

	#get required data from primatives
        a1 = p1["A"]
	a2 = p2["A"]
	
	R1 = p1["R"]
	R2 = p2["R"]
	
	n1 = p1["N"]
	n2 = p2["N"]

	#init repeated values
	m = a1 * a2
	a = a1 + a2
	d = m / a
	dist = [ abs( R1[dim] - R2[dim]  ) ** 2.0 for dim in range(3) ]	

	#proportionality constant 
	#Szabo pg. 411 equ A.3
	K = [ np.exp( -d * dist[dim] ) for dim in range(3) ]
	
	#gpt integral constant
	#pg. 412 equ A.9
	constant = pow((np.pi / a), (3.0/2.0))  	

	#calculate primative overlap
	#pg. 412 equ A.9 
	overlap = constant * K[0] * K[1] *K[2] * p1["N"] * p2["N"] * p1["CC"] * p2["CC"]
	
        print("a" + str(a) )
        print("m " + str(m))
        print("D "+str(d) + "\n")
        print("K " + str(K))
	print("Constant "  +  str(constant))
        print("Normalized Overlap " + str(overlap)) 
	
	#create gpt dictionary		
	gpt = {
	
	"m":m,
	"a":a,
	"d":d,
	"dist":dist,
	"K":K,
	"constant":constant,
	"overlap":overlap,
	"C12": p1["CC"] * p2["CC"],
	"N12": p1["N"] * p2["N"]

}
	return gpt

########################
    def overlap(self, basisSet):
    
	#number of basis functions
	nbf = len(basisSet)
	
	#init empty overlap matrix
	S = np.asmatrix(np.zeros([nbf, nbf]))

	print(nbf)

	#loop over basis functions twice
	for b1 in range(nbf):

	#    print("Iterate over b1")

	    for b2 in range(nbf):
	#	print("iterate over b2 " + str(b2))
           		
		#loop over primatives twice
		for p1 in basisSet[b1]:
	#	    print("iterate over p1")
		    print(basisSet[b2])
		    for p2 in basisSet[b2]:
	#		print("iterate over p2")
			
                        print("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                        print("P1")
                        print(p1)
                        #print("\n")
                        print("P2")
                        print(p2)
                        #print("\n")

                        #perform guassian product theorem			
			gpt = self.GPT(p1, p2)
	#		print(gpt["overlap"])	
			#build overlap matrix
			S[b1,b2] += gpt["overlap"]
                        print(S[b1,b2])
			#print("B12 " + str(b1) + " " + str(b2)) 
			
                       # print("\n S:")
                        #print(S[b1,b2])
                        #print("\n")
	

	return S

	
	
