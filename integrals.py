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

        print("!!!!!!!!!!!!!")
	a = a1 + a2
        print(a)
	d = m / a
	dist = [ abs( R1[dim] - R2[dim]  ) ** 2.0 for dim in range(3) ]	
        
        mid = []
        for dim in range(3):
            mid.append( ( (a1 * R1[dim]) + (a2 * R2[dim]) / a ) ) 

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
        "mid":mid,
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
			
                       # print("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                       # print("P1")
                       # print(p1)
                        #print("\n")
                       # print("P2")
                        #print(p2)
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

########################
    def KE(self, basisSet):
        #uses equation A.11 on page 427 of Szabo

        print(basisSet)

        #get number of basis functions
        nbf = len(basisSet)
        
        print("$$$$$$$$$$$$$$$$$")
        print( "nbf" + str(nbf))
        
        print("$$$$$$$$$$$$$$$$$")

        #init emppty KE matrix
        T = np.asmatrix( np.zeros( [nbf, nbf]) )

        #loop over basis set twice
        for b1 in range(45):

            print("$$$$$$$$$$$$$$$$$")
            print("basis 1 = ")
            print(b1)
            print("$$$$$$$$$$$$$$$$$")
            for b2 in range(2):
                

                print("$$$$$$$$$$$$$$$$$")
                print("basis 2 =")
                print(b2)
                print("$$$$$$$$$$$$$$$$$")
                #loop over primatives twice
                for p1 in basisSet[b1]:
                    print("$$$$$$$$$$$$$$$$$")
                    print(b1)
                    print("$$$$$$$$$$$$$$$$$")
                    for p2 in basisSet[b2]:

                        #begin calculating the KE matrix
                        
                        gpt = self.GPT(p1, p2)

                        
                        constant = 1.0 / ( 2.0 * gpt["a"] )
                        
                        termC = []
                        for dim in range(3):
                            termC.append( (gpt["mid"][dim] - p2["R"][dim]) * (gpt["mid"][dim] - p2["R"][dim]) + constant )
                        

                        term0 = 3.0 * p2["A"] * gpt["overlap"]
                        term1 = 2.0 * p2["A"] * p2["A"] * termC[0] * gpt["overlap"] 
                        term2 = 2.0 * p2["A"] * p2["A"] * termC[1] * gpt["overlap"] 
                        term3 = 2.0 * p2["A"] * p2["A"] * termC[2] * gpt["overlap"] 

                        print("@@@@@@@@@@@@@@@@@@@@")
                        print("overlap")
                        print(gpt["overlap"])
                        
                        print("c12")
                        print(gpt["C12"])

                        print("p1")
                        print(p1)
                        print("p2")
                        print(p2)

                        print("TermC \n")
                        print(termC)
                        #print("\n")

                        print("Term0 \n")
                        print(term0)
                        #print("\n")

                        print("Term1 \n")
                        print(term1)
                        #print("\n")
                        
                        print("Term2 \n")
                        print(term2)
                        #print("\n")
                        
                        print("Term3 \n")
                        print(term3)
                        #print("\n") 
                        
                        T[b1,b2] += term0 - term1 - term2 - term3
                        

                print("\n T")
                print(T)
                print("\n\n")
	
        return T

########################
