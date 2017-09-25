#python implementation of Hartree Fock integrals as described in Szabo

import numpy as np
from scipy import special as advMath

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

        #print("!!!!!!!!!!!!!")
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
	
#        print("a" + str(a) )
#        print("m " + str(m))
#        print("D "+str(d) + "\n")
#        print("K " + str(K))
#	print("Constant "  +  str(constant))
#        print("Normalized Overlap " + str(overlap)) 
#	
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

	#print(nbf)

	#loop over basis functions twice
	for b1 in range(nbf):

	#    print("Iterate over b1")

	    for b2 in range(nbf):
	#	print("iterate over b2 " + str(b2))
           		
		#loop over primatives twice
		for p1 in basisSet[b1]:
	#	    print("iterate over p1")
		    #print(basisSet[b2])
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
                        #print(S[b1,b2])
			#print("B12 " + str(b1) + " " + str(b2)) 
			
                       # print("\n S:")
                        #print(S[b1,b2])
                        #print("\n")
	

	return S

########################
    def KE(self, basisSet):
        #uses equation A.11 on page 427 of Szabo

        #print(basisSet)

        #get number of basis functions
        nbf = len(basisSet)
        
#        print("$$$$$$$$$$$$$$$$$")
#       print( "nbf" + str(nbf))
        
#        print("$$$$$$$$$$$$$$$$$")

        #init emppty KE matrix
        T = np.asmatrix( np.zeros( [nbf, nbf]) )

        #loop over basis set twice
        for b1 in range(nbf):

 #           print("$$$$$$$$$$$$$$$$$")
  #          print("basis 1 = ")
   #         print(b1)
    #        print("$$$$$$$$$$$$$$$$$")
            for b2 in range(nbf):
                

     #           print("$$$$$$$$$$$$$$$$$")
     #           print("basis 2 =")
     #           print(b2)
     #           print("$$$$$$$$$$$$$$$$$")
                #loop over primatives twice
                for p1 in basisSet[b1]:
      #              print("$$$$$$$$$$$$$$$$$")
      #              print(b1)
      #              print("$$$$$$$$$$$$$$$$$")
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

#                        print("@@@@@@@@@@@@@@@@@@@@")
#                        print("overlap")
#                        print(gpt["overlap"])
#                       
#                        print("c12")
#                        print(gpt["C12"])
#
#                        print("p1")
#                        print(p1)
#                        print("p2")
#                        print(p2)
# 
#                        print("TermC \n")
#                        print(termC)
#                        #print("\n")
# 
#                        print("Term0 \n")
#                        print(term0)
#                        #print("\n")
# 
#                        print("Term1 \n")
#                        print(term1)
#                        #print("\n")
#                        
#                        print("Term2 \n")
#                        print(term2)
#                        #print("\n")
#                        
#                        print("Term3 \n")
#                        print(term3)
#                        #print("\n") 
#                         
                        T[b1,b2] += term0 - term1 - term2 - term3
                        

                #print("\n T")
                #print(T)
                #print("\n\n")
	
        return T

########################
    def nucAttract(self, basisSet, Z, R):
        #builds nuclear attraction/external potential matrix
        #using eq. A.33 on pg. 430 of Szabo

        #get number of basis functions
        nbf = len(basisSet)

        #init empty nuclear attraction matrix
        Vext = np.asmatrix( np.zeros( [nbf,nbf]) )

        #loop over basis set twice
        for b1 in range(nbf):
            for b2 in range(nbf):

                #loop over primatives twice
                for p1 in basisSet[b1]:
                    for p2 in basisSet[b2]:
                        
                        #loop over atomic numbers:
                        for atom in range(len(Z)):

                            #begin building nuclear attraction matrix
                        
                            constant = -2.0 *  (np.pi / (p1["A"] + p2["A"]) ) * Z[atom]
                            
                            k = self.GPT(p1,p2)["K"]
                            K = k[0] * k[1] * k [2] * p1["N"] * p2["N"]
                         
                            #init values for Boys function
                            a1 = p1["A"]
                            a2 = p2["A"]
                     
                            R1 = p1["R"]
                            R2 = p2["R"]                            

                            a = a1 + a2
                     
                            mid2 = []
                            for dim in range(3):
                               mid2.append( (a1 * R1[dim]) + (a2 * R2[dim]))
                     
                            print("\n$$$$$$$$$$$$$$$")
                            
                     
                            RPA = np.asarray(R[atom]) - np.asarray(mid2)
                            print("OLD RPA")
                            print(RPA)
                            for dim in range(3):
                                if(RPA[dim] < 0.0):
                                    RPA[dim] += 0.0j
                                RPA[dim] = ( RPA[dim] ** 0.2 ) * a  
                            
                            RPA2 = sum(RPA)
                            
                            print("RPA")
                            print(RPA)
                            print("constant")
                            print(constant)
                            print("K")
                            print(K)
                            print("boys")
                            print(self.Boys(0, RPA2))
                            print("RPA2")
                            print(RPA2)
                            print("\n")

                
                            Vext[b1,b2] += constant * K * self.Boys(0, RPA2) * p1["CC"] * p2["CC"]
        return Vext 
                            
########################
    def Boys(self, n, x):
       #python implementation of the any order Boys function       
    
       print("BOYS")
       print("x")
       print(x)
       
       if(x == 0.0):
           return 1.0 / ( (2.0 * n) + 1.0 )
    
       N = n + 0.5
       C = 2.0 * ( x ** N)
       
       print("C")
       print(C)
       print("gamma")
       print(advMath.gamma(N))
       print("gammainc")
       print(advMath.gammainc(N,x))

       return advMath.gamma( N ) * advMath.gammainc( N, x ) / C

########################
    def elecRepulsion(self, basisSet):
        
        #get number of basis functions
        nbf = len(basisSet)

        #init empty electron repulsion matrix
        G = np.zeros( [nbf, nbf, nbf, nbf] )

        #loop over basis sets and primatives four times
        for b1 in range(nbf):
            for b2 in range(nbf):

                for p1 in basisSet[b1]:
                    for p2 in basisSet[b2]:

                        #perform guassian product theory for first basis iteration
                        gpt1 = self.GPT(p1, p2)

                        #get partial overlap for first basis iteration
                        O1 = gpt1["K"][0] * gpt1["K"][1] * gpt1["K"][2] * gpt1["C12"]

                        for b3 in range(nbf):
                            for b4 in range(nbf):
                                
                                for p3 in basisSet[b3]:
                                    for p4 in basisSet[b4]:
                                        
                                        #perform gpt for second basis iteration
                                        gpt2 = self.GPT(p3, p4)
                                        
                                        #get partial overlap for first basis iteration
                                        O2 = gpt2["K"][0] * gpt2["K"][1] * gpt2["K"][2] * gpt2["C12"]

                                        #get sum of all primative alphas
                                        a3 = gpt1["a"] + gpt2["a"]
                                        
                                        #calculate values needed for boys function
                                        A = ( gpt1["a"] + gpt2["a"]  ) / a3
                                        boy = 0.0
                                        for dim in range(3):
                                            boy += (gpt1["mid"][dim] - gpt2["mid"][dim]) ** 2.0  
                                        
                                        #calculate values needed for integral                                        
                                        C1 = self.Boys(0, A * boy ) * O1 * O2
                                        C2 = (2.0 * pow(np.pi, 2.5) ) / ( gpt1["a"] * gpt2["a"] * np.sqrt(a3) ) 
#                                        print("***************")
#                                        print("C1")
#                                        print(C1)
#                                        print("\n C2")
#                                        print(C2)
#                                        print("\n Val")
#                                        print (C1 * C2 * gpt1["N12"] * gpt2["N12"])
#                                        print("\n")
#

                                        #start building integral matrix
                                        G[b1][b2][b3][b4] += C1 * C2 * gpt1["N12"] * gpt2["N12"]
        
        #print("GGGG")
        ##print(G)
        return G

########################
    def HTExchange(self, D, G):
        #builds Hartree Exchange energy from Density and Electron Repulsion Matrix
        
        #get number of basis functions
        nbf = len(G)

        #init empty exchange energy matrix
        HTEx = np.asmatrix( np.zeros([nbf,nbf]) )

        #loop over basis functions four times
        for b1 in range(nbf):
            for b2 in range(nbf):
                for b3 in range(nbf):
                    for b4 in range(nbf):

                        #init term for HTEx
                        term = (2.0*G[b1][b2][b3][b4]) - G[b1][b3][b2][b4]

                        #start building hartree exchange energy
                        HTEx[b1, b2] += D[b3,b4] * term

        
        print("XXXXXXXXX")
        print(HTEx)
        print("XXXXXXXXXXXXX")
        return HTEx



