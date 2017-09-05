#test integral values from Crawdad Programming: http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3
#used to test that main SCF loop is functioning properly

import numpy as np
import os

######################
class testSCF:

######################
    def testerParse(self):
       
       os.chdir("./tests")
        
       #open test data files
       #build test matricies
       
       s = open("S.dat")
       vext = open("Vext.dat")
       ke = open("KE.dat")
       
       S = np.asmatrix(np.zeros([7,7]))
       Vext = np.asmatrix(np.zeros([7,7])) 
       KE = np.asmatrix(np.zeros([7,7]))
######################
       #adds to respective matrix given file line
       def matrixBuild(matrix, line):
            
           #split line and clean it 
           line = line.split("\n")[0].split(" ")
           
           #parse value and index from line
           value = float(line[len(line)-1])
           index = [float(line[4])-1,float(line[9])-1]
       
       
           #add value to correct position on matrix
           matrix[index[0],index[1]] = value
           matrix[index[1],index[0]] = value
       
           return matrix
       
######################
       #loop through all files and lines in each file
       for line in s:
           S = matrixBuild(S, line)
       
       for line in vext:
           Vext = matrixBuild(Vext, line)
       
       for line in ke:
           KE = matrixBuild(KE, line)
       
       return S, Vext, KE
       


