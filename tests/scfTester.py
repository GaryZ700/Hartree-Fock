#tests SCF part of HartreeFock for any errors

import numpy as np
from scf import SCF
scf = SCF()

########################
#test getTransform using equation 3.170 on page 159

#create mock overlap matrix
S = np.random.rand([5,5])

print("S:")
print(S)
print("\n")

#create transform matrix
X = scf.getTransform(S)

print("X:")
print(X)
print("\n")

#run test case
test0 = X.conjugate().transpose() * S * X

print(test0)
print("\n")

