import numpy as np
from numpy import linalg as LA
import sys
import os.path
sys.path.insert(0, "./src")
from Func import *

params_Nplus1 = ParamModel(Num_level = 5, Delta = 1.0, Num_part = 3, Scat_ampl = 1.0)



print("#####################################################################################")
print("N+1 PARTICLE SUBSPACE")
print("#####################################################################################")
# Setting the parameters of the model
params_Nplus1 = ParamModel(Num_level=5, Delta=1.0, Num_part=3, Scat_ampl=1.0)

# Formation of an array 'energy_arr' which contains the single-particle energy spectrum

energy_arr = Kin_Energy(params_Nplus1)
print("Array of single-particle energies:")
for i in range(len(energy_arr)):
    print("Level index: ", i, " single particle energy: ", energy_arr[i])

# energy_arr[3] = 10. <- way to control the energies manually

# Formation of an array 'scat_elem' which contains scattering elements


scat_elem = Find_Scat_Elem(params_Nplus1, restr=True)

print(scat_elem)
print(scat_elem[0][1])

print("List of non-zero UNIQUE scattering elements")
print("(in the format [[i,j,k,l],V_ijkl]):")
for i in range(len(scat_elem)):
    print(i, scat_elem[i])

# Changing a scattering element

index = 0  # number of scat. elem. in 'scat_elem' list
value = 10.  # new value of scat. elem.

Set_Scat_Element(scat_elem, index, value)



print("List of non-zero UNIQUE scattering elements after changing V_ijkl:")

# Changing the same scattering element back

index = 0
value = 1.

for i in range(len(scat_elem)):
    print(i, scat_elem[i])

Set_Scat_Element(scat_elem, index, value)

print("List of non-zero UNIQUE scattering elements after changing V_ijkl:")

for i in range(len(scat_elem)):
    print(i, scat_elem[i])

for i in range(len(scat_elem)):
    print(i, scat_elem[i])

basis_Nplus1, eigen_vectors_Nplus1, myMatrix = SubspaceInfo(params_Nplus1, energy_arr, scat_elem, True, True)

print("Printing a given eigen vector:")
print(eigen_vectors_Nplus1[:, 1])

print("Basis vectors with corresponding energy:")
print(basis_Nplus1)
print("Number of basis vectors:")
print(len(basis_Nplus1))
print("Eigen vectors of Hamiltonian matrix (columns):")
print(eigen_vectors_Nplus1)


print("#####################################################################################")
print("N PARTICLE SUBSPACE")
print("#####################################################################################")

params_N = ParamModel(Num_level=4, Delta=1.0, Num_part=1, Scat_ampl=1.0)

# Formation of an array 'energy_arr' which contains the single-particle energy spectrum
energy_arr = Kin_Energy(params_N)
print("Array of single-particle energies:")
for i in range(len(energy_arr)):
    print("Level index: ", i, " single particle energy: ", energy_arr[i])
# energy_arr[3] = 10. <- way to control the energies manually

# Formation of an array 'scat_elem' which contains scattering elements


scat_elem = Find_Scat_Elem(params_N, restr=True)

print("List of non-zero scattering elements")
print("(in the format [[i,j,k,l],V_ijkl]):")
for i in range(len(scat_elem)):
    print(i, scat_elem[i])

basis_N, eigen_vectors_N, my_Matrix = SubspaceInfo(params_N, energy_arr, scat_elem, True, True)

print("Basis vectors with corresponding energy:")
print(basis_N)
print("Number of basis vectors:")
print(len(basis_N))
print("Eigen vectors of Hamiltonian matrix:")
print(eigen_vectors_N)


print("#####################################################################################")
print("MATRIX ELEMENTS FOR MASTER EQUATION")
print("#####################################################################################")
MastEqPrep(basis_Nplus1, eigen_vectors_Nplus1, params_Nplus1, basis_N, eigen_vectors_N, False)