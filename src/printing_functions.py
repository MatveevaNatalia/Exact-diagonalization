import numpy as np


def print_eigen_values(eigen_values):
    print("Eigen values: ")
    print(eigen_values)


def print_eigen_vectors(eigen_vectors):
    print("Eigen vectors: ")
    print(eigen_vectors)


def print_kin_energy(energy_arr):
    print("Array of single-particle energies:")
    for i in range(len(energy_arr)):
        print("Level index: ", i, " single particle energy: ", energy_arr[i])


def print_scat_elem(scat_elem):
    print("List of non-zero UNIQUE scattering elements")
    print("(in the format [[i,j,k,l],V_ijkl]):")
    for i in range(len(scat_elem)):
        print(i, scat_elem[i])

def print_Fock_basis(basis_Fock):
    print("Fock basis")
    print(basis_Fock)
    #print("Printing a given eigen vector:")
    #print(eigen_system_Nplus1.eigen_vectors[:, 1])

def print_matrix(mat):
    # Prints the matrix for Hamiltonian
    #np.set_printoptions(precision=3)
    #print (" ","   ".join([str(x) for x in range(len(mat))]))
    print("  ", "   ".join([str((round(x,3))) for x in range(len(mat))]))
    for i,x in enumerate(mat):
        print (i," ".join([str(round(y,3)) for y in x]))