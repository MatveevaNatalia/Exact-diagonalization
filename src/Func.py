from numpy import linalg as LA
from bit_operations import *
from quant_operators import *
from state_operations import *
from quant_operators import *


class HData():
    def __init__(self, basisFock, eigen_values, eigen_vectors, the_matrix, size):
        self.basisFock = basisFock
        self.eigen_values = eigen_values
        self.eigen_vectors = eigen_vectors
        self.the_matrix = the_matrix
        self.size = size


class HCalc():
    def __init__(self, num_level, delta, scat_ampl, num_part):
        self.num_level = num_level
        self.delta = delta
        self.scat_ampl = scat_ampl
        self.num_part = num_part

    def kin_energy(self):
        # Returns the array which contains single-particle energies
        # energy = [0]*params.Num_level
        # energy = [i * params.Delta for i, __ in enumerate(energy)]
        energy = [i * self.delta for i in range(self.num_level)]
        return energy

    def find_scat_elem(self, restr=True):
        # Finds all nonzero matrix elements of scattering operator.
        # If restr = True, assumes the condition:
        # V_{ijkl} != 0 only if i+j = k+l
        # If restr = False, doesn't assume this condition
        # Forming the array which contains unique pairs of indexes:
        # pairs = [[0,1],[0,2],.., [1, 2], [1,3],...,[3,4]]
        pairs = []
        for i in range(0, self.num_level):
            for j in range(i + 1, self.num_level):
                pairs.append([i, j])
        print("Pairs array:")
        print(pairs)

        list_main = []
        for p, [i, j] in enumerate(pairs):
            # can enumerate be used in the second loop also ?
            for q in range(p + 1, len(pairs)):
                [k, l] = pairs[q]
                # i!=j and k!=l are already accounted in 'pairs'
                cond = (i != k and i != l and j != k and j != l)
                if restr:
                    cond = (cond and i + j == k + l)
                if cond:
                    list_main.append([[i, j, k, l], self.scat_ampl])
                    if restr:
                        if (not verification_find_scat_elem([i, j, k, l])):
                            raise RuntimeError("More then one unique index in scattering element!")

        return list_main

    def get_Fock_basis(self, energy_arr, verbose=True):
        v_ground = ground_state(self.num_level, self.num_part)
        state_list = list(perm_unique(v_ground))
        final_list = final_list_create(energy_arr, state_list)
        basisFock = sorted(final_list, key=get_key)

        if verbose:
            print("Ground state vector:")
            print(v_ground)
            print("List of all Fock states:")
            print(state_list)
            print("Number of all Fock states")
            print(len(state_list))
            print("Fock states states written as integers with corresponding kinetic energy:")
            print(final_list)

        return basisFock

    def fill_matrix(self, scat_elem, basis_Fock):
        # Fills the matrix for Hamiltonian. The electrostatic part is not included.
        num_state = len(basis_Fock)
        the_matrix = [[0.0 for x in range(num_state)] for y in range(num_state)]

        for [scat_index, scat_ampl] in scat_elem:
            for i_stat, [__, state_initial] in enumerate(basis_Fock):

                [i, j, k, l] = scat_index

                # Action of V_ijkl on basis state 'state_initial'
                state_temp, coeff = fill_matrix_helper(state_initial, scat_index, self.num_level)

                if coeff != 0:
                    i_col = find_state_index(basis_Fock, state_temp)
                    the_matrix[i_stat][i_col] += coeff * scat_ampl

                    # Action of V_klij on basis state 'state_initial'

                state_temp1, coeff1 = fill_matrix_helper(state_initial, [k, l, i, j], self.num_level)

                if coeff1 != 0:
                    i_col1 = find_state_index(basis_Fock, state_temp1)
                    the_matrix[i_stat][i_col1] += coeff1 * scat_ampl

        for i in range(num_state):
            for j in range(num_state):
                if i == j:
                    the_matrix[i][j] = basis_Fock[i][0]

        return the_matrix

    def h_solver(self, energy_arr, scat_elem, verb_basis=True):
        # This function returns the information about
        #  N-subparticle space: the list "final_list_sorted"
        # contains the basis vectors in binary representation
        # and corresponding energies, list "eigen_vectors"
        # contains eigen vectors.

        basis_Fock = self.get_Fock_basis(energy_arr, verb_basis)
        the_matrix = self.fill_matrix(scat_elem, basis_Fock)
        eigen_values, eigen_vectors = LA.eig(the_matrix)
        size = len(eigen_vectors)

        data = HData(basis_Fock, eigen_values, eigen_vectors, the_matrix, size)

        return data




