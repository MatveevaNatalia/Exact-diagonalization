import numpy as np
import sys
import os.path
sys.path.insert(0, "./src")
from Func import *
from quant_operators import *


def wf_amp(index_lead, index_state):
    #Returns the wave function amplitude at the point of the contact
    #with electrode l (l = left or right)

    if index_lead == 'left':
        amplitude = 1
    elif index_lead == 'right':
        if index_state % 2 == 1:
            amplitude = 1 # ORTOHOGONALITY IN DANGER!
        else:
            amplitude = -1
    else:
        raise ValueError("Unrecognized index is given to WaveFunctionAmplitude !")

    return amplitude


def dens_state(x, gap):
    if x > gap:
        return 2.*x / np.sqrt(x * x - gap * gap)
    else:
        return 0.


class Rate():
    def __init__(self, gamma, vol, gap, num_lev, mu, wf, H_data_np1, H_data_n):
        self.gamma = gamma
        self.vol = vol
        self.gap = gap
        self.num_lev = num_lev
        self.mu = mu
        self.wf = wf
        self.H_data_np1 = H_data_np1
        self.H_data_n = H_data_n

    def removal_rate(self, index_left, index_right):

        #This function calculates
        #R_{ab} = gamma|\sum_{m=0}^{M-1}wf_m*<\phi^N_a|\hat{c}_m|\psi^{N+1}_b>|^2*nu(x),
        # x = E_b + μ − E_a + V
        #for given N and N+1.

        num_lev = self.num_lev
        basisFock_n = self.H_data_n.basisFock
        basisFock_np1 = self.H_data_np1.basisFock
        eig_vec_np1 = self.H_data_np1.eigen_vectors
        eig_vec_n = self.H_data_n.eigen_vectors
        eig_val_np1 = self.H_data_np1.eigen_values
        eig_val_n = self.H_data_n.eigen_values

        sum_total = 0

        # Loop over indexes for annihilation operator \hat{c}_index
        for index in range(self.num_lev):

            ket_np1 = eigen_vector_ket(basisFock_np1, eig_vec_np1, index_right)
            ket_final = kill_eigen_state(ket_np1, num_lev, index)

            # The complex conjugation in the left bra vector is not considered here,
            # because all elements of eigen states are chosen to be real.
            ket_n = eigen_vector_ket(basisFock_n, eig_vec_n, index_left)

            # Now we can calculate <phi_{index_left}|c_index|phi_{index_right}>
            scal_prod = 0

            for [ff, ee] in ket_n:
                for [ff1, ee1] in ket_final:
                    if ee == ee1:
                        scal_prod += ff * ff1

            sum_total += scal_prod * self.wf[index]

        Ea = eig_val_n[index_left]
        Eb = eig_val_np1[index_right]

        temp = Eb + self.mu - Ea + self.vol

        ds = dens_state(temp, self.gap)

        #return sum_total*sum_total

        return 0.5*self.gamma * sum_total * sum_total * ds

    def addition_rate(self, index_left, index_right):

        # This function calculates
        # A_{ba} = gamma|\sum_{m=0}^{M-1}wf_m*<\phi^N+1_b|\hat^+{c}_m|\phi^N_a>|^2*nu(x),
        # x = E_a − V − E_b − μ
        # for given N and N+1.

        num_lev = self.num_lev
        basisFock_n = self.H_data_n.basisFock
        basisFock_np1 = self.H_data_np1.basisFock
        eig_vec_np1 = self.H_data_np1.eigen_vectors
        eig_vec_n = self.H_data_n.eigen_vectors
        eig_val_np1 = self.H_data_np1.eigen_values
        eig_val_n = self.H_data_n.eigen_values

        sum_total = 0
        # Loop over indexes for annihilation operator \hat{c}_index
        for index in range(self.num_lev):

            ket_np1 = eigen_vector_ket(basisFock_np1, eig_vec_np1, index_left)
            ket_n = eigen_vector_ket(basisFock_n, eig_vec_n, index_right)
            ket_final = create_eigen_state(ket_n, num_lev, index)

            # The complex conjugation in the left bra vector is not considered here,
            # because all elements of eigen states are chosen to be real.

            # Now we can calculate <phi_{index_left}|c_index|phi_{index_right}>
            scal_prod = 0

            for [ff, ee] in ket_np1:
                for [ff1, ee1] in ket_final:
                    if ee == ee1:
                        scal_prod += ff * ff1

            sum_total += scal_prod * self.wf[index]

        Ea = eig_val_n[index_right]
        Eb = eig_val_np1[index_left]

        temp = Ea - self.vol - Eb - self.mu

        ds = dens_state(temp, self.gap)

        #return sum_total*sum_total

        return 0.5 * self.gamma * sum_total * sum_total * ds


class Tunneling():
    def __init__(self, rate_l, rate_r):
        self.rate_l = rate_l
        self.rate_r = rate_r

    def matrix(self):
        # Returns the matrix for the lhs of the ODE system
        # for probabilities {p} and {q}
        size_n = self.rate_l.H_data_n.size
        size_np1 = self.rate_l.H_data_np1.size
        size_tot = size_n + size_np1

        the_matrix = [[0.0 for x in range(size_tot)] for y in range(size_tot)]

        print("Removal rate, left")
        for a in range(size_n):
            for b in range(size_np1):
                r_ab_l = self.rate_l.removal_rate(a, b)
                print("a= ", a, "b=", b, r_ab_l)

        print("Addition rate, left")
        for a in range(size_n):
            for b in range(size_np1):
                a_ba_l = self.rate_l.addition_rate(b, a)
                print("b= ", b, "a=", a, a_ba_l)

        print("Removal rate, right")
        for a in range(size_n):
            for b in range(size_np1):
                r_ab_r = self.rate_r.removal_rate(a, b)
                print("a= ", a, "b=", b, r_ab_r)

        print("Addition rate, right")
        for a in range(size_n):
            for b in range(size_np1):
                a_ba_r = self.rate_r.addition_rate(b, a)
                print("b= ", b, "a=", a, a_ba_r)

        # Filling the rows that correspond to p_alpha
        s1 = 0
        for i_row in range(size_n):

            s1 = self.sum_helper(i_row, "add")

            the_matrix[i_row][i_row] = s1

            for i_col in range(size_n, size_tot):
                r_ab_l = self.rate_l.removal_rate(i_row, i_col - size_n)
                r_ab_r = self.rate_r.removal_rate(i_row, i_col - size_n)

                the_matrix[i_row][i_col] = r_ab_l + r_ab_r

        # Filling the rows that correspond to q_beta

        s2 = 0
        for i_row in range(size_n, size_tot):

            s2 = self.sum_helper(i_row - size_n, "remove")

            the_matrix[i_row][i_row] = s2

            for i_col in range(size_n):
                a_ab_l = self.rate_l.addition_rate(i_row - size_n, i_col)
                a_ab_r = self.rate_r.addition_rate(i_row - size_n, i_col)

                the_matrix[i_row][i_col] = a_ab_l + a_ab_r

        return the_matrix

    def sum_helper(self, index, key):
        # For key = "add" returns the sum_b(A_{b, index}),
        # for key = "remove" returns the sum_{a}(R_{a, index}).
        # Here A is the addition rate, and
        # R is the removal rate.

        s = 0.
        if key == "add":
            size_np1 = self.rate_l.H_data_np1.size
            for b in range(size_np1):
                a_ba_l = self.rate_l.addition_rate(b, index)
                a_ba_r = self.rate_r.addition_rate(b, index)

                s -= a_ba_l + a_ba_r

        elif key == "remove":
            size_n = self.rate_l.H_data_n.size
            for a in range(size_n):
                r_ab_l = self.rate_l.removal_rate(a, index)
                r_ab_r = self.rate_r.removal_rate(a, index)
                s -= r_ab_l + r_ab_r
        else:
            raise ValueError("Wrong key is given to sum_helper !")
        return s

    def der_func(self, y, t, the_matrix):
        # This function returns derivatives for odeint
        size_n = self.rate_l.H_data_n.size
        size_np1 = self.rate_l.H_data_np1.size

        # y[0]= dp0/dt = list_der[0], ...,
        # y[size_n -1]= dp_{size_n-1} = list_der[size_n-1]
        # y[size_n] = dq0_dt = list_der[size_n], ...,
        # y[size_n+size_np-1] = dq_{size_np -1} = list_der[size_n + size_np-1]

        size_tot = size_n + size_np1
        list_der = [0] * size_tot

        for i_row in range(size_tot):
            for i_col in range(size_tot):
                list_der[i_row] += the_matrix[i_row][i_col] * y[i_col]

        return list_der

