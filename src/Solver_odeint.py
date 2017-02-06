import numpy as np
import sys
import os.path
sys.path.insert(0, "./src")
from scipy.integrate import odeint
from Rate_equation import *
from Plotting import *


class SolverRateEquation():
    def __init__(self, my_tunnel_1, my_tunnel_2, param_particle, param_rate):
        # Instance of the Tunneling class (Rate_equation.py) for the 1st half of period
        self.my_tunnel_1 = my_tunnel_1
        # Instance of the Tunneling class (Rate_equation.py) for the 2nd half of period
        self.my_tunnel_2 = my_tunnel_2
        # Dictionary containing single particle parameters
        self.param_particle = param_particle
        # Dictionary containing parameters for the rate equations
        self.param_rate = param_rate
        # Dimension of the Fock space for N particles
        self.size_n = self.my_tunnel_1.rate_l.H_data_n.size
        # Dimension of the Fock space for N+1 particles
        self.size_np1 = self.my_tunnel_1.rate_l.H_data_np1.size
        self.size_tot = self.size_n + self.size_np1


    def calc_odeint(self, init_vec_1, i):
        # Solves the system of rate equations at given period using standard 'odeint()' solver
        # 'init_vec_1' are the initial conditions at the begining of a given period
        # 'i' indicates the given period.

        num_T = self.param_rate['num_T']
        period = self.param_rate['T']
        num_points = self.param_rate['num_points']
        init_vec_2 = np.zeros(self.size_tot)

        t = []
        y = []

        the_matrix_1 = self.my_tunnel_1.matrix()
        the_matrix_2 = self.my_tunnel_2.matrix()

        t1 = np.linspace(i * period, i * period + period / 2., num_points)
        y1 = odeint(self.my_tunnel_1.der_func, init_vec_1, t1, args=(the_matrix_1,))

        init_vec_2 = y1[num_points - 1, :]

        t2 = np.linspace(i * period + period / 2., (i + 1) * period, num_points)

        y2 = odeint(self.my_tunnel_2.der_func, init_vec_2, t2, args=(the_matrix_2,))

        return t1, y1, t2, y2

    def calc_diag(self, init_vec_1, i):
        # Solves the system of rate equations at given period using
        # diagonalization method (the method itself contains in 'calc_diag_solver()' method).
        # 'init_vec' are the initial conditions at the begining of a given period
        # 'i' indicates the given period.

        period = self.param_rate['T']
        num_points = self.param_rate['num_points']

        the_matrix_1 = self.my_tunnel_1.matrix()
        the_matrix_2 = self.my_tunnel_2.matrix()

        t1 = np.linspace(0, period / 2., num_points)
        y1 = self.calc_diag_solver(the_matrix_1, init_vec_1, t1)

        init_vec_2 = y1[num_points - 1, :]

        t2 = np.linspace(0, period/2., num_points)

        y2 = self.calc_diag_solver(the_matrix_2, init_vec_2, t2)

        t1 = np.linspace(i * period, i * period + period / 2., num_points)
        t2 = np.linspace(i * period + period / 2., (i + 1) * period, num_points)

        return t1, y1, t2, y2, init_vec_2

    def calc_diag_solver(self, the_matrix, init_vec, t):
        # Solver of a system of linear differential equations
        # using diagonalization method.
        # 'the matrix' is the matrix for the rhs of the system,
        # 'init_vec' is the array containing the initial conditions,
        # 't' is the array containing time at which the solution needs to be found.

        lam, evec = LA.eig(the_matrix)
        coef = LA.solve(evec, init_vec)
        sol = []

        for i in t:
            sol.append(np.dot(evec, coef * np.exp(lam * i)))

        sol1 = np.array(sol).real
        return sol1


    def get_labels(self):
        # Prepares the labels for the 'populations_plot()' function (Plotting.py)
        q = [0] * self.size_np1
        p = [0] * self.size_n

        num_level = self.param_particle['num_level']
        for i in range(self.size_n):
            p[i] = str(to_bitfield(self.my_tunnel_1.rate_l.H_data_n.basisFock[i][1], num_level))

        for i in range(self.size_np1):
            q[i] = str(to_bitfield(self.my_tunnel_1.rate_l.H_data_np1.basisFock[i][1], num_level))

        return p, q

    def plotting(self, t1, y1, t2, y2):
        # Calls 'populations_plot()' with the proper labels.
        # Receives arrays containg the populations of the Fock levels for the first and the
        # second half of period.
        # Also receives arrays 't1' and 't2' that contains the corresponding time.
        p, q = self.get_labels()
        populations_plot(t1, t2, y1, y2, p, q, self.param_rate)

    def single_level_calc(self, y1, y2):
        # Calculates the populations of single particle levels.
        # Receives the arrays 'y1' and 'y2' that contains the populations
        # of Fock levels for the first and second half of the period, correspondingly.
        num_level = self.param_particle['num_level']
        num_points = self.param_rate['num_points']
        single_level1 = np.zeros((num_points,num_level))
        single_level2 = np.zeros((num_points,num_level))

        for i_sl in range(num_level):
            for i_fl in range(self.size_n):
                state = self.my_tunnel_1.rate_l.H_data_n.basisFock[i_fl][1]
                bit = to_bitfield(state, num_level)[i_sl]
                single_level1[:, i_sl] += y1[:, i_fl] * bit
                #print(i_sl, i_fl, y1[0, i_fl], bit)
                single_level2[:, i_sl] += y2[:, i_fl] * bit

            for i_fl in range(self.size_np1):
                state = self.my_tunnel_1.rate_l.H_data_np1.basisFock[i_fl][1]
                bit = to_bitfield(state, num_level)[i_sl]
                single_level1[:, i_sl] += y1[:, self.size_n + i_fl] * bit
                single_level2[:, i_sl] += y2[:, self.size_n + i_fl] * bit

        return single_level1, single_level2

