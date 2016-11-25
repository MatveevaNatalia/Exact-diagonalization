import numpy as np
import sys
import os.path
sys.path.insert(0, "./src")
from scipy.integrate import odeint
from Rate_equation import *
from Plotting import *


class SolverRateEquation():
    def __init__(self, my_tunnel_1, my_tunnel_2, param_particle, param_rate):
        self.my_tunnel_1 = my_tunnel_1
        self.my_tunnel_2 = my_tunnel_2
        self.param_particle = param_particle
        self.param_rate = param_rate
        self.size_n = self.my_tunnel_1.rate_l.H_data_n.size
        self.size_np1 = self.my_tunnel_1.rate_l.H_data_np1.size
        self.size_tot = self.size_n + self.size_np1

    def calc_odeint(self, init_vec_1):
        num_T = self.param_rate['num_T']
        T = self.param_rate['T']
        num_points = self.param_rate['num_points']
        init_vec_2 = [0] * self.size_tot

        t = []
        y = []

        the_matrix_1 = self.my_tunnel_1.matrix()
        the_matrix_2 = self.my_tunnel_2.matrix()

        for i in range(num_T):

            if i != 0:
                for i in range(self.size_tot):
                    init_vec_1[i] = y2[num_points - 1, i]

            t1 = np.linspace(i * T, i * T + T / 2., num_points)
            y1 = odeint(self.my_tunnel_1.der_func, init_vec_1, t1, args=(the_matrix_1,))

            for j in range(self.size_tot):
                init_vec_2[j] = y1[num_points - 1, j]

            t2 = np.linspace(i * T + T / 2., (i + 1) * T, num_points)
            #print(i * T + T / 2., (i + 1) * T, t2[num_points - 1])

            y2 = odeint(self.my_tunnel_2.der_func, init_vec_2, t2, args=(the_matrix_2,))

        return t1, y1, t2, y2

    def get_labels(self):
        q = [0] * self.size_np1
        p = [0] * self.size_n

        num_level = self.param_particle['num_level']
        for i in range(self.size_n):
            p[i] = str(to_bitfield(self.my_tunnel_1.rate_l.H_data_n.basisFock[i][1], num_level))

        for i in range(self.size_np1):
            q[i] = str(to_bitfield(self.my_tunnel_1.rate_l.H_data_np1.basisFock[i][1], num_level))

        return p, q

    def plotting(self, t1, y1, t2, y2):
        p, q = self.get_labels()
        populations_plot(t1, t2, y1, y2, p, q, self.param_rate)

