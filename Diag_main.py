import numpy as np
import time
from numpy import linalg as LA
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys
import os.path
sys.path.insert(0, "./src")
from Func import *
from Rate_equation import * # To write what is imported
from printing_functions import *
from Plotting import populations_plot
from Solver_odeint import *
from Current import *

num_level = 4
delta = 1.0 # single-particle level spacing
scat_ampl = 0.0
num_part = 2

param_particle = {
    'num_level': num_level,
    'delta': delta,
    'scat_ampl': scat_ampl,
    'num_part': num_part
}


print("#####################################################################")
print("N+1 PARTICLE SUBSPACE")
print("#####################################################################")


H_Np1 = HCalc(num_level, delta, scat_ampl, num_part)

#energy_arr_Np1 = H_Np1.kin_energy()

energy_arr_Np1 = [-1.5*delta, -0.5*delta, 0.5*delta, 1.5*delta]

print_kin_energy(energy_arr_Np1)

scat_elem_Np1 = H_Np1.find_scat_elem(restr=True)
print_scat_elem(scat_elem_Np1)

#Set_Scat_Element(scat_elem, index, value)

H_data_Np1 = H_Np1.h_solver(energy_arr_Np1, scat_elem_Np1, verb_basis=True)

print_Fock_basis(H_data_Np1.basisFock)
print_matrix(H_data_Np1.the_matrix)
print_eigen_values(H_data_Np1.eigen_values)
print_eigen_vectors(H_data_Np1.eigen_vectors)



print("######################################################################")
print("N PARTICLE SUBSPACE")
print("######################################################################")

H_N = HCalc(num_level, delta, scat_ampl, num_part-1)

#energy_arr_N = H_N.kin_energy()

energy_arr_N = [-1.5*delta, -0.5*delta, 0.5*delta, 1.5*delta]


print_kin_energy(energy_arr_N)

scat_elem_N = H_N.find_scat_elem(restr=True)
print_scat_elem(scat_elem_N)

# possibility to manyally set scattering elements:
#Set_Scat_Element(scat_elem, index, value)

H_data_N = H_N.h_solver(energy_arr_N, scat_elem_N, verb_basis=True)

print_Fock_basis(H_data_N.basisFock)
print_matrix(H_data_N.the_matrix)
print_eigen_values(H_data_N.eigen_values)
print_eigen_vectors(H_data_N.eigen_vectors)

print("######################################################################")
print("RATE EQUATION PREPARATION")
print("######################################################################")

# Setting the model parameters

gamma_L = 100. # escape rate from the normal island into left electrode
gamma_R = 100. # escape rate from the normal island into right electrode
V_L = -2. # voltage on the left electrode
V_R = 2. # voltage on the right electrode

gap = 3. # superconducting gap

mu_minus = -2.
mu_plus = 2.

T = 1.
num_points = 1000
num_T = 20

# Tolerance for finding the zero values for eigen values of rate equations matrix
tolerance = 1.e-16

i_current = True


param_rate = {
    'gamma_L': gamma_L,
    'gamma_R': gamma_R,
    'V_L': V_L,
    'V_R': V_R,
    'gap': gap,
    'mu_minus': mu_minus,
    'mu_plus': mu_plus,
    'scat_ampl': scat_ampl,
    'T': T,
    'num_points': num_points,
    'num_T': num_T,
    'tolerance': tolerance
}


wf_L = [0.5, 0.5, 0.5, 0.5]
wf_R = [0.5, -0.5, 0.5, -0.5]



print("wf_L= ", wf_L)
print("wf_R= ", wf_R)

# Constructing the matrix for first half of period

rate_L_1 = Rate(gamma_L, V_L, gap, num_level, mu_minus, wf_L, H_data_Np1, H_data_N)
rate_R_1 = Rate(gamma_R, V_R, gap, num_level, mu_minus, wf_R, H_data_Np1, H_data_N)

my_tunnel_1 = Tunneling(rate_L_1, rate_R_1)


the_matrix_1 = my_tunnel_1.matrix()
print_matrix(the_matrix_1)

# Constructing the matrix for second half of period

rate_L_2 = Rate(gamma_L, V_L, gap, num_level, mu_plus, wf_L, H_data_Np1, H_data_N)
rate_R_2 = Rate(gamma_R, V_R, gap, num_level, mu_plus, wf_R, H_data_Np1, H_data_N)

my_tunnel_2 = Tunneling(rate_L_2, rate_R_2)

the_matrix_2 = my_tunnel_2.matrix()
print_matrix(the_matrix_2)


print("######################################################################")
print("RATE EQUATION SOLVING")
print("######################################################################")

solver_rate_eq = SolverRateEquation(my_tunnel_1, my_tunnel_2, param_particle, param_rate)

size_tot = H_data_N.size + H_data_Np1.size

init_vec = [0.25, 0.25, 0.25, 0.25, 0, 0, 0, 0, 0,0]
print("init_vec= ", init_vec)

I_array_left = [0]*num_T
I_array_right = [0]*num_T

thefile = open('current.txt', 'w')
for i in range(num_T):

    t1, y1, t2, y2, init_vec_interm = solver_rate_eq.calc_diag(init_vec, i)
    # 'init_vec_interm' is the vector of initial conditions for the 2nd half of period.
    current = Current(t1, y1, t2, y2)

    # The new way to calculate current - using vector operations
    start_time = time.time()
    I_av_left_new_1, I_av_right_new_1 = \
        current.calc_average(my_tunnel_1, init_vec, param_rate['T'], param_rate['tolerance'])
    I_av_left_new_2, I_av_right_new_2 = \
        current.calc_average(my_tunnel_2, init_vec_interm, param_rate['T'], param_rate['tolerance'])
    print(("--- %s seconds ---" % (time.time() - start_time)))
    I_av_left_new = (I_av_left_new_1 + I_av_left_new_2)/2
    I_av_right_new = (I_av_right_new_1 + I_av_right_new_2)/2

    thefile.write("{0} {1} {2}\n".format(i, I_av_left_new, I_av_right_new))
    # The old way to calculate current
    ###
    I_left_1, I_left_2, I_av_left = current.calc_left(my_tunnel_1, my_tunnel_2)
    I_right_1, I_right_2, I_av_right = current.calc_right(my_tunnel_1, my_tunnel_2)
    ###

    print("i, I_av_left_old, I_av_right_old, I_av_left_new=, I_av_right_new=  ", i, I_av_left, I_av_right,
        I_av_left_new, I_av_right_new)

    if i == num_T - 1:
        solver_rate_eq.plotting(t1, y1, t2, y2)
        single_level1, single_level2 = solver_rate_eq.single_level_calc(y1, y2)
        single_level_plot(t1, t2, single_level1, single_level2, param_rate)
        current_plot(I_left_1, I_left_2, I_right_1, I_right_2, t1, t2, param_rate, I_av_left_new, I_av_right_new)

    init_vec = y2[num_points - 1, :]
    print("i= ", i, "summ of all probabilities:", np.sum(init_vec))

thefile.close()
current_average_plot_file("current.txt", param_rate)



#####
'''

for i in range(num_T):
  thefile.write("{0} {1} {2}\n".format(i, I_array_left[i], I_array_right[i]))



'''



