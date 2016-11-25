import numpy as np
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

gamma_L = 5. # escape rate from the normal island into left electrode
gamma_R = 5. # escape rate from the normal island into right electrode
V_L = -2. # voltage on the left electrode
V_R = 2. # voltage on the right electrode

gap = 3. # superconducting gap

mu_minus = -2.
mu_plus = 2.

T = 1.
num_points = 1000
num_T = 50

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
    'num_T': num_T
}

#wf_L = [wf_amp("left", i) for i in range(num_level)]
#wf_R = [wf_amp("right", i) for i in range(num_level)]

wf_L = [1, np.sqrt(3), np.sqrt(5), np.sqrt(7)]
wf_R = [np.sqrt(2), np.sqrt(4), np.sqrt(6), np.sqrt(8)]

print("wf_L= ", wf_L)
print("wf_R= ", wf_R)

# Constructing the matrix for first half of period

rate_L_1 = Rate(gamma_L, V_L, gap, num_level, mu_minus, wf_L, H_data_Np1, H_data_N)
rate_R_1 = Rate(gamma_R, V_R, gap, num_level, mu_minus, wf_R, H_data_Np1, H_data_N)

my_tunnel_1 = Tunneling(rate_L_1, rate_R_1)


#the_matrix_1 = my_tunnel_1.matrix()
#print_matrix(the_matrix_1)

# Constructing the matrix for second half of period

rate_L_2= Rate(gamma_L, V_L, gap, num_level, mu_plus, wf_L, H_data_Np1, H_data_N)
rate_R_2 = Rate(gamma_R, V_R, gap, num_level, mu_plus, wf_R, H_data_Np1, H_data_N)

my_tunnel_2 = Tunneling(rate_L_2, rate_R_2)

#the_matrix_2 = my_tunnel_2.matrix()
#print_matrix(the_matrix_2)


print("######################################################################")
print("RATE EQUATION SOLVING")
print("######################################################################")

#init_vec = [0,0,0,0, 1/6., 1/6., 1/6., 1/6., 1/6., 1/6.]
#init_vec = [0.25, 0.25, 0.25, 0.25, 0, 0, 0, 0, 0, 0]

solver_rate_eq = SolverRateEquation(my_tunnel_1, my_tunnel_2, param_particle, param_rate)

size_tot = H_data_N.size + H_data_Np1.size
init_vec = np.arange(size_tot, dtype=float) + 1

t1, y1, t2, y2 = solver_rate_eq.calc_odeint(init_vec)

solver_rate_eq.plotting(t1, y1, t2, y2)





