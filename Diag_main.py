import numpy as np
from numpy import linalg as LA
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys
import os.path
sys.path.insert(0, "./src")
from Func import *
from Rate_equation import *
from printing_functions import *

num_level = 4
delta = 1.0
scat_ampl = 1.0
num_part = 2

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


#print("######################################################################")
#print("MATRIX ELEMENTS FOR MASTER EQUATION")
#print("######################################################################")
#RemovalRate(H_data_Np1.basisFock, H_data_Np1.eigen_vectors, params_Nplus1, basis_N, eigen_system_N.eigen_vectors, False)

print("######################################################################")
print("RATE EQUATION PREPARATION")
print("######################################################################")

gamma_L = 1.
gamma_R = 1.
V_L = 5
V_R = 3
gap = 1.5

mu_minus = -3.
mu_plus = 0.1
period = 1.

#mu = Mu(mu_minus,mu_plus, period)

#wf_L = [wf_amp("left", i) for i in range(num_level)]

wf_L = [1, np.sqrt(3), np.sqrt(5), np.sqrt(7)]

print(wf_L)

rate_L = Rate(gamma_L, V_L, gap, num_level, mu_minus, wf_L, H_data_Np1, H_data_N)

#wf_R = [wf_amp("right", i) for i in range(num_level)]

wf_R = [0,0,0,0]

print(wf_R)

rate_R = Rate(gamma_R, V_R, gap, num_level, mu_minus, wf_R, H_data_Np1, H_data_N)

my_tunnel = Tunneling(rate_L, rate_R)

the_matrix = my_tunnel.matrix()
print_matrix(the_matrix)

init_vec = [0,0,0,0, 1/6., 1/6., 1/6., 1/6., 1/6., 1/6.]

#init_vec = [0.25, 0.25, 0.25, 0.25, 0, 0, 0, 0, 0, 0]

t = np.linspace(0, period/2.0, 1000)

y = odeint(my_tunnel.der_func, init_vec, t, args=(the_matrix,))

plt.plot(t, y[:, 0], 'r-', label='p0 ')
plt.plot(t, y[:, 1], '-', color='hotpink', label='p1 ')
plt.plot(t, y[:, 2], '-', color='orange', label='p2 ')
plt.plot(t, y[:, 3], 'b-', label='p3')
plt.plot(t, y[:, 4], '-', color='black', label='q0 ')
plt.plot(t, y[:, 8], '-', color='magenta', label='q4 ')
plt.plot(t, y[:, 9], 'g-', label='q5 ')
plt.legend(loc='upper left', fontsize=14, numpoints=1)
plt.show()

