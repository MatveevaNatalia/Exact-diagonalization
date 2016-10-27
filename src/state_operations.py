from bit_operations import *
from quant_operators import *

def ground_state(num_level, num_part):
    # Fills the vector for ground state
    if num_level < 0:
        raise ValueError("Number of levels must be non negative.")
    if num_part < 0:
        raise ValueError("Number of particles must be non negative.")
    if num_level < num_part:
        raise ValueError("Number of levels is smaller than the number of particles.")

    v_ground = [0] * num_level
    for i in range(num_part):
        v_ground[i] = 1
    return v_ground


def final_list_create(energy_arr, state_list):
    # Creates the final list where 2 numbers correspond
    # to each state. The first one is the kinetic energy
    # and the second one is the integer number which
    # corresponds to the binary representation of the state.
    num_level = len(state_list[0])
    list_energy = []

    for ss in state_list:
        state_bin = to_integer(list(ss))
        energy = energy_calc(energy_arr, state_bin, num_level)
        list_energy.append([energy, state_bin])

    return list_energy


def energy_calc(energy_arr, state_bin, num_level):
    # Calculates the kinetic energy for a given state
    energy = 0
    for i in range(num_level):
        # smaller bits correspond to higher energies
        mask_res = (state_bin >> i) & 1
        energy += mask_res * energy_arr[num_level - 1 - i]
    return energy


def verification_find_scat_elem(v_list):
    # This function checks if there are the same elements
    # in every sublist of list_main.
    # If yes, it returns check = 0
    check = True

    [i, j, k, l] = v_list

    cond_1 = (i == j or i == k or i == l)
    cond_2 = (j == k or j == l)
    cond_3 = (k == l)

    if cond_1 or cond_2 or cond_3:
        check = False

    return check


def set_scat_element(scat_elem, index, ampl):
    scat_elem[index][1] = ampl
    return scat_elem


def get_key(item):
    # Supplementary function for sorted() function.
    # This function is the standard library function
    # which performs the sorting of an array.
    return item[0]


def fill_matrix_helper(state_initial, scat_index, num_level):
    # This function helpes to fill the main matrix.
    # It performs action of operators  C^+_iC^+_jC_kC_l
    # on a given basis state.
    coeff = 1.0
    state_temp = state_initial

    for ii in range(3, -1, -1):
        if ii >= 2:
            action = "Kill"
        else:
            action = "Create"
        state_temp, coeff_temp = operator_c(action, scat_index[ii], state_temp, num_level)
        coeff *= coeff_temp

    return state_temp, coeff


def find_state_index(final_list_sorted, state_final):
    # Given an integer number finds it's corresponding index.
    # This operation is necessary when filling the matrix -
    # in order to put the matrix element in the correct
    # raw and column.
    for index, [__, state_basis] in enumerate(final_list_sorted):
        if state_basis == state_final:
            return index
