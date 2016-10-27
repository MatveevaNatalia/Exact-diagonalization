from bit_operations import *

def kill_state(index, state_bin_init, num_level):
    # Annihilation operator c_{index}. Kills the particle
    # at energy level which is marked by variable index.
    if index >= num_level:
        raise ValueError("Index is bigger than number of levels")

    if index < 0 or state_bin_init < 0 or num_level < 0:
        raise ValueError("Index or state integer or number of levels is negative!")

        # To make also a check if
    # number of bits in state_bin_init <= Num_levels

    state_bin_fin = state_bin_init ^ (1 << (num_level - 1 - index))
    if state_bin_init < state_bin_fin:
        state_bin_fin = 0
        alive = 0.0
    else:
        alive = 1.0

    number_perm = find_number_bits(index, state_bin_init, num_level)
    coeff = alive * (-1) ** number_perm

    return state_bin_fin, coeff


def create_state(index, state_bin_init, num_level):
    # Creation operator c^{+}_{index}. Creates the particle
    # at energy level which is marked by variable index.
    if index >= num_level:
        raise ValueError("Index is bigger than number of levels")

    if index < 0 or state_bin_init < 0 or num_level < 0:
        raise ValueError("Index or state integer or number of levels is negative!")

        # Make also check if number of bits in state_bin_init <= Num_levels

    state_bin_fin = state_bin_init ^ (1 << (num_level - 1 - index))
    if state_bin_init > state_bin_fin:
        state_bin_fin = 0
        coeff = 0.0
    else:
        coeff = 1.0

    number_perm = find_number_bits(index, state_bin_init, num_level)
    coeff *= (-1) ** number_perm

    return state_bin_fin, coeff


def operator_c(action, index, bin_init, num_level):
    # This is the unique function for create_state and kill_state
    if action == "Create":
        bin_fin, coeff = create_state(index, bin_init, num_level)
    elif action == "Kill":
        bin_fin, coeff = kill_state(index, bin_init, num_level)
    else:
        raise ValueError("Unrecognized 'action' is given to operator_c function !")
    return bin_fin, coeff


def kill_eigen_state(ket_vector, num_level, index):
    # Performs an action of annihilation operator c_index on
    # an eigen state marked with variable index_state.
    # Ket vector is the following list:
    # [[ff_0, e_0], [ff_1, e_1] ...[ff_N-1,e_N-1]],
    # where ff_0, ff_1, ..., ff_N-1 are the components
    # of eigen state vector,
    # e_0, e_1, ..., e_N-1 are the basis states in
    # binary representation.
    ket_vector_final = []

    for [ff, ee] in ket_vector:
        state_fin, coeff = kill_state(index, ee, num_level)
        ket_vector_final.append([ff * coeff, state_fin])

    return ket_vector_final

def create_eigen_state(ket_vector, num_level, index):
    # Performs an action of annihilation operator c_index on
    # an eigen state marked with variable index_state.
    # Ket vector is the following list:
    # [[ff_0, e_0], [ff_1, e_1] ...[ff_N-1,e_N-1]],
    # where ff_0, ff_1, ..., ff_N-1 are the components
    # of eigen state vector,
    # e_0, e_1, ..., e_N-1 are the basis states in
    # binary representation.
    ket_vector_final = []

    for [ff, ee] in ket_vector:
        state_fin, coeff = create_state(index, ee, num_level)
        ket_vector_final.append([ff * coeff, state_fin])

    return ket_vector_final


def eigen_vector_ket(basis, eigen_vectors, index):
    # Index indicates the eigen state vector
    # which we want to represent as ket vector.
    # The function returns the list
    # [[ff_0, e_0], [ff_1, e_1] ...[ff_N-1,e_N-1]],
    # where ff_0, ff_1, ..., ff_N-1 are the components
    # of eigen state vector,
    # e_0, e_1, ..., e_N-1 are the basis states
    # in binary representation.

    ket_vector = []

    for i, [__, ee] in enumerate(basis):
        ket_vector.append([eigen_vectors[i, index], ee])

    return ket_vector


