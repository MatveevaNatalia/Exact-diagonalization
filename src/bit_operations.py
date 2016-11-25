import numpy as np

#Calculates all possible unique permutations
#of ground state

    #Beginning

class UniqueElement:
    def __init__(self, value, occurrences):
        self.value = value
        self.occurrences = occurrences


def perm_unique(elements):
    eset = set(elements)
    listunique = [UniqueElement(i, elements.count(i)) for i in eset]
    u = len(elements)
    return perm_unique_helper(listunique, [0] * u, u - 1)


def perm_unique_helper(listunique, result_list, d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d] = i.value
                i.occurrences -= 1
                for g in perm_unique_helper(listunique, result_list, d - 1):
                    yield g
                i.occurrences += 1

    # End


def to_integer(bitlist):
    # Transform a vector of zeros and ones to
    # an integer number
    # Ex: (0,1,0,1) -> 5
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit
    return out


def to_bitfield(n, size):
    # Transform an integer into numpy array, which contains its
    # binary representation. The dimension of the array equals to
    # 'size'.
    # Ex: size = 4
    # 5 -> (0,1,0,1)

    bitfield = [1 if digit == '1' else 0 for digit in bin(n)[2:]]
    if len(bitfield) < size:
        temp = [0] * (size - len(bitfield))
        bitfield = np.append(temp, bitfield)

    return np.array(bitfield)


def find_number_bits(index, state_bin, num_level):
    # Returns number of non-zero bits before a bit
    # with number 'index'.
    # Ex.:
    # For index = 5, state_bin = [1,1,0,1,0,1,1] returns 3.
    count = 0
    for i in range(index):
        mask_result = (state_bin >> (num_level-1-i)) & 1
        count += mask_result
    return count