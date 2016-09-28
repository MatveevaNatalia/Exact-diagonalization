from numpy import linalg as LA

class ParamModel():
    def __init__(self, Num_level, Num_part, Delta, Scat_ampl):
        self.Num_level = Num_level
        self.Num_part = Num_part
        self.Delta = Delta
        self.Scat_ampl = Scat_ampl
        

def Kin_Energy(params):
    '''
    Returns the array which contains single-particle energies
    '''
    energy = [0]*params.Num_level    
    energy = [i * params.Delta for i, __ in enumerate(energy)]
    
    return energy


def Ground_State(Num_level, Num_part):
    '''
    Fills the vector for ground state
    '''
    if Num_level < 0:
        raise ValueError("Number of levels must be non negative.")
    if Num_part < 0:
        raise ValueError("Number of particles must be non negative.")
    if Num_level < Num_part:
        raise ValueError("Number of levels is smaller than the number of particles.")
        
    v_ground = [0]*Num_level
    for i in range(Num_part):
        v_ground[i] = 1    
    return v_ground

'''
Calculates all possible unique permutations 
of ground state
'''
# Begining

class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

        
def perm_unique(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)


def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1

#End
                
def To_integer(bitlist):
    '''
    Transform a vector of zeros and ones to an integer number
    Ex: (0,1,0,1) -> 5
    '''
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit
    return out
 
def To_bitfield(n):
    '''
    Transform an integer into list which contains its binary representation.
    Ex: 5 -> (0,1,0,1)
    '''
    return [1 if digit=='1' else 0 for digit in bin(n)[2:]]
    
def Final_list_create(energy_arr, state_list, Delta):
    '''
    Creates the final list where 2 numbers correspond to each state.
    The first one is the kinetic energy and the second one is the 
    integer number which corresponds to the binary representation of the state.
    '''
    Num_level = len(state_list[0])
    
    list_energy = []
          
    for ss in state_list:
        state_bin = To_integer(list(ss))
        Energy = Energy_Calc(energy_arr, state_bin, Num_level, Delta) 
        list_energy.append([Energy, state_bin])
               
    return list_energy


def Energy_Calc(energy_arr, state_bin, Num_level, Delta):
    '''
    Calculates the kinetic energy for a given state
    '''
    Energy = 0
    for i in range(Num_level):
        # smaller bits correspond to higher energies 
        mask_res = (state_bin >> i)&1 
        Energy += mask_res * energy_arr[Num_level-1-i]
    return Energy

def Find_Number_Bits(index, state_bin, Num_level):
    '''
    Returns number of non-zero bits before a bit with number 'index'.
    Ex.: For index = 5, state_bin = [1,1,0,1,0,1,1] returns 3.
    '''
    count = 0
    for i in range(index):
        mask_result = (state_bin >> (Num_level-1-i))&1
        count += mask_result
    return count


def Kill_State(index, state_bin_init, Num_level):
    '''
    Annihilation operator c_{index}. Kills the particle
    at energy level which is marked by variable index.
    '''
    if(index >= Num_level):
        raise ValueError("Index is bigger than number of levels")

    if(index < 0 or state_bin_init < 0 or Num_level < 0):
        raise ValueError("Index or state integer or number of levels is negative!")   
        
    # Make also check if number of bits in state_bin_init <= Num_levels
    
    state_bin_fin = state_bin_init^(1<<(Num_level-1-index))
    if(state_bin_init < state_bin_fin):
        state_bin_fin = 0
        alive = 0.0
    else:
        alive = 1.0
        
    number_perm = Find_Number_Bits(index, state_bin_init, Num_level)
    coeff  = alive * (-1)**number_perm
        
    return (state_bin_fin, coeff)



def Create_State(index, state_bin_init, Num_level):
    '''
    Creation operator c^{+}_{index}. Creates the particle
    at energy level which is marked by variable index.
    '''
    
    if(index >= Num_level):
        raise ValueError("Index is bigger than number of levels")
        
    if(index < 0 or state_bin_init < 0 or Num_level < 0):
        raise ValueError("Index or state integer or number of levels is negative!")   
        
    # Make also check if number of bits in state_bin_init <= Num_levels
    
    state_bin_fin = state_bin_init^(1<<(Num_level-1-index))
    if(state_bin_init > state_bin_fin):
        state_bin_fin = 0
        coeff = 0.0
    else:
        coeff = 1.0
        
    number_perm = Find_Number_Bits(index, state_bin_init, Num_level)
    coeff  *= (-1)**number_perm
    
    return (state_bin_fin, coeff)


def Find_Scat_Elem(params, restr = True):
    '''
    Finds all nonzero matrix elements of scattering operator.
    If restr = True, assumes the condition: V_{ijkl} != 0 only if i+j = k+l
    If restr = False, doesn't assume this condition
    '''    

    list_main = []
    for i in range(0, params.Num_level):
        for j in range(i, params.Num_level):
            for k in range(0, params.Num_level):
                for l in range(k, params.Num_level):
                                        
                    cond = (i!=j and i!=k and i!=l and j!=k and j!=l and k!=l)
                    
                    if(restr):
                        cond = (cond and i+j == k+l)
                        
                    if(cond):
                                                
                        list_main.append([[i, j, k, l], [params.Scat_ampl] ])
                        
                        if(restr):
                            if(not Verification_Find_Scat_Elem([i, j, k, l])):
                                raise RuntimeError("More then one unique index in scattering element!")
    return list_main


def Verification_Find_Scat_Elem(V_list):
    '''
    This function checks if there are the same elements
    in every sublist of list_main.
    If yes, it returns check = 0
    '''
    
    check = True
    
    [i,j,k,l] = V_list
    
    cond_1 = (i == j or i == k or i == l)
    cond_2 = (j == k or j == l)
    cond_3 = (k == l)
        
    if(cond_1 or cond_2 or cond_3):
        check = False

    return check
        
def Set_Scat_Element(scat_elem, index, ampl):
    scat_elem[index][1][0] = ampl
    return scat_elem
    
    
    
def getKey(item):
    '''
    Supplimentary function for sorted() function. 
    This function is the standard library function 
    which performs the sorting of an array.
    '''
    return item[0]

    
    
def Fill_myMatrix(scat_elem, final_list_sorted, Num_level):
    '''
    #Fills the matrix for Hamiltonian. The electrostatic part is not included.
    '''
    
    Num_state = len(final_list_sorted)
    myMatrix =  [[0.0 for x in range(Num_state)] for y in range(Num_state)] 
   
    
    for [[i,j,k,l],[scat_ampl]] in scat_elem:        
        for i_stat, [__, state_initial] in enumerate(final_list_sorted):
                        
            coeff = 1.0
            (state_temp, coeff_temp) = Kill_State(l, state_initial, Num_level)
            coeff *= coeff_temp 

            if(coeff == 0):
                continue

            (state_temp, coeff_temp) = Kill_State(k, state_temp, Num_level)
            coeff *= coeff_temp 

            if(coeff == 0):
                continue

            (state_temp, coeff_temp) = Create_State(j, state_temp, Num_level)
            coeff *= coeff_temp    
            
            if(coeff == 0):
                continue
                
            (state_temp, coeff_temp) = Create_State(i, state_temp, Num_level)
            coeff *= coeff_temp  
            
            if(coeff == 0):
                continue
            
            i_col = Find_State_Index(final_list_sorted, state_temp)
            
            myMatrix[i_stat][i_col] += coeff * scat_ampl    
    
    for i in range(Num_state):
        for j in range(Num_state):
            if(i == j):
                myMatrix[i][j] = final_list_sorted[i][0]
    

    return myMatrix
    
        

def Find_State_Index(final_list_sorted, state_final):
    '''
    Given an integer number finds it's corresponding index.
    This operation is necessary when filling the matrix - 
    in order to put the matrix element in the correct raw and coulmn.
    '''

    for index, [__, state_basis] in enumerate(final_list_sorted):    
        if state_basis == state_final:
            return index

        
def printMatrix(mat): 
    '''
    Prints the matrix for Hamiltonian
    '''
    print("Matrix representation of Hamiltonian")
    print (" ","   ".join([str(x) for x in range(len(mat))]))
    for i,x in enumerate(mat):
        print (i," ".join([str(y) for y in x])) 
        
        
def SubspaceInfo(params, energy_arr, scat_elem, verbose_states = True, verbose_matrix = True):
    '''
    This function returns the information about N-subparticle space:
    the list "final_list_sorted" contains the basis vectors in binary repersentation 
    and corresponding energies, list "eigen_vectors" contains eigen vectors. 
    '''
    v_ground = Ground_State(params.Num_level, params.Num_part)
    state_list = list(perm_unique(v_ground))
    final_list = Final_list_create(energy_arr, state_list, params.Delta)
    final_list_sorted = sorted(final_list, key=getKey)

    if(verbose_states):
        print("Ground state vector:")
        print(v_ground)
        print("List of all states:")
        print(state_list)
        print("Number of all states")
        print(len(state_list))
        print("Final list of all states written as integers with corresponding kinetic energy:")
        print(final_list)
        print("Sorted final list:")
        print(final_list_sorted)
        print("List of all nonzero scattering elements:")
        print(scat_elem)
        
    myMatrix = Fill_myMatrix(scat_elem, final_list_sorted, params.Num_level)   
    if(verbose_matrix):    
        printMatrix(myMatrix)
    
    eigen_values, eigen_vectors = LA.eig(myMatrix)
    
    return final_list_sorted, eigen_vectors, myMatrix

def Eigen_vector_ket(basis, eigen_vectors, index):
    '''
    Index indicates the eigen state vector
    which we want to represent as ket vector.
    The function returns the list [[ff_0, e_0], [ff_1, e_1] ...[ff_N-1,e_N-1]],
    where ff_0, ff_1, ..., ff_N-1 are the components of eigen state vector,
    e_0, e_1, ..., e_N-1 are the basis states in binary representation.
    '''
    ket_vector = []
    
    for i, [__, ee] in enumerate(basis):
        ket_vector.append([eigen_vectors[i, index], ee])

    return ket_vector

def Kill_Eigen_State(ket_vector, params, index):
    '''
    Performs an action of annihilation operator c_index on an eigen state marked 
    with variable index_state.
    Ket vector is the following list: 
    [[ff_0, e_0], [ff_1, e_1] ...[ff_N-1,e_N-1]],
    where ff_0, ff_1, ..., ff_N-1 are the components of eigen state vector,
    e_0, e_1, ..., e_N-1 are the basis states in binary representation.
    
    '''            
    ket_vector_final = []
            
    for [ff, ee] in ket_vector:
        state_fin, coeff = Kill_State(index, ee, params.Num_level)
        ket_vector_final.append([ff*coeff, state_fin])
        
    return ket_vector_final    

def MastEqPrep(basis_Nplus1, eigen_Nplus1, params, basis_N, eigen_N, verbose):
    '''
    This function calculates 
    A_{ij} = \sum_{m=0}^{M-1}<\phi^N_i|\hat{c}_m|\psi^{N+1}_j>
    for given N and N+1.
    '''   
    
    print("N= ", params.Num_part -1, "N+1= ", params.Num_part) 
    
    # Loop over indexes for the left state
    for index_left in range(len(basis_N)):
         # Loop over indexes for the right state
        for index_right in range(len(basis_Nplus1)): 
            sum_total = 0
             # Loop over indexes for annihilation operator \hat{c}_index
            for index in range(params.Num_level):
                  
                ket_Nplus1 = Eigen_vector_ket(basis_Nplus1, eigen_Nplus1, index_right)        
                ket_final = Kill_Eigen_State(ket_Nplus1, params, index)
                
                # The complex conjugation in the left bra vector is not considered here,
                # because all elements of eigen states are choosen to be real.
                ket_N = Eigen_vector_ket(basis_N, eigen_N, index_left) 

                # Now we can calculate <phi_{index_left}|c_index|phi_{index_right}>
                scal_prod = 0
                                             
                for [ff, ee] in ket_N:
                    for [ff1, ee1] in ket_final:
                        if ee == ee1:
                            scal_prod += ff*ff1
                            
                sum_total += scal_prod

                if(verbose):                    
                    print("index= ", index, "i= ", index_left, "j= ", index_right)
                    print(ket_N)
                    print(ket_Nplus1)
                    print(ket_final)    
                    print("scal_prod= ", scal_prod)
            print("i= ", index_left, "j=", index_right, sum_total)
