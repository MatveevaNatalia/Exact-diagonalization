def Kin_Energy(level_index, Delta):
    '''
    Calculates the kinetic energy 
    for a give energy level
    '''
    return Delta*level_index


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
    for i in range(0, Num_part):
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
    Transform a vector of zeros and ones 
    to an integer number
    Ex: (0,1,0,1) -> 5
    '''
    out = 0
    for bit in bitlist:
        out = (out << 1) | bit
    return out
 
    
def Final_list_create(state_list, Delta):
    '''
    Creates the final list where 2 numbers correspond to each state.
    The first one is the kinetic energy and the second one is the 
    integer number which corresponds to the binary representation of the state.
    '''
    Num_level = len(state_list[0])
    final_list = [0]*len(state_list)
    for i_stat in range(0, len(state_list)):
        state_bin = To_integer(list(state_list[i_stat]))
        Energy = Energy_Calc(state_bin, Num_level, Delta)    
        temp = [Energy, state_bin]
        final_list[i_stat] = temp
    return final_list


def Energy_Calc(state_bin, Num_level, Delta):
    '''
    Calculates the kinetic energy for a given state
    '''
    Energy = 0
    for i in range(Num_level):
        mask_res = (state_bin >> i)&1
        Energy += mask_res * Kin_Energy(Num_level-1-i, Delta)
    return Energy


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
        coeff = 0.0
    else:
        coeff = 1.0
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
    return (state_bin_fin, coeff)


def Find_Scat_Elem(Num_level):
    '''
    Finds all nonzero matrix elements of scattering operator.
    Assumes the condition: V_{ijkl} != 0 only if i+j = k+l
    '''
    
    list_main=[[]]*0
    for i in range(0, Num_level):
        for j in range(i, Num_level):
            for k in range(0, Num_level):
                for l in range(k, Num_level):
                    if(i!=j and i!=k and i!=l and j!=k and j!=l and k!=l and (i+j)==(k+l)):
                        index_list = [i, j, k, l]
                        list_temp = [0]*4
                        for index in range(4):
                            list_temp[index] = index_list[index]
                        list_main.append(list_temp) 
                        
    if(Verification_Find_Scat_Elem(list_main) == 0):
        raise RuntimeError("More then one unique index in scattering element!")
    return list_main


def Verification_Find_Scat_Elem(list_main):
    '''
    This function checks if there are the same elements
    in every sublist of list_main.
    If yes, it returns check = 0
    '''
    
    check = 1
    for i in range(len(list_main)):
        el_0 = list_main[i][0]
        el_1 = list_main[i][1]
        el_2 = list_main[i][2]
        el_3 = list_main[i][3]
        cond_1 = (el_0 == el_1 or el_0 == el_2 or el_0 == el_3)
        cond_2 = (el_1 == el_2 or el_1 == el_3)
        cond_3 = (el_2 == el_3)
        
        if(cond_1 or cond_2 or cond_3):
            check = 0
    return check
        
        


def getKey(item):
    '''
    Supplimentary function for sorted() function. 
    This function is the standard library function 
    which performs the sorting of an array.
    '''
    return item[0]


def Fill_myMatrix(scat_elem, final_list_sorted, Num_level, scat_ampl):
    '''
    Fills the matrix for Hamiltonian. The electrostatic part is not included.
    '''
    
    Num_state = len(final_list_sorted)
    myMatrix =  [[0.0 for x in range(Num_state)] for y in range(Num_state)] 
   
    
    for i_elem in range(len(scat_elem)):
       
        for i_stat in range(Num_state):
            coeff = 1.0
            (state_temp, coeff_temp) = Kill_State(scat_elem[i_elem][3], final_list_sorted[i_stat][1], Num_level)
            coeff *= coeff_temp 
            #print(i_stat,'state_temp_1= ', state_temp, 'coeff_temp_1 =', coeff_temp, coeff)
            if(coeff == 1):
                (state_temp, coeff_temp) = Kill_State(scat_elem[i_elem][2], state_temp, Num_level)
                coeff *= coeff_temp 
                #print(i_stat, 'state_temp_2= ', state_temp,'coeff_temp_2 =', coeff_temp, coeff)
                if(coeff == 1):
                    (state_temp, coeff_temp) = Create_State(scat_elem[i_elem][1], state_temp, Num_level)
                    coeff *= coeff_temp 
                    #print(i_stat, 'state_temp_3= ', state_temp,'coeff_temp_3 =', coeff_temp, coeff)
                    if(coeff == 1):
                        (state_temp, coeff_temp) = Create_State(scat_elem[i_elem][0], state_temp, Num_level)
                        coeff *= coeff_temp 
                        #print(i_stat, 'state_temp_4= ', state_temp,'coeff_temp_4 =', coeff_temp, coeff)
                        if(coeff == 1):
                            i_col = Find_State_Index(final_list_sorted, state_temp)
                            myMatrix[i_stat][i_col] += coeff * scat_ampl

    # Fill Energy
    
    for i in range(Num_state):
        for j in range(Num_state):
            if(i == j):
                myMatrix[i][j] = final_list_sorted[i][0]
    
    printMatrix(myMatrix)
    return myMatrix


def Find_State_Index(final_list_sorted, state_bin):
    '''
    Given an integer number finds it's corresponding index.
    This operation is necessary when filling the matrix - 
    in order to put the matrix element in the correct raw and coulmn.
    '''
    for i in range(len(final_list_sorted)):
        if(final_list_sorted[i][1] == state_bin):
            return i

        
def printMatrix(mat): 
    '''
    Prints the matrix for Hamiltonian
    '''
    print (" ","   ".join([str(x) for x in range(len(mat))]))
    for i,x in enumerate(mat):
        print (i," ".join([str(y) for y in x])) 