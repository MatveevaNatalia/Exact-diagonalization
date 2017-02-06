import numpy as np
import sys
import os.path
sys.path.insert(0, "./src")
from numpy import linalg as LA

class Current():
    def __init__(self, t1, y1, t2, y2):
        self.dt1 = t1[1]-t1[0] # step on t for the 1st half of the period
        self.y1 = y1           # solution of rate equations for the 1st half of period
        self.dt2 = t2[1]-t2[0] # step on t for the 2nd half of the period
        self.y2 = y2           # solution of rate equations for the 2nd half of period

    def calc_left(self, my_tunnel_1, my_tunnel_2):
        # Calls the 'helper_left_vector()' method, which calculates
        # the LEFT current as the function of time (using vectorized formula)
        # and the average current (from rectangle method),
        # for the 1st and the 2nd half of the period.
        # Receives 'my_tunnel_1' and 'my_tunnel_2', which are the instances
        # of Tunneling class (Rate_equation.py)
        # for the 1st and the 2nd half of the period
        current_1, average_1 = self.helper_left_vector(my_tunnel_1, self.y1, self.dt1)
        current_2, average_2 = self.helper_left_vector(my_tunnel_2, self.y2, self.dt2)
        average = average_1 + average_2
        return current_1, current_2, average

    def calc_right(self, my_tunnel_1, my_tunnel_2):
        # Calls the 'helper_right_vector()' method, which calculates
        # the RIGHT current as the function of time (using vectorized formula)
        # and the average current (from rectangle method),
        # for the 1st and the 2nd half of the period.
        # Receives 'my_tunnel_1' and 'my_tunnel_2', which are the instances
        # of Tunneling class (Rate_equation.py)
        # for the 1st and the 2nd half of the period
        current_1, average_1 = self.helper_right_vector(my_tunnel_1, self.y1, self.dt1)
        current_2, average_2 = self.helper_right_vector(my_tunnel_2, self.y2, self.dt2)
        average = average_1 + average_2
        return current_1, current_2, average

    def helper_left_vector(self, my_tunnel, y, dt):
        # Calculates the left current as a function of time
        # as scalar product of J^L_i and y_i.
        # Also calculates the average current
        # NOT using the vectorized formula,
        # but just using the rectangle method of numerical integration
        # Input parameters: 'my_tunnel' is the instance of Tunneling class (Rate_equation.py),
        # 'y' is the solution of the rate equations at a given half of period,
        # 'dt' is the time step.

        num_points = len(y)
        current = [0] * num_points
        average = 0.
        for i in range(num_points):
            current[i] = np.dot(my_tunnel.cur_vec_left, y[i, :])
            average += current[i]
        return current, average * dt

    def helper_right_vector(self, my_tunnel, y, dt):
        # Calculates right current as a function of time
        # as scalar product of J^R_i and y_i.
        # Also calculates the average current
        # NOT using vectorized formula,
        # but just using rectangle method of numerical integration
        # Input parameters: 'my_tunnel' is the instance of Tunneling class (Rate_equation.py),
        # 'y' is the solution of the rate equations at a given half of period,
        # 'dt' is the time step.
        num_points = len(y)
        current = [0] * num_points
        average = 0.
        for i in range(num_points):
            current[i] = np.dot(my_tunnel.cur_vec_right, y[i, :])
            average += current[i]
        return current, average * dt

    def calc_average(self, my_tunnel, init_vec, period, tolerance):
        # Calculates the average current (both the left and right)
        # for a given half of period.
        # It uses the vectorized formula:
        # I_av = (2/period)*sum_{\lambda}((exp(\lambda*period/2) -1)/lambda)*c_{\lambda}*J_i*m_{i\lambda},
        # here J_i is the vector for the current,
        # \lambda are the eigen values of the rate equations matrix
        # c_{\lambda} is the vector of probabilities at t=0,
        # m_{i\lambda} is the matrix which consists of eigen vectors of the rate equations matrix

        # Input parameters: 'my_tunnel' is the instance of Tunneling class (Rate_equation.py),
        # 'init_vec' is the array containing the initial conditions,
        # 'period'  is the value of time period.



        the_matrix = my_tunnel.the_matrix
        lam, evec = LA.eig(the_matrix)
        coef = LA.solve(evec, init_vec)

        # Treating the zero eigen_values
        zero_val = np.abs(lam) < tolerance
        lam_vec = (np.exp(lam * period / 2) - 1) * coef
        lam_vec[np.logical_not(zero_val)] /= lam[np.logical_not(zero_val)]
        lam_vec[zero_val] = coef[zero_val]

        J_left = my_tunnel.cur_vec_left
        I_av_left = (2. / period) * np.dot(lam_vec, np.dot(J_left, evec))

        J_right = my_tunnel.cur_vec_right
        I_av_right = (2./period) * np.dot(lam_vec, np.dot(J_right, evec))

        return I_av_left, I_av_right

    '''
    # Below there are the functions for the current calculation
    # for the old approach (using only the loops and not the vectorized operations)

    def calc_left(self, rate_1, rate_2):
        # Calls the 'helper_left()' method, which calculates
        # the LEFT current as the function of time (using only loops)
        # and the average current (from rectangle method),
        # for the 1st and the 2nd half of the period.
        # Receives 'rate_1' and 'rate_2', which are the instances
        # of 'Rate' class (Rate_equation.py)
        # for the 1st and the 2nd half of the period

        current_1, average_1 = self.helper_left(rate_1, self.y1, self.dt1)
        current_2, average_2 = self.helper_left(rate_2, self.y2, self.dt2)
        average = average_1 + average_2
        return current_1, current_2, average

    def calc_right(self, rate_1, rate_2):
        # Calls the 'helper_right()' method, which calculates
        # the RIGHT current as the function of time (using only loops)
        # and the average current (from rectangle method),
        # for the 1st and the 2nd half of the period.
        # Receives 'rate_1' and 'rate_2', which are the instances
        # of 'Rate' class (Rate_equation.py)
        # for the 1st and the 2nd half of the period
        current_1, average_1 = self.helper_right(rate_1, self.y1, self.dt1)
        current_2, average_2 = self.helper_right(rate_2, self.y2, self.dt2)
        average = average_1 + average_2
        return current_1, current_2, average

    def helper_left(self, rate, y, dt):
        # Calculates the left current as a function of time (using only loops).
        # Also calculates the average current
        # NOT using the vectorized formula,
        # but just using the rectangle method of numerical integration
        # Input parameters: 'rate' is the instance of Rate class (Rate_equation.py),
        # 'y' is the solution of the rate equations at a given half of period,
        # 'dt' is the time step.

        size_n = rate.H_data_n.size
        size_np1 = rate.H_data_np1.size
        num_points = len(y)
        current = [0]*num_points
        average = 0.
        for i in range(num_points):
            for alpha in range(size_n):
                for beta in range(size_np1):
                    a_ba = rate.addition_rate(beta, alpha)
                    r_ab = rate.removal_rate(alpha, beta)
                    current[i] += a_ba * y[i, alpha] - r_ab * y[i, size_n + beta]
                    #average += current[i]  #was before, it is an error
            average += current[i]
        return current, average*dt

    def helper_right(self, rate, y, dt):
        # Calculates the right current as a function of time (using only loops).
        # Also calculates the average current
        # NOT using the vectorized formula,
        # but just using the rectangle method of numerical integration
        # Input parameters: 'rate' is the instance of Rate class (Rate_equation.py),
        # 'y' is the solution of the rate equations at a given half of period,
        # 'dt' is the time step.
        size_n = rate.H_data_n.size
        size_np1 = rate.H_data_np1.size
        num_points = len(y)
        current = [0] * num_points
        average = 0.
        for i in range(num_points):
            for alpha in range(size_n):
                for beta in range(size_np1):
                    r_ab = rate.removal_rate(alpha, beta)
                    a_ba = rate.addition_rate(beta, alpha)
                    current[i] += r_ab * y[i, size_n + beta] - a_ba * y[i, alpha]
                    #average += current[i] #was before, it is an error
            average += current[i]
        return current, average*dt
    '''

