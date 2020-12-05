'''
do not directly run this script, you should execute the unit test
by launching the "run_test.sh"
'''

import libqpsolver
import os
import numpy as np
from cvxopt import matrix, solvers

sol_diff_cnt = 0
curr_test_num = 0

verbose = False

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def quadprog_cvxopt(P, q, A=None, b=None, A_eq=None, b_eq=None, options=None):
    """
    qp solver provided by cvxopt package:

    Minimize    (1/2)(x.' * P * x) + (q.' * x)
    Subject to  A * x < b
    and         A_eq * x = b_eq
    """

    #verbose option
    options = solvers.options['show_progress'] = verbose

    #objective function
    P, q = matrix(P), matrix(q)

    #inequality constraint
    if (A is not None) and (b is not None):
        A, b = matrix(A), matrix(b)
    else:
        A, b = None, None

    #equality constraint
    if (A_eq is not None) and (b_eq is not None):
        A_eq, b_eq = matrix(A_eq), matrix(b_eq)
    else:
        A_eq, b_eq = None, None

    #solve qp
    sol = solvers.qp(P, q, A, b, A_eq, b_eq, options);
    return np.array(sol['x'])

def matrix_compare(mat1, mat2):
    if len(mat1) != len(mat2):
        print('dimension of matrices are not equal for comparasion')
        return False;

    for i in range(len(mat1)):
        error_percentage = abs(mat1[i] - mat2[i]) / abs(mat1[i])

        #tolerate 5% of error
        if error_percentage > 0.05:
            return False

    return True

def generate_symmetric_positive_definite_matrix(row, column, max_val):
    vec = np.random.rand(row, column)
    vec_transposed = np.transpose(vec)

    spd_matrix = np.matmul(vec, vec_transposed)
    spd_matrix = np.multiply(max_val, spd_matrix)

    return spd_matrix

def generate_random_row_vector(row, max_val):
    vec = np.random.rand(row, 1)
    vec = np.multiply(max_val, vec)

    return vec

def test_random_2x2_qp_problem(cost_func_max_val):
    P = generate_symmetric_positive_definite_matrix(2, 2, cost_func_max_val);
    q = generate_random_row_vector(2, cost_func_max_val)

    #randomlly turn on / off the inequality constraints
    A = None
    b = None
    if np.random.rand(1, 1) < 0.5:
        A = np.array([[+1.0, +1.0],
                      [-1.0, +2.0],
                      [+2.0, +1.0]])
        b = np.array([[2.0],
                      [2.0],
                      [3.0]])

    #randomlly turn on / off the equality constraints
    A_eq = None;
    b_eq = None;
    if np.random.rand(1, 1) < 0.5:
        A_eq = np.array([[1.0, 1.0]])
        b_eq = np.array([[0.0]])

    if verbose == True:
        print('[Test input matrices]')
        print('P = \n%s' %(P))
        print('q = \n%s' %(q))
        print('A = \n%s' %(A))
        print('b = \n%s' %(b))
        print('A_eq = \n%s' %(A_eq))
        print('b_eq = \n%s\n' %(b_eq))

    if verbose == True:
        print('[debug message from CVXOPT]')

    cvxopt_sol = quadprog_cvxopt(P, q, A, b, A_eq, b_eq, None)

    if verbose == True:
        print('\n[debug message from libqpsolver]')

    libqpsolver_sol = libqpsolver.quadprog(P, q, A, b, A_eq, b_eq, None, None)

    if verbose == True:
        print('\n[Optimal solution by CVXOPT]')
        print(cvxopt_sol)

        print('\n[Optimal solution by libqpsolver]')
        print(libqpsolver_sol)

    test_result = matrix_compare(cvxopt_sol, libqpsolver_sol)

    global sol_diff_cnt
    global curr_test_num
    curr_test_num = curr_test_num + 1

    if test_result == True:
        if verbose == True:
            print(f"{bcolors.OKGREEN}\n[unit test of #%d is passed]{bcolors.ENDC}" %(curr_test_num))
    else:
        if verbose == True:
            print(f"{bcolors.FAIL}\n[unit test of #%d is failed]{bcolors.ENDC}" %(curr_test_num))
        sol_diff_cnt = sol_diff_cnt + 1

    if verbose == True:
        print('===============================================================')

    return test_result

def test_random_NxN_qp_problem(N, cost_func_max_val):
    P = generate_symmetric_positive_definite_matrix(N, N, cost_func_max_val);
    q = generate_random_row_vector(N, cost_func_max_val)

    #randomlly turn on / off the inequality constraints
    A = None
    b = None
    if np.random.rand(1, 1) < 0.5:
        A = np.identity(N)
        b = np.zeros((N, 1))

        for i in range(N):
            while True:
                ran = np.random.rand(1, 1)
                if(ran > 0.1 and ran < 0.9):
                    b[i, 0] =  N * ran
                    break

    #randomlly turn on / off the equality constraints
    A_eq = None
    b_eq = None
    if np.random.rand(1, 1) < 0.5:
        A_eq = np.ones((1, N));
        b_eq = np.array([[0.0]])

    if verbose == True:
        print('[Test input matrices]')
        print('P = \n%s' %(P))
        print('q = \n%s' %(q))
        print('A = \n%s' %(A))
        print('b = \n%s' %(b))
        print('A_eq = \n%s' %(A_eq))
        print('b_eq = \n%s\n' %(b_eq))

    if verbose == True:
        print('[debug message from CVXOPT]')

    cvxopt_sol = quadprog_cvxopt(P, q, A, b, A_eq, b_eq, None)

    if verbose == True:
        print('\n[debug message from libqpsolver]')

    libqpsolver_sol = libqpsolver.quadprog(P, q, A, b, A_eq, b_eq, None, None)

    if verbose == True:
        print('\n[Optimal solution by CVXOPT]')
        print(cvxopt_sol)

        print('\n[Optimal solution by libqpsolver]')
        print(libqpsolver_sol)

    test_result = matrix_compare(cvxopt_sol, libqpsolver_sol)

    global sol_diff_cnt
    global curr_test_num
    curr_test_num = curr_test_num + 1

    if test_result == True:
        if verbose == True:
            print(f"{bcolors.OKGREEN}\n[unit test of #%d is passed]{bcolors.ENDC}" %(curr_test_num))
    else:
        if verbose == True:
            print(f"{bcolors.FAIL}\n[unit test of #%d is failed]{bcolors.ENDC}" %(curr_test_num))
        sol_diff_cnt = sol_diff_cnt + 1

    if verbose == True:
        print('===============================================================')

    return test_result


def test_libqpsolver():
    test_suite_exec_times = 1000

    for i in range(0, test_suite_exec_times):
        #test my specified 2x2 problems
        test_random_2x2_qp_problem(100)
        test_random_2x2_qp_problem(500)
        test_random_2x2_qp_problem(1000)
        test_random_2x2_qp_problem(10000)

        #test 3x3 problems with simple constraints
        test_random_NxN_qp_problem(3, 1000)
        test_random_NxN_qp_problem(3, 500)
        test_random_NxN_qp_problem(3, 1000)
        test_random_NxN_qp_problem(3, 10000)

        #test 15x15 problems with simple constraints
        test_random_NxN_qp_problem(15, 1000)
        test_random_NxN_qp_problem(15, 500)
        test_random_NxN_qp_problem(15, 1000)
        test_random_NxN_qp_problem(15, 10000)

        #test 50x50 problems with simple constraints
        test_random_NxN_qp_problem(50, 1000)

    total_test_times = test_suite_exec_times * 13

    correctness = (1.0 - (sol_diff_cnt / total_test_times)) * 100.0

    print(f"{bcolors.BOLD}unit test total run times = %d{bcolors.ENDC}" %(total_test_times))
    print(f"-> %d of %d failed" %(sol_diff_cnt, total_test_times))

    #if error count exceed 1% of total test count, than the solver is not stable
    if (sol_diff_cnt / total_test_times) > 0.01:
        print(f"{bcolors.FAIL}[failed] correctness = %f%%{bcolors.ENDC}" %(correctness))
        print(f"{bcolors.FAIL}the solver is unstable due to the correctness is lower than 99%{bcolors.ENDC}")
        exit(1)
    else:
        print(f"{bcolors.OKGREEN}[passed] correctness = %f%%{bcolors.ENDC}" %(correctness))
        exit(0)

def main():
    test_libqpsolver()

if __name__ == "__main__": main()
