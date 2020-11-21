#!/bin/bash

'''
do not directly run this script, you should execute the unit test
by launching the "run_test.sh"
'''

import libqpsolver
import os
import numpy as np
from cvxopt import matrix, solvers

os.system('soruce ./LD_PRELOAD.sh')

def quadprog_cvxopt(P, q, A=None, b=None, A_eq=None, b_eq=None, options=None):
    """
    qp solver provided by cvxopt package:

    Minimize    (1/2)(x.' * P * x) + (q.' * x)
    Subject to  A * x < b
    and         A_eq * x = b_eq
    """

    #verbose option
    options = solvers.options['show_progress'] = False

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
    return solvers.qp(P, q, A, b, A_eq, b_eq, options);

def matrix_compare(mat1, mat2):
    epsilon = 1e-3;

    if len(mat1) != len(mat2):
        print('dimension of matrices are not equal for comparasion')
        return False;

    for i in range(len(mat1)):
        abs_diff = abs(mat1[i] - mat2[i])
        print(abs_diff)
        if abs_diff > epsilon:
            return False

    return True

def test_qp():
    P = np.array([[+1.0, -1.0],
                  [-1.0, +2.0]])
    q = np.array([[-2.0],
                  [-6.0]])
    A = np.array([[+1.0, +1.0],
                  [-1.0, +2.0],
                  [+2.0, +1.0]])
    b = np.array([[2.0],
                  [2.0],
                  [3.0]])
    A_eq = np.array([[1.0, 1.0]])
    b_eq = np.array([[0.0]])

    cvxopt_sol = quadprog_cvxopt(P, q, A, b, A_eq, b_eq, None)
    cvxopt_sol = np.array(cvxopt_sol['x'])

    print('Test matrices:')

    print('P = \n%s' %(P))
    print('q = \n%s' %(q))
    print('A = \n%s' %(A))
    print('b = \n%s' %(b))
    print('A_eq = \n%s' %(A_eq))
    print('b_eq = \n%s' %(b_eq))

    print('\nOptimal solution given by cvxopy:')
    print(cvxopt_sol)

    libqpsolver_sol = libqpsolver.quadprog(P, q, A, b, A_eq, b_eq, None, None)
    print(libqpsolver_sol)

    if matrix_compare(cvxopt_sol, libqpsolver_sol) == True:
        print('\nunit test passed')
    else:
        print('\nerror, unit test did not passed');

test_qp()
