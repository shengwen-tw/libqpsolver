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

def quadprog_cvxopt(P, q, G=None, h=None, A=None, b=None, options=None):
    """
    qp solver provided by cvxopt package:

    Minimize    (1/2)(x.' * P * x) + (q.' * x)
    Subject to  Gx < h
    and         Ax = b
    """

    #verbose option
    options = solvers.options['show_progress'] = False

    #objective function
    P, q = matrix(P), matrix(q)

    #inequality constraint
    if (G is not None) and (h is not None):
        G, h = matrix(G), matrix(h)
    else:
        G, h = None, None

    #equality constraint
    if (A is not None) and (b is not None):
        A, b = matrix(A), matrix(b)
    else:
        A, b = None, None

    #solve qp
    return solvers.qp(P, q, G, h, A, b, options);

def test_qp():
    P = np.array([[+1.0, -1.0],
                  [-1.0, +2.0]])
    q = np.array([[-2.0],
                  [-6.0]])
    G = np.array([[+1.0, +1.0],
                  [-1.0, +2.0],
                  [+2.0, +1.0]])
    h = np.array([[2.0],
                  [2.0],
                  [3.0]])
    A = np.array([[1.0, 1.0]])
    b = np.array([[0.0]])

    cvxopt_sol = quadprog_cvxopt(P, q, G, h, A, b, None)

    print('Test matrices:')

    print('P = \n%s' %(P))
    print('q = \n%s' %(q))
    print('G = \n%s' %(G))
    print('h = \n%s' %(h))
    print('A = \n%s' %(A))
    print('b = \n%s' %(b))

    print('\nOptimal solution given by cvxopy:')
    print(cvxopt_sol['x'])

    libqpsolver_sol = libqpsolver.quadprog(P, q, G, h, A, b, None, None)
    print(libqpsolver_sol)

test_qp()
