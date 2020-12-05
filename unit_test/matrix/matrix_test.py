'''
do not directly run this script, you should execute the unit test
by launching the "run_test.sh"
'''

import matrix_wrapper
import os
import numpy as np

def test_matrix_functions():
    ########################
    # solver linear system #
    ########################

    print('solve AX=B:')

    A = np.array([[3.0, 5.0, 2.0],
                  [2.0, 1.0, 3.0],
                  [4.0, 3.0, 2.0]])
    B = np.array([[57, 23],
                  [22, 12],
                  [41, 84]])

    X = matrix_wrapper.solve_linear_system(A, B);

    print('A = \n%s' %(A))
    print('B = \n%s' %(B))
    print('X = \n%s' %(X))

    ####################
    # matrix inversion #
    ####################
    print('\nsolve matrix inversion inv(A):')

    A = np.array([[ 1.0, 0.0,  2.0],
                  [-1.0, 5.0,  0.0],
                  [ 0.0, 3.0, -9.0]])

    A_inv = matrix_wrapper.matrix_inverse(A)

    print('A = \n%s' %(A))
    print('A_inv = \n%s' %(A_inv))

    #########################
    # matrix multiplication #
    #########################

    print('\nsolve matrix multiplication B = AX:')

    A = np.array([[3, 5, 2],
                  [2, 1, 3],
                  [4, 3, 2]])
    X = np.array([[2.0,   38.39],
                  [9.00, -11.30],
                  [3.00, -17.82]])
    B = matrix_wrapper.matrix_multiply(A, B)

    print('A = \n%s' %(A))
    print('X = \n%s' %(X))
    print('B = \n%s' %(B))

    ####################
    # QR factorization #
    ####################

    print("\nsolve QR factorization:")

    A = np.array([[+1, +4, +0, +1, -3, +2],
                  [+2, +8, +1, +1, -4, +6],
                  [-1, -4, -1, +0, +1, -2],
                  [+1, +4, +0, +1, -3, +1]])

    Q = matrix_wrapper.qr_factorization_Q(A)
    R = matrix_wrapper.qr_factorization_R(A)

    print('A = \n%s' %(A))
    print('Q = \n%s' %(Q))
    print('R = \n%s' %(R))

    ###########################
    # solve matrix null space #
    ###########################

    print('\nRank of matrix A:');

    A = np.array([[10, 0,  0,     0],
                  [0, 25,  0,     0],
                  [0,  0, 34,     0],
                  [0,  0,  0, 1e-15]])
    rank = matrix_wrapper.matrix_rank(A)

    print('A = \n%s' %(A))
    print('The rank of matrix A is %d' %(rank));

def main():
    test_matrix_functions()

if __name__ == "__main__": main()
