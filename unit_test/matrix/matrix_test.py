'''
do not directly run this script, you should execute the unit test
by launching the "run_test.sh"
'''

import matrix_wrapper
import os
import numpy as np

def test_matrix_functions():
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

def main():
    test_matrix_functions()

if __name__ == "__main__": main()
