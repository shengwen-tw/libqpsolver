'''
do not directly run this script, you should execute the unit test
by launching the "run_test.sh"
'''

import matrix_wrapper
import os
import numpy as np

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

def matrix_compare(mat1, mat2):
    if (mat1.shape[0] != mat2.shape[0]) or (mat1.shape[1] != mat2.shape[1]):
        print('dimension of matrices are not equal for comparasion')
        return False;

    for i in range(mat1.shape[0]):
        for j in range(mat1.shape[1]):
            error_percentage = abs(mat1[i][j] - mat2[i][j]) / abs(mat1[i][j])

        #tolerate 5% of error
        if error_percentage > 0.05:
            return False

    return True

def generate_random_matrix(M, N, magnitude):
    rand_mat = np.random.rand(M, N);

    for r in range(0, M):
        for c in range(0, N):
            rand_mat[r, c] = magnitude * rand_mat[r, c]

    return rand_mat

def unit_test_matrix_add(M, N, magnitude):
    mat1 = generate_random_matrix(M, N, magnitude);
    mat2 = generate_random_matrix(M, N, magnitude);

    np_result = np.add(mat1, mat2)
    my_result = matrix_wrapper.matrix_add(mat1, mat2);

    if verbose == True:
        print('mat1 = \n%s' %(mat1))
        print('mat2 = \n%s' %(mat2))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_add{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_add{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_add_by(M, N, magnitude):
    mat1 = generate_random_matrix(M, N, magnitude);
    mat2 = generate_random_matrix(M, N, magnitude);

    np_result = np.add(mat1, mat2)
    my_result = matrix_wrapper.matrix_add_by(mat1, mat2);

    if verbose == True:
        print('mat1 = \n%s' %(mat1))
        print('mat2 = \n%s' %(mat2))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_add_by{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_add_by{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_sub(M, N, magnitude):
    mat1 = generate_random_matrix(M, N, magnitude);
    mat2 = generate_random_matrix(M, N, magnitude);

    np_result = np.subtract(mat1, mat2)
    my_result = matrix_wrapper.matrix_sub(mat1, mat2);

    if verbose == True:
        print('mat1 = \n%s' %(mat1))
        print('mat2 = \n%s' %(mat2))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_sub{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_sub{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_multiply(M, N, magnitude):
    mat1 = generate_random_matrix(M, N, magnitude);
    mat2 = generate_random_matrix(M, N, magnitude);

    np_result = np.matmul(mat1, mat2)
    my_result = matrix_wrapper.matrix_multiply(mat1, mat2);

    if verbose == True:
        print('mat1 = \n%s' %(mat1))
        print('mat2 = \n%s' %(mat2))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_multiply{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_multiply{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_scale_by(M, N, magnitude):
    mat = generate_random_matrix(M, N, magnitude);

    scaler = 50

    np_result = scaler * mat 
    my_result = matrix_wrapper.matrix_scale_by(scaler, mat);

    if verbose == True:
        print('mat = \n%s' %(mat))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_scale_by{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_scale_by{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_scaling(M, N, magnitude):
    mat = generate_random_matrix(M, N, magnitude);

    scaler = 50

    np_result = scaler * mat 
    my_result = matrix_wrapper.matrix_scaling(scaler, mat);

    if verbose == True:
        print('mat = \n%s' %(mat))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_scaling{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_scaling{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_transpose(M, N, magnitude):
    mat = generate_random_matrix(M, N, magnitude);

    np_result = np.transpose(mat)
    my_result = matrix_wrapper.matrix_transpose(mat);

    if verbose == True:
        print('mat = \n%s' %(mat))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        print(f"{bcolors.OKGREEN}[passed] matrix_transpose{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_transpose{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_rank(M, N, magnitude):
    mat = generate_random_matrix(M, N, magnitude);

    np_result = np.linalg.matrix_rank(mat)

    my_result = matrix_wrapper.matrix_rank(mat);

    if verbose == True:
        print('mat = \n%s' %(mat))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if np_result == my_result:
        print(f"{bcolors.OKGREEN}[passed] matrix_rank{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_rank{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_all_functions():
    unit_test_matrix_add(3, 3, 100)
    unit_test_matrix_add_by(3, 3, 100)
    unit_test_matrix_sub(3, 3, 100)
    unit_test_matrix_multiply(3, 3, 100)
    unit_test_matrix_scaling(3, 3, 100)
    unit_test_matrix_scale_by(3, 3, 100)
    unit_test_matrix_transpose(3, 3, 100)
    unit_test_matrix_rank(3, 3, 100)

def main():
    #test_matrix_functions()
    unit_test_matrix_all_functions()

if __name__ == "__main__": main()
