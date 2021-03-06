'''
do not directly run this script, you should execute the unit test
by launching the "run_test.sh"
'''

import matrix_wrapper
import os
import time
import progressbar
import numpy as np

#show detailed unit test message
verbose = False

#unit test run time and test item numbers
test_suite_exec_times = 50
test_suite_items = 5

progress_cnt_max = test_suite_exec_times * test_suite_items
progress =  \
    progressbar.ProgressBar(maxval=progress_cnt_max, \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])

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
    rand_mat = np.random.rand(M, N)

    for r in range(0, M):
        for c in range(0, N):
            rand_mat[r, c] = magnitude * rand_mat[r, c]

    return rand_mat

def generate_random_psd_matrix(M, N, magnitude):
    A = None
    while True:
        #randomly generate a matrix until it is full rank
        A = generate_random_matrix(M, N, magnitude)
        if np.linalg.matrix_rank(A) == M:
            break;

    A_transposed = np.transpose(A)

    psd_mat = np.matmul(A, A_transposed)

    return psd_mat

def unit_test_matrix_inverse(M, N, magnitude):
    mat = generate_random_psd_matrix(M, N, magnitude);

    np_result = np.linalg.inv(mat)
    my_result = matrix_wrapper.matrix_inverse(mat);

    if verbose == True:
        print('mat = \n%s' %(mat))
        print('np_result = \n%s' %(np_result))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(np_result, my_result) == True:
        if verbose == True:
            print(f"{bcolors.OKGREEN}[passed] matrix_inverse{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_inverse{bcolors.ENDC}")
        exit(1)

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
        if verbose == True:
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
        if verbose == True:
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
        if verbose == True:
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
        if verbose == True:
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
        if verbose == True:
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
        if verbose == True:
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
        if verbose == True:
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
        if verbose == True:
            print(f"{bcolors.OKGREEN}[passed] matrix_rank{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] matrix_rank{bcolors.ENDC}")
        exit(1)

def unit_test_solve_linear_system(M, N, magnitude):
    A = generate_random_psd_matrix(M, N, magnitude);
    X = generate_random_matrix(M, N, magnitude);
    B = np.matmul(A, X)

    my_result = matrix_wrapper.solve_linear_system(A, B);
    B_check = np.matmul(A, my_result)

    if verbose == True:
        print('np_result = \n%s' %(X))
        print('my_result = \n%s' %(my_result))

    if matrix_compare(B, B_check) == True:
        if verbose == True:
            print(f"{bcolors.OKGREEN}[passed] solve_linear_system{bcolors.ENDC}")
        return
    else:
        print(f"{bcolors.FAIL}[failed] solve_linear_system{bcolors.ENDC}")
        exit(1)

def unit_test_matrix_all_functions(N, magnitude):
    if verbose == True:
        print(f"{bcolors.BOLD}unit test with %dx%d matrices, scale=%d{bcolors.ENDC}" \
              %(N, N, magnitude));
    unit_test_matrix_inverse(N, N, magnitude)
    unit_test_matrix_add(N, N, magnitude)
    unit_test_matrix_add_by(N, N, magnitude)
    unit_test_matrix_sub(N, N, magnitude)
    unit_test_matrix_multiply(N, N, magnitude)
    unit_test_matrix_scaling(N, N, magnitude)
    unit_test_matrix_scale_by(N, N, magnitude)
    unit_test_matrix_transpose(N, N, magnitude)
    unit_test_matrix_rank(N, N, magnitude)
    unit_test_solve_linear_system(N, N, magnitude)

def main():
    print(f"{bcolors.BOLD}start the unit test of linear algebra functions{bcolors.ENDC}")
    print('test items: %d' %(test_suite_exec_times * test_suite_items))

    #progress bar
    progress.start()
    progress.update(0)

    time_start = time.time()
    time_last = time_start

    #test_matrix_functions()
    for i in range(test_suite_exec_times):
        time_now = time.time()

        #update the progress bar every 10 seconds
        if (time_now - time_last) > 10:
            time_last = time_now
            progress.update(i * test_suite_items)
            print('\nelapsed time: %d seconds' %(time_now - time_start))

        for N in range(2, 51): #test from 2x2 to 50x50
            unit_test_matrix_all_functions(N, 10)
            unit_test_matrix_all_functions(N, 100)
            unit_test_matrix_all_functions(N, 1000)
            unit_test_matrix_all_functions(N, 10000)
            unit_test_matrix_all_functions(N, 100000)

    progress.finish()
    time_now = time.time()
    print('elapsed time: %d seconds' %(time_now - time_start))

    print(f"{bcolors.OKGREEN}[all passed]{bcolors.ENDC}")

if __name__ == "__main__": main()
