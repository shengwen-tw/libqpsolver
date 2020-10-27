#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

void test_solve_linear_system(void)
{
	printf("solve AX=B:\n");

	/* solve Ax = b */
	matrix_t A, B, X;
	matrix_construct(&A, 3, 3, ELEMENTS(3, 5, 2,
	                                    2, 1, 3,
	                                    4, 3, 2));

	matrix_construct(&B, 3, 2, ELEMENTS(57, 23,
	                                    22, 12,
	                                    41, 84));

	matrix_construct(&X, 3, 2, ELEMENTS(0, 0,
	                                    0, 0,
	                                    0, 0));

	solve_linear_system(&A, &X, &B);

	PRINT_MATRIX(A);
	PRINT_MATRIX(B);
	PRINT_MATRIX(X);
}

void test_matrix_inversion(void)
{
	printf("\nsolve matrix inversion inv(A):\n");

	matrix_t A, A_inv;
	matrix_construct(&A, 3, 3, ELEMENTS( 1,  0,  2,
	                                     -1,  5,  0,
	                                     0,  3, -9));

	matrix_construct(&A_inv, 3, 3, ELEMENTS(0, 0, 0,
	                                        0, 0, 0,
	                                        0, 0, 0));

	matrix_inverse(&A, &A_inv);
	PRINT_MATRIX(A);
	PRINT_MATRIX(A_inv);
}

void test_matrix_multiplication(void)
{
	printf("\n3.solve matrix multiplication B = AX:\n");

	matrix_t A, B, X;
	matrix_construct(&A, 3, 3, ELEMENTS(3, 5, 2,
	                                    2, 1, 3,
	                                    4, 3, 2));

	matrix_construct(&B, 3, 2, ELEMENTS(0));

	matrix_construct(&X, 3, 2, ELEMENTS(2.0,   38.39,
	                                    9.00, -11.30,
	                                    3.00, -17.82));

	matrix_multiply(&A, &X, &B);
	PRINT_MATRIX(A);
	PRINT_MATRIX(X);
	PRINT_MATRIX(B);
}

void test_qr_factorization(void)
{
	printf("\n4.solve QR factorization:\n");

	matrix_t A, Q, R;
	matrix_construct(&A, 4, 4, ELEMENTS(16,  2,  3, 13,
	                                    5, 11, 10,  8,
	                                    9,  7,  6, 12,
	                                    4, 14, 15,  1));
	matrix_construct(&Q, 4, 4, ELEMENTS(0));
	matrix_construct(&R, 4, 4, ELEMENTS(0));

	matrix_qr_factorization(&A, &Q, &R);

	PRINT_MATRIX(A);
	PRINT_MATRIX(Q);
	PRINT_MATRIX(R);
}

int main(void)
{
	test_solve_linear_system();
	test_matrix_inversion();
	test_matrix_multiplication();
	test_qr_factorization();

	return 0;
}
