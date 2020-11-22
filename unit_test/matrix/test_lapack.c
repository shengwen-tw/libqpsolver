#include <stdio.h>
#include <math.h>
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

	matrix_t A; //A is m-by-n matrix
	matrix_construct(&A, 4, 6, ELEMENTS(+1, +4, +0, +1, -3, +2,
	                                    +2, +8, +1, +1, -4, +6,
	                                    -1, -4, -1, +0, +1, -2,
	                                    +1, +4, +0, +1, -3, +1));
	matrix_t *Q = matrix_zeros(4, 4); //Q is m-by-m orthogonal matrix
	matrix_t *R = matrix_zeros(4, 6); //R is m-by-n upper triangular matrix

	matrix_qr_factorization(&A, &Q, &R);

	PRINT_MATRIX(A);
	PRINT_MATRIX(*Q);
	PRINT_MATRIX(*R);

	matrix_delete(Q);
	matrix_delete(R);
}

void test_solve_null_space(void)
{
	printf("\nsolve null space of A with QR decomposition:\n");

	matrix_t A;  //m-by-n matrix
	matrix_construct(&A, 4, 6, ELEMENTS(+1, +4, +0, +1, -3, +2,
	                                    +2, +8, +1, +1, -4, +6,
	                                    -1, -4, -1, +0, +1, -2,
	                                    +1, +4, +0, +1, -3, +1));
	matrix_t *At = matrix_new(6, 4); //n-by-m matrix
	matrix_transpose(&A, At);

	matrix_t *Q; //Q is n-n orthogonal matrix
	matrix_t *R; //R is n-m upper triangular matrix
	matrix_qr_factorization(At, &Q, &R); //decompose transpose(A) using QR factorization

	/* test the accuracy of QR factorization */
	matrix_t *tmp = matrix_new(6, 4);
	matrix_t *A_test = matrix_new(4, 6);
	matrix_multiply(Q, R, tmp);    //calculate QR
	matrix_transpose(tmp, A_test); //A_test = transpose(QR)

	/* calculate the zeros row count of R matrix */
	int n_zero_cols = 0;
	int r, c;
	for(r = (R->row - 1); r >= 0; r--) {
		FLOAT norm = 0;
		for(c = 0; c < R->column; c++) {
			norm += matrix_at(R, r, c) * matrix_at(R, r, c);
		}

		if(fabs(norm) < 1e-8) {
			n_zero_cols++;
		}
	}

	/* copy null space bias vectors from Q matrix */
	matrix_t *N_A = matrix_new(Q->row, n_zero_cols);
	for(r = 0; r < N_A->row; r++) {
		for(c = 0; c < n_zero_cols; c++) {
			matrix_at(N_A, r, c) = matrix_at(Q, r, (Q->column - n_zero_cols + c));
		}
	}

	printf("if A is very close to A_test then the result of Null(A) is correct:\n");
	PRINT_MATRIX(A);
	PRINT_MATRIX(*A_test);
	printf("QR decomposition of transpose(A):\n");
	PRINT_MATRIX(*Q);
	PRINT_MATRIX(*R);
	printf("Null space of A:\n");
	PRINT_MATRIX(*N_A);

	matrix_delete(A_test);
	matrix_delete(Q);
	matrix_delete(R);
}

void test_matrix_rank(void)
{
	printf("\nRank of matrix A:\n");

	matrix_t A;
	matrix_construct(&A, 4, 4, ELEMENTS(10,  0,  0,     0,
	                                    0, 25,  0,     0,
	                                    0,  0, 34,     0,
	                                    0,  0,  0, 1e-15));

	int rank = matrix_rank(&A);

	PRINT_MATRIX(A);
	printf("The rank of matrix A is %d\n", rank);
}

int main(void)
{
	test_solve_linear_system();
	test_matrix_inversion();
	test_matrix_multiplication();
	test_qr_factorization();
	test_solve_null_space();
	test_matrix_rank();

	return 0;
}
