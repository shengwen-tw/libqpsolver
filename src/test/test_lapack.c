#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

int main(void)
{
	printf("1.solve linear system AX=B:\n");

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

	/* solve matrix inversion */
	printf("\n2.solve matrix inversion inv(M):\n");

	matrix_t M, M_inv;
	matrix_construct(&M, 3, 3, ELEMENTS( 1,  0,  2,
	                                     -1,  5,  0,
	                                     0,  3, -9));

	matrix_construct(&M_inv, 3, 3, ELEMENTS(0, 0, 0,
	                                        0, 0, 0,
	                                        0, 0, 0));

	matrix_inverse(&M, &M_inv);
	PRINT_MATRIX(M);
	PRINT_MATRIX(M_inv);

	/* solve matrix multiplication */
	printf("\n3.solve matrix multiplication B = AX:\n");
	matrix_multiply(&A, &X, &B);
	PRINT_MATRIX(A);
	PRINT_MATRIX(X);
	PRINT_MATRIX(B);

	/* solve QR factorization */
	matrix_t H;
	matrix_construct(&H, 4, 4, ELEMENTS(16,  2,  3, 13,
	                                    5, 11, 10,  8,
	                                    9,  7,  6, 12,
	                                    4, 14, 15,  1));
	matrix_qr_factorization(&H);

	return 0;
}
