#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

int main(void)
{
	printf("1.solve linear system AX=B:\n");

	/* solve Ax = b */
	DECLARE_MATRIX(A, 3, 3,
		       (3, 5, 2,
		        2, 1, 3,
		        4, 3, 2));

	DECLARE_MATRIX(B, 3, 2,
		       (57, 23,
		        22, 12,
		        41, 84));

	DECLARE_MATRIX(X, 3, 2,
		       (0, 0,
		        0, 0,
			0, 0));

	int pivots1[3] = {0};

	solve_linear_system(&A, &X, &B, pivots1);

	PRINT_MATRIX(A);
	PRINT_MATRIX(B);
	PRINT_MATRIX(X);

	/* solve matrix inversion */
	printf("\n2.solve matrix inversion inv(M):\n");

	DECLARE_MATRIX(M, 3, 3,
		       ( 1,  0,  2,
		        -1,  5,  0,
		         0,  3, -9));

	DECLARE_MATRIX(M_inv, 3, 3,
		       (0, 0, 0,
		        0, 0, 0,
		        0, 0, 0));

	int pivots2[3] = {0};

	matrix_inverse(&M, &M_inv, pivots2);
	PRINT_MATRIX(M);
	PRINT_MATRIX(M_inv);

	/* solve matrix multiplication */
	printf("\n3.solve matrix multiplication B = AX:\n");
	matrix_multiply(&A, &X, &B);
	PRINT_MATRIX(A);
	PRINT_MATRIX(X);
	PRINT_MATRIX(B);

	return 0;
}
