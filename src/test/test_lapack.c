#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

int main(void)
{
	printf("(1)solve linear system:\n");

	/* solve Ax = b */
	DECLARE_MATRIX(A, 3, 3,
		       (3, 5, 2,
		        2, 1, 3,
		        4, 3, 2));
	PRINT_MATRIX(A);

	DECLARE_MATRIX(B, 3, 2,
		       (57, 23,
		        22, 12,
		        41, 84));
	PRINT_MATRIX(B);

	DECLARE_MATRIX(X, 3, 2,
		       (0, 0,
		        0, 0,
			0, 0));

	int pivots1[3] = {0};

	solver_linear_system(&A, &X, &B, pivots1);
	PRINT_MATRIX(X);

	/* solve matrix inversion */
	printf("\nsolve matrix inversion:\n");

	DECLARE_MATRIX(M, 3, 3,
		       ( 1,  0,  2,
		        -1,  5,  0,
		         0,  3, -9));

	PRINT_MATRIX(M);

	DECLARE_MATRIX(M_inv, 3, 3,
		       (0, 0, 0,
		        0, 0, 0,
		        0, 0, 0));

	int pivots2[3] = {0};

	matrix_inverse(&M, &M_inv, pivots2);
	PRINT_MATRIX(M_inv);

	return 0;
}
