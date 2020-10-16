#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

int main(void)
{
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

	int pivots[3] = {0};

	solver_linear_system(&A, &X, &B, pivots);
	PRINT_MATRIX(X);

	return 0;
}
