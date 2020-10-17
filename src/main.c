#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"
#include "qpsolver.h"

int main(void)
{
	/* quadratic programming */
	DECLARE_QP_PROBLEM(qp);

	//optimization variable
	DECLARE_VECTOR(x, 2, 1,
                       (0,
                        0));

	//objective function
	DECLARE_MATRIX(P, 2, 2,
                       (+1, -1,
                        -1, +2));
	DECLARE_VECTOR(q, 2, 1,
                       (-2,
                        -6));
	PRINT_MATRIX(P);
	PRINT_MATRIX(q);

	//equaility constraint
	DECLARE_MATRIX(A, 1, 2,
                       (1, 1));
	DECLARE_MATRIX(b, 2, 1,
                       (0,
                        0));
	PRINT_MATRIX(A);
	PRINT_MATRIX(b);

	qp_solve_set_optimization_variable(&qp, &x);
	qp_solve_set_cost_function(&qp, &P, &q, NULL);
	qp_solve_set_equality_constraints(&qp, &A, &b);
	qp_solve_start(&qp);

	printf("the optimal solution of the problem is:\n");
	PRINT_MATRIX(x);

	return 0;
}
