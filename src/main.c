#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"
#include "qpsolver.h"

int main(void)
{
	/* quadratic programming */
	qp_t qp;

	DECLARE_VECTOR(x, 2, 1,
		       (0,
			0));

	DECLARE_MATRIX(P, 2, 2,
		       (+1, -1,
		        -1, +2));
	DECLARE_VECTOR(q, 2, 1,
		       (-2,
			-6));

	qp_solve_set_optimization_variable(&qp, &x);
	qp_solve_set_cost_function(&qp, &P, &q, NULL);
	qp_solve_start(&qp);

	return 0;
}
