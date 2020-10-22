#include <stdio.h>
#include <sys/time.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"
#include "qpsolver.h"

double time(void)
{
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, NULL);

	if (sec < 0) sec = tv.tv_sec;
	return (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

int main(void)
{
	/* quadratic programming */
	DECLARE_QP_PROBLEM(qp);

	//optimization variable
	DECLARE_VECTOR(x, 2, 1,
                       (0,
                        0));
	PRINT_MATRIX(x);

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
	DECLARE_MATRIX(A_eq, 1, 2,
                       (1, 1));
	DECLARE_MATRIX(b_eq, 2, 1,
                       (0,
                        0));
	PRINT_MATRIX(A_eq);
	PRINT_MATRIX(b_eq);

	//inequality constraints
	DECLARE_MATRIX(lb, 2, 1,
                       (-2,
                        -2));
	DECLARE_MATRIX(ub, 2, 1,
                       (+2,
                        +2));

	DECLARE_MATRIX(A, 3, 2,
                       (+1, +1,
                        -1, +2,
                        +2, +1));
	DECLARE_MATRIX(b, 3, 1,
                       (2,
                        2,
                        3));
	PRINT_MATRIX(A);
	PRINT_MATRIX(b);

	PRINT_MATRIX(lb);
	PRINT_MATRIX(ub);
	PRINT_MATRIX(A);
	PRINT_MATRIX(b);

	qp_solve_set_optimization_variable(&qp, &x);
	qp_solve_set_cost_function(&qp, &P, &q, NULL);
//	qp_solve_set_equality_constraints(&qp, &A_eq, &b_eq);
	qp_solve_set_lower_bound_inequality_constraints(&qp, &lb);
	qp_solve_set_upper_bound_inequality_constraints(&qp, &ub);
	qp_solve_set_affine_inequality_constraints(&qp, &A, &b);

	double start_time = time();
	qp_solve_start(&qp);
	double end_time = time();

	printf("the optimal solution of the problem is:\n");
	PRINT_MATRIX(x);

	printf("run time: %lf seconds\n"
               "optimization took %d iterations\n",
               end_time - start_time, qp.iters);

	return 0;
}
