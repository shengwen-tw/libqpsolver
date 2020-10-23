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
	qp_t qp;
	qp_init(&qp);

	//optimization variable
	matrix_t x;
	matrix_construct(&x, 2, 1, ELEMENTS(1.1,
				            1.1));

	//objective function
	matrix_t P, q;
	matrix_construct(&P, 2, 2, ELEMENTS(+1, -1,
                                            -1, +2));
	matrix_construct(&q, 2, 1, ELEMENTS(-2,
                                            -6));

	//equaility constraint
	matrix_t A_eq, b_eq;
	matrix_construct(&A_eq, 1, 2, ELEMENTS(1, 2,
                                               1, 1));
	matrix_construct(&b_eq, 2, 1, ELEMENTS(0,
                                               0));

	//inequality constraints
	vector_t lb, ub;
	vector_construct(&lb, 2, 1, ELEMENTS(-1,
                                             -1));

	vector_construct(&ub, 2, 1, ELEMENTS(3.2,
                                             3.3));

	matrix_t A, b;
	matrix_construct(&A, 3, 2, ELEMENTS(+1, +1,
                                            -1, +2,
                                            +2, +1));
	matrix_construct(&b, 3, 1, ELEMENTS(2,
                                            2,
                                            3));

	PRINT_MATRIX(x);
	PRINT_MATRIX(P);
	PRINT_MATRIX(q);
	PRINT_MATRIX(A_eq);
	PRINT_MATRIX(b_eq);
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
