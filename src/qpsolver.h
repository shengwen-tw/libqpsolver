#ifndef __QPSOLVER_H__
#define __QPSOLVER_H__

#include <stdbool.h>
#include "matrix.h"

#if VERBOSE == 0
#define VERBOSE_PRINT(...)
#define VERBOSE_PRINT_MATRIX(...)
#else
#define VERBOSE_PRINT(...) printf(__VA_ARGS__)
#define VERBOSE_PRINT_MATRIX PRINT_MATRIX
#endif

enum {
	QP_SUCCESS_SOLVED,
	QP_ERROR_NO_OPTIMIZATION_VARIABLE,
	QP_ERROR_NO_OBJECTIVE_FUNCTION,
	QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT
} QP_RETURN_STATE;

typedef struct {
	/* optimization variable */
	vector_t *x;

	/* cost function */
	matrix_t *P;
	vector_t *q;
	vector_t *r;

	/* equality constraint */
	matrix_t *A_eq;
	vector_t *b_eq;

	/* inequality constraints */
	vector_t *lb; //lower bound inequality
	vector_t *ub; //upper bound inequality
	matrix_t *A;  //affine inequality matrix
	vector_t *b;  //affine inequality vector

	/* stop criterions */
	FLOAT eps;
	int max_iters;
	int iters;
} qp_t;

void qp_init(qp_t *qp);
void qp_solve_set_optimization_variable(qp_t *qp, vector_t *x);
void qp_solve_set_cost_function(qp_t *qp, matrix_t *P, vector_t *q, vector_t *r);
void qp_solve_set_equality_constraints(qp_t *qp, matrix_t *A, vector_t *b);
void qp_solve_set_upper_bound_inequality_constraints(qp_t *qp, vector_t *ub);
void qp_solve_set_lower_bound_inequality_constraints(qp_t *qp, vector_t *lb);
void qp_solve_set_affine_inequality_constraints(qp_t *qp, matrix_t *A, vector_t *b);
int qp_solve_start(qp_t *qp);

#endif
