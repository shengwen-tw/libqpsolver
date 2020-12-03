#ifndef __QPSOLVER_H__
#define __QPSOLVER_H__

#include <stdbool.h>
#include "matrix.h"

#if VERBOSE_MESSAGE == 0
#define VERBOSE_PRINT(...)
#else
#define VERBOSE_PRINT(...) printf(__VA_ARGS__)
#endif

#if DEBUG_MESSAGE == 0
#define DEBUG_PRINT_MATRIX(...)
#define DEBUG_PRINT_VAR(...)
#define DEBUG_PRINT(...)
#else
#define DEBUG_PRINT_MATRIX PRINT_MATRIX
#define DEBUG_PRINT_VAR PRINT_VAR
#define DEBUG_PRINT(...) printf(__VA_ARGS__)
#endif

enum {
	QP_SUCCESS_SOLVED,
	QP_ERROR_NO_OPTIMIZATION_VARIABLE,
	QP_ERROR_NO_OBJECTIVE_FUNCTION,
	QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT,
	QP_ERROR_INCOMPLETE_INEQUAILITY_CONSTRAINT
} QP_RETURN_STATE;

enum {
	QP_PHASE1_FEASIBLE,
	QP_PHASE1_INFEASIBLE
} QP_PHASE1_RETVAL;

typedef struct {
	/* slack variable parameters for feasibility problem */
	double s_margin;
	double beta;

	/* log barrier parameters*/
	double t_init;
	double t_max;
	double mu;

	/* gradient descent */
	double step_size; //not used, replaced with backtracking line search

	/* backtracking line search parameters*/
	double backtracking_alpha;
	double backtracking_beta;

	/* stop criterions */
	double eps;
	int max_iters;
	int iters;
} phase1_param;

typedef struct {
	/* log barrier parameters*/
	double t_init;
	double t_max;
	double mu;

	/* stop criterions */
	double eps;
	int max_iters;
	int iters;
} phase2_param;

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

	/* parameters of phase1 (feasibility) solver */
	phase1_param phase1;

	/* parameters of phase2 (quadratic programming) solver */
	phase2_param phase2;
} qp_t;

void qp_set_default(qp_t *qp);
void qp_config_phase1(qp_t *qp, phase1_param *phase1_config);
void qp_config_phase2(qp_t *qp, phase2_param *phase2_config);
bool qp_start_point_feasibility_check(qp_t *qp);
void qp_solve_set_optimization_variable(qp_t *qp, vector_t *x);
void qp_solve_set_cost_function(qp_t *qp, matrix_t *P, vector_t *q, vector_t *r);
void qp_solve_set_equality_constraints(qp_t *qp, matrix_t *A, vector_t *b);
void qp_solve_set_upper_bound_inequality_constraints(qp_t *qp, vector_t *ub);
void qp_solve_set_lower_bound_inequality_constraints(qp_t *qp, vector_t *lb);
void qp_solve_set_affine_inequality_constraints(qp_t *qp, matrix_t *A, vector_t *b);
int qp_solve_start(qp_t *qp);

#endif
