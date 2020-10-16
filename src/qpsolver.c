#include <stdbool.h>
#include "qpsolver.h"

void qp_solve_set_optimization_variable(qp_t *qp, vector_t *x)
{
	qp->x = x;
}

void qp_solve_set_cost_function(qp_t *qp, matrix_t *P, vector_t *q, vector_t *r)
{
	qp->P = P;
	qp->q = q;
	qp->r = r;
}

void qp_solve_set_equality_constraints(qp_t *qp, matrix_t *A, vector_t *b)
{
	qp->A = A;
	qp->b = b;
}

void qp_solve_set_upper_bound_inequality_constraints(qp_t *qp, vector_t *ub)
{
	qp->ub = ub;
}

void qp_solve_set_lower_bound_inequality_constraints(qp_t *qp, vector_t *lb)
{
	qp->lb = lb;
}

/*static*/ void qp_solve_no_constraint_problem(qp_t *qp)
{
	/* x = H \ (-f) */
}

/*static*/ void qp_solve_all_constraints_problem(qp_t *qp)
{
}

/*static*/ void qp_solve_equality_constraint_problem(qp_t *qp)
{
}

/*static*/ void qp_solve_inequality_constraint_problem(qp_t *qp)
{
}

int qp_solve_start(qp_t *qp)
{
	if(qp->x == NULL) return QP_ERROR_NO_OPTIMIZATION_VARIABLE;
	if(qp->P == NULL) return QP_ERROR_NO_OBJECTIVE_FUNCTION;
	if(qp->A == NULL && qp->b != NULL) return QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT;
	if(qp->A != NULL && qp->b == NULL) return QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT;

	/* no constraint optimization */
	if((qp->A == NULL) && (qp->lb == NULL) && (qp->ub == NULL)) {
		qp_solve_no_constraint_problem(qp);
	}

	/* equality constrained optmization */
	if((qp->A != NULL) && (qp->lb == NULL) && (qp->ub == NULL)) {
		qp_solve_equality_constraint_problem(qp);
	}

	/* inequality constrained optimization */
	if((qp->A == NULL) && ((qp->lb != NULL) || (qp->ub != NULL))) {
		qp_solve_inequality_constraint_problem(qp);
	} 

	/* equality-inequality constrained optimization */
	if((qp->A != NULL) && ((qp->lb != NULL) || (qp->ub != NULL))) {
		qp_solve_all_constraints_problem(qp);
	}

	return QP_SUCCESS_SOLVED;
}
