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

bool qp_solve_start(qp_t *qp)
{
	return true;
}
