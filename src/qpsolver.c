#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "qpsolver.h"

void qp_init(qp_t *qp)
{
	qp->x = NULL;
	qp->P = NULL;
	qp->q = NULL;
	qp->r = NULL;
	qp->A = NULL;
	qp->b = NULL;
	qp->lb = NULL;
	qp->ub = NULL;

	qp->eps = 1e-3;
	qp->max_iters = 25;
	qp->iters = 0;
}

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

static void qp_solve_no_constraint_problem(qp_t *qp)
{
	VERBOSE_PRINT("identify quadratic programming problem without any constraints\n");

	/* the closed form solution is given by setting the first derivative equal
	 * to zero, i.e: Px = -q */

	/* construct -q vector */
	int r;
	for(r = 0; r < qp->q->row; r++) {
		qp->q->data[r] *= -1;
	}

	/* solve Px = -q */
	int *pivots = (int *)malloc(sizeof(int) * qp->P->row);
	solve_linear_system(qp->P, qp->x, qp->q, pivots);
	free(pivots);
}

/*static*/ void qp_solve_all_constraints_problem(qp_t *qp)
{
	VERBOSE_PRINT("identify qudratic programming problem with equality and inequality constraints\n");
}

static void qp_solve_equality_constraint_problem(qp_t *qp)
{
	VERBOSE_PRINT("identify qudratic programming problem with equality constraint\n");

	/* the closed form solution of the problem can be obtained by solving the *
	 * KKT system, i.e: [P  A.'][ x*] = [-q]                                  *
	 *                  [A   0 ][nu*]   [ b]                                  *
	 * where x* is the optimal solution and nu* is the optimal dual variable  */

	int r, c;

	/* construct the [-q; b] matrix */
	int qb_vec_row = qp->q->row + qp->b->row;
	int qb_vec_column = 1;
	FLOAT *qb_vec_data = (FLOAT *)malloc(sizeof(FLOAT) * qb_vec_row * qb_vec_column);
	vector_t qb_vec = {
		.data = qb_vec_data,
		.row = qb_vec_row,
		.column = qb_vec_column
	};

	//copy -q
	for(r = 0; r < qp->q->row; r++) {
		//copy -q
		MATRIX_DATA(&qb_vec, r, 0) = -MATRIX_DATA(qp->q, r, 0);
	}
	//copy b
	for(r = 0; r < qp->b->row; r++) {
		//copy b
		MATRIX_DATA(&qb_vec, r + qp->q->row, 0) =
			MATRIX_DATA(qp->b, r, 0);
	}
	VERBOSE_PRINT_MATRIX(qb_vec);

	/* construct the KKT matrix */
	int kkt_row = qp->P->row + qp->A->row;
	int kkt_column = qp->P->column + qp->A->row;
	FLOAT *KKT_data = (FLOAT *)calloc(kkt_row * kkt_column, sizeof(FLOAT));
	matrix_t KKT = {
		.data = KKT_data,
		.row = kkt_row,
		.column = kkt_column
	};

	//copy P
	for(r = 0; r < qp->P->row; r++) {
		for(c = 0; c < qp->P->column; c++) {
			MATRIX_DATA(&KKT, r, c) = MATRIX_DATA(qp->P, r, c);
		}
	}
	//copy A.'
	for(r = 0; r < qp->A->column; r++) {
		for(c = 0; c < qp->A->row; c++) {
			MATRIX_DATA(&KKT, r, (c + qp->A->column)) =
				MATRIX_DATA(qp->A, c, r);
		}
	}
	//copy A
	for(r = 0; r < qp->A->row; r++) {
		for(c = 0; c < qp->A->column; c++) {
			MATRIX_DATA(&KKT, (r + qp->P->row), c) =
				MATRIX_DATA(qp->A, r, c);
		}
	}
	//set zero matrix
	for(r = 0; r < qp->A->row; r++) {
		for(c = 0; c < qp->A->row; c++) {
			MATRIX_DATA(&KKT, (r + qp->P->row), (c + qp->A->column)) = 0;
		}
	}
	VERBOSE_PRINT_MATRIX(KKT);

	/* construct kkt solution vector */
	int kkt_sol_row = qp->x->row + qp->b->row;
	int kkt_sol_column = 1;
	FLOAT *kkt_sol_data = (FLOAT *)malloc(sizeof(FLOAT) * kkt_sol_row * kkt_sol_column);
	vector_t kkt_sol = {
		.data = kkt_sol_data,
		.row = kkt_sol_row,
		.column = kkt_sol_column
	};

	/* solve the KKT system */
	int *pivots = (int *)malloc(sizeof(int) * (qp->P->row + qp->A->row));
	solve_linear_system(&KKT, &kkt_sol, &qb_vec, pivots);
	VERBOSE_PRINT_MATRIX(kkt_sol);

	/* copy the optimal solution back to x */
	for(r = 0; r < qp->x->row; r++) {
		MATRIX_DATA(qp->x, r, 0) = MATRIX_DATA(&kkt_sol, r, 0);
	}

	/* free mallocs */
	free(KKT.data);
	free(kkt_sol.data);
	free(qb_vec.data);
}

static void qp_solve_inequality_constraint_problem(qp_t *qp)
{
	VERBOSE_PRINT("identify qudratic programming problem with inequality constraint\n");

	int r;
	int *pivots = (int *)malloc(sizeof(int) * qp->P->row);

	/* log barrier's parameter */
	float t = 100.0f;

	/* save previous optimization result */
	FLOAT *x_last_data = (FLOAT *)malloc(sizeof(FLOAT) * qp->x->row * qp->x->column);
	vector_t x_last = {
		.data = x_last_data,
		.row = qp->x->row,
		.column = qp->x->column
	};
	vector_copy(&x_last, qp->x);

	/* calculate the inverted second derivative of the objective function */
	FLOAT *D2_inv_data = (FLOAT *)malloc(sizeof(FLOAT) * qp->P->row * qp->P->column);
	matrix_t D2_inv = {
		.data = D2_inv_data,
		.row = qp->P->row,
		.column = qp->P->column
	};
	matrix_inverse(qp->P, &D2_inv, pivots);
	
	FLOAT *D1_data = (FLOAT *)malloc(sizeof(FLOAT) * qp->x->row * qp->x->column);
	matrix_t D1 = {
		.data = D1_data,
		.row = qp->x->row,
		.column = qp->x->column
	};

	FLOAT *newton_step_data =
		(FLOAT *)malloc(sizeof(FLOAT) * qp->x->row * qp->x->column);
	vector_t newton_step = {
		.data = newton_step_data,
		.row = qp->x->row,
		.column = qp->x->column
	};

	/* log barrier function of inequality constraints */
	int log_barriers_row = qp->lb->row + qp->ub->row;
	int log_barriers_column = 1;
	FLOAT *log_barriers_data =
		(FLOAT *)malloc(sizeof(FLOAT) * log_barriers_row * log_barriers_column);
	vector_t log_barriers = {
		.data = log_barriers_data,
		.row = log_barriers_row,
		.column = log_barriers_column,
	};

	int D1_phi_row = qp->lb->row + qp->ub->row;
	int D1_phi_column = qp->x->row;
	FLOAT *D1_phi_data = (FLOAT *)malloc(sizeof(FLOAT) * D1_phi_row * D1_phi_column);
	matrix_t D1_phi = {
		.data = D1_phi_data,
		.row = D1_phi_row,
		.column = D1_phi_column
	};

	/* inequality constraints with the form of Fx < c */
	int F_row = qp->lb->row + qp->ub->row;
	int F_column = qp->x->row;
	FLOAT *F_data = (FLOAT *)malloc(sizeof(FLOAT) * F_row * F_column);
	matrix_t F = {
		.data = F_data,
		.row = F_row,
		.column = F_column
	};
	//lower bound
	for(r = 0; r < qp->lb->row; r++) {
		MATRIX_DATA(&F, r, r) = 1;
	}
	//upper bound
	for(r = 0; r < qp->ub->row; r++) {
		int _r = r + qp->lb->row;
		MATRIX_DATA(&F, _r, r) = 1;
	}
	//affine inequality
	//TODO
	VERBOSE_PRINT_MATRIX(F);

	/* F*trans(T), we use this to calculate the second derivative of the log barrier function */
	FLOAT *Ft_data = (FLOAT *)malloc(sizeof(FLOAT) * F_row * F_column);
	matrix_t Ft = {.data = Ft_data};
	matrix_transpose(&F, &Ft);
	VERBOSE_PRINT_MATRIX(Ft);

	FLOAT *FFt_data = (FLOAT *)malloc(sizeof(FLOAT) * F_row * F_row);
	matrix_t FFt = {
		.data = FFt_data,
		.row = F_row,
		.column = F_row
	};
	matrix_multiply(&F, &Ft, &FFt);
	VERBOSE_PRINT_MATRIX(FFt);

	FLOAT *D2_phi_data = (FLOAT *)malloc(sizeof(FLOAT) * F_row * F_row);
	matrix_t D2_phi = {
		.data = D2_phi_data,
		.row = F_row,
		.column = F_row
	};

	FLOAT *D2_phi_inv_data = (FLOAT *)malloc(sizeof(FLOAT) * F_row * F_row);
	matrix_t D2_phi_inv = {
		.data = D2_phi_inv_data,
		.row = F_row,
		.column = F_row
	};

	while(qp->iters < qp->max_iters) {
		VERBOSE_PRINT("iteration %d\n", qp->iters + 1);

		//preseve for checking convergence
		vector_copy(&x_last, qp->x);	

		/* calculate the log barriers */
		//lower bound
		for(r = 0; r < qp->lb->row; r++) {
			float log_barrier_i =
				-log10f(MATRIX_DATA(qp->x, r, 0) - MATRIX_DATA(qp->lb, r, 0));

			MATRIX_DATA(&log_barriers, r, 0) = log_barrier_i;
		}
		//upper bound
		for(r = 0; r < qp->ub->row; r++) {
			int _r = (r + qp->lb->row);

			float log_barrier_i =
				-log10f(-MATRIX_DATA(qp->x, r, 0) + MATRIX_DATA(qp->ub, r, 0));

			MATRIX_DATA(&log_barriers, _r, 0) = log_barrier_i;
		}
		//affine inequality constraints
		//TODO
		VERBOSE_PRINT_MATRIX(log_barriers);

		/* calculate the first derivative of the log barrier function */
		for(r = 0; r < qp->lb->row; r++) {
			MATRIX_DATA(&D1_phi, r, r) = (1 / MATRIX_DATA(&log_barriers, r, 0));
		}
		for(r = 0; r < qp->ub->row; r++) {
			int _r = r + qp->lb->row;
			MATRIX_DATA(&D1_phi, _r, r) = (1 / MATRIX_DATA(&log_barriers, _r, 0));
		}
		VERBOSE_PRINT_MATRIX(D1_phi);

		/* calculate the inverted second derivative of the log barrier function */
		float div_sum_of_squared_log_barriers = 0;
		for(r = 0; r < log_barriers.row; r++) {
			float squared_log_barrier =
				MATRIX_DATA(&log_barriers, r, 0) * MATRIX_DATA(&log_barriers, r, 0);

			div_sum_of_squared_log_barriers += squared_log_barrier;
		}
		div_sum_of_squared_log_barriers = 1 / div_sum_of_squared_log_barriers;

		matrix_copy(&D2_phi, &FFt);
		matrix_scaling(div_sum_of_squared_log_barriers, &D2_phi);
		VERBOSE_PRINT_MATRIX(D2_phi);

		int *pivots = (int *)malloc(sizeof(int) * D2_phi_inv.row);
		matrix_inverse(&D2_phi, &D2_phi_inv, pivots);
		free(pivots);
		VERBOSE_PRINT_MATRIX(D2_phi_inv);

		/* calculate the firt derivative of the objective function *
		 * D[f(x)] = Px + r                                        */
		matrix_multiply(qp->P, qp->x, &D1);
		//TODO: implment vector addition
		for(r = 0; r < qp->x->row; r++) {
			MATRIX_DATA(&D1, r, 0) += MATRIX_DATA(qp->q, r, 0);
		}

		/* effected by the new objective function t*f(x) + sum(i=1,m){phi(x)} */
		vector_scaling(t, &D1);
		matrix_scaling(1.0 / t, &D2_inv);

		/* calculate the  newton's step          * 
		 * newton_step = -D^2[f(x)]^-1 * D[f(x)] */
		matrix_multiply(&D2_inv, &D1, &newton_step);
		vector_negate(&newton_step);
		VERBOSE_PRINT_MATRIX(D2_inv);
		//VERBOSE_PRINT_MATRIX(D1);
		VERBOSE_PRINT_MATRIX(newton_step);

		/* update the optimization variable *
		 * x(k+1) = x(k) + newton_step      */
		for(r = 0; r < qp->x->row; r++) {
			MATRIX_DATA(qp->x, r, 0) += MATRIX_DATA(&newton_step, r, 0);
		}
		VERBOSE_PRINT_MATRIX(*qp->x);

		qp->iters++;

		FLOAT resid =  vector_residual(qp->x, &x_last);
		VERBOSE_PRINT("residual: %f\n", resid);
		VERBOSE_PRINT("---\n");

		/* exit if already converged */
		if(resid < qp->eps) {
			break;
		}
	}

	free(D2_inv.data);
	free(D1.data);
	free(D1_phi.data);
	free(F.data);
	free(Ft.data);
	free(newton_step.data);
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
