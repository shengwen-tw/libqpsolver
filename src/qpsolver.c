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
	qp->A_eq = NULL;
	qp->b_eq = NULL;
	qp->lb = NULL;
	qp->ub = NULL;

	qp->eps = 1e-6;
	qp->max_iters = 1000000;
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
	qp->A_eq = A;
	qp->b_eq = b;
}

void qp_solve_set_upper_bound_inequality_constraints(qp_t *qp, vector_t *ub)
{
	qp->ub = ub;
}

void qp_solve_set_lower_bound_inequality_constraints(qp_t *qp, vector_t *lb)
{
	qp->lb = lb;
}

void qp_solve_set_affine_inequality_constraints(qp_t *qp, matrix_t *A, vector_t *b)
{
	qp->A = A;
	qp->b = b;
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
	solve_linear_system(qp->P, qp->x, qp->q);
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
	int qb_vec_row = qp->q->row + qp->b_eq->row;
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
		matrix_at(&qb_vec, r, 0) = -matrix_at(qp->q, r, 0);
	}
	//copy b
	for(r = 0; r < qp->b_eq->row; r++) {
		//copy b
		matrix_at(&qb_vec, r + qp->q->row, 0) =
			matrix_at(qp->b_eq, r, 0);
	}
	VERBOSE_PRINT_MATRIX(qb_vec);

	/* construct the KKT matrix */
	int kkt_row = qp->P->row + qp->A_eq->row;
	int kkt_column = qp->P->column + qp->A_eq->row;
	FLOAT *KKT_data = (FLOAT *)calloc(kkt_row * kkt_column, sizeof(FLOAT));
	matrix_t KKT = {
		.data = KKT_data,
		.row = kkt_row,
		.column = kkt_column
	};

	//copy P
	for(r = 0; r < qp->P->row; r++) {
		for(c = 0; c < qp->P->column; c++) {
			matrix_at(&KKT, r, c) = matrix_at(qp->P, r, c);
		}
	}
	//copy A.'
	for(r = 0; r < qp->A_eq->column; r++) {
		for(c = 0; c < qp->A_eq->row; c++) {
			matrix_at(&KKT, r, (c + qp->A_eq->column)) =
				matrix_at(qp->A_eq, c, r);
		}
	}
	//copy A
	for(r = 0; r < qp->A_eq->row; r++) {
		for(c = 0; c < qp->A_eq->column; c++) {
			matrix_at(&KKT, (r + qp->P->row), c) =
				matrix_at(qp->A_eq, r, c);
		}
	}
	//set zero matrix
	for(r = 0; r < qp->A_eq->row; r++) {
		for(c = 0; c < qp->A_eq->row; c++) {
			matrix_at(&KKT, (r + qp->P->row), (c + qp->A_eq->column)) = 0;
		}
	}
	VERBOSE_PRINT_MATRIX(KKT);

	/* construct kkt solution vector */
	int kkt_sol_row = qp->x->row + qp->b_eq->row;
	int kkt_sol_column = 1;
	FLOAT *kkt_sol_data = (FLOAT *)malloc(sizeof(FLOAT) * kkt_sol_row * kkt_sol_column);
	vector_t kkt_sol = {
		.data = kkt_sol_data,
		.row = kkt_sol_row,
		.column = kkt_sol_column
	};

	/* solve the KKT system */
	solve_linear_system(&KKT, &kkt_sol, &qb_vec);
	VERBOSE_PRINT_MATRIX(kkt_sol);

	/* copy the optimal solution back to x */
	for(r = 0; r < qp->x->row; r++) {
		matrix_at(qp->x, r, 0) = matrix_at(&kkt_sol, r, 0);
	}

	/* free mallocs */
	free(KKT.data);
	free(kkt_sol.data);
	free(qb_vec.data);
}

static void qp_solve_inequality_constraint_problem(qp_t *qp)
{
	VERBOSE_PRINT("identify qudratic programming problem with inequality constraint\n");

	//log barrier's parameter
	float t = 2000.0f;

	//save previous optimization result
	matrix_t *x_last = matrix_new(qp->x->row, qp->x->column);
	vector_copy(x_last, qp->x);

	//first derivative of the objective function
	matrix_t *D1_f0 = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the objective function
	matrix_t *D2_f0 = matrix_new(qp->P->row, qp->P->column);
	//inverted second derivative of the objective function
	matrix_t *D2_f0_inv = matrix_new(qp->x->row, qp->x->row);

	//first derivative of the i-th inequality constraint function
	matrix_t *D1_fi = matrix_zeros(qp->x->row, qp->x->column);
	//transposed first derivative of the i-th inequality constraint function
	matrix_t *D1_fi_t = matrix_zeros(qp->x->column, qp->x->row);
	//D[fi(x)] * D[fi(x)].'
	matrix_t *D1_fi_D1_fi_t = matrix_new(qp->x->row, qp->x->row);

	//first derivative of the sumation of log barrier functions
	matrix_t *D1_phi = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the summation of the log barrier functions
	matrix_t *D2_phi = matrix_new(qp->x->row, qp->x->row);

	//newton step's vector
	vector_t *newton_step = vector_new(qp->x->row, qp->x->column);

	/* i-th inenquality constraint function */
	float fi = 0;
	float div_fi = 0;
	float div_fi_squared = 0;

	int r, c;

	while(qp->iters < qp->max_iters) {
		VERBOSE_PRINT("iteration %d\n", qp->iters + 1);

		/* preseve last x for checking convergence */
		vector_copy(x_last, qp->x);	

		t *= 1.001;

		matrix_reset_zeros(D1_phi);
		matrix_reset_zeros(D2_phi);

#if (ENABLE_LOWER_BOUND_INEQUALITY == 1)
		/*===================================================================*
		 * calculate first and second derivative of lower bound inequalities *
		 *===================================================================*/

		//first derivative
		for(r = 0; r < qp->lb->row; r++) {
			fi = -matrix_at(qp->x, r, 0) + matrix_at(qp->lb, r, 0);
			div_fi = 1 / fi;
			div_fi_squared = 1 / (fi * fi);

			matrix_at(D1_phi, r, 0) += div_fi;

			int i, j;
			for(i = 0; i < D1_fi->row; i++) {
				for(j = 0; j < D1_fi->column; j++) {
					matrix_at(D1_fi, i, j) = 0;
				}
			}
			matrix_at(D1_fi, r, 0) = -1;

			//second derivative
			matrix_transpose(D1_fi, D1_fi_t);
			matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
			for(i = 0; i < D2_phi->row; i++) {
				for(j = 0; j < D2_phi->column; j++) {
					matrix_at(D2_phi, i, j) += 
						(div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j));
				}
			}
		}

#endif

#if (ENABLE_UPPER_BOUND_INEQUALITY == 1)
		/*===================================================================*
		 * calculate first and second derivative of upper bound inequalities *
		 *===================================================================*/

		//first derivative
		for(r = 0; r < qp->ub->row; r++) {
			fi = matrix_at(qp->x, r, 0) -  matrix_at(qp->ub, r, 0);
			div_fi = -1 / fi;
			div_fi_squared = 1 / (fi * fi);

			matrix_at(D1_phi, r, 0) += div_fi;

			int i, j;
			for(i = 0; i < D1_fi->row; i++) {
				for(j = 0; j < D1_fi->column; j++) {
					matrix_at(D1_fi, i, j) = 0;
				}
			}
			matrix_at(D1_fi, r, 0) = 1;

			//second derivative
			matrix_transpose(D1_fi, D1_fi_t);
			matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
			for(i = 0; i < D2_phi->row; i++) {
				for(j = 0; j < D2_phi->column; j++) {
					matrix_at(D2_phi, i, j) += 
						(div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j));
				}
			}
		}
#endif

#if (ENABLE_AFFINE_INEQUALITY == 1)
		/*==============================================================*
		 * calculate first and second derivative of affine inequalities *
		 *==============================================================*/
		for(r = 0; r < qp->A->row; r++) {
			int i, j;

			/* calculate constraint function value */
			fi = 0;
			for(j = 0; j < qp->A->column; j++) {
				fi += matrix_at(qp->A, r, j) * matrix_at(qp->x, j, 0);
			}
			fi -= matrix_at(qp->b, r, 0);
			div_fi = -(1 / fi);

			/* calculate first derivative */
			for(i = 0; i < qp->A->row; i++) {
				for(j = 0; j < qp->A->column; j++) {
						matrix_at(D1_fi, i, 0) =
							div_fi * matrix_at(qp->A, r, j);
						matrix_at(D1_fi_t, 0, i) =
							matrix_at(D1_fi, i, 0);
				}
			}
			for(i = 0; i < D1_phi->row; i++) {
				for(j = 0; j < D1_phi->column; j++) {
					matrix_at(D1_phi, i, j) +=
						matrix_at(D1_fi, i, j);
				}
			}

			/* calculate second derivative */
			matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
			div_fi_squared = 1 / (fi * fi);

			for(i = 0; i < D2_phi->row; i++) {
				for(j = 0; j < D2_phi->column; j++) {
					matrix_at(D2_phi, i, j) +=
						div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j);
				}
			}
		}
#endif
		VERBOSE_PRINT_MATRIX(*D1_phi);
		VERBOSE_PRINT_MATRIX(*D2_phi);

		/*============================================================*
		 * 3. calculate the firt derivative of the objective function *
		 *============================================================*/
		matrix_multiply(qp->P, qp->x, D1_f0);
		for(r = 0; r < qp->x->row; r++) {
			matrix_at(D1_f0, r, 0) += matrix_at(qp->q, r, 0);
		}
		vector_scaling(t, D1_f0);

		/*============================================================*
		 * 4. calculate the second derivate of the objective function *
		 *============================================================*/
		matrix_copy(D2_f0, qp->P);
		matrix_scaling(t, D2_f0);

		/*=====================================================================*
		 * combine derivatives of objective function and log barrier functions *
		 *=====================================================================*/

		/* first derivative */
		for(r = 0; r < D1_f0->row; r++) {
			for(c = 0; c < D1_f0->column; c++) {
				matrix_at(D1_f0, r, c) += matrix_at(D1_phi, r, c);
			}
		}

		/* second derivative */
		for(r = 0; r < D2_f0->row; r++) {
			for(c = 0; c < D2_f0->column; c++) {
				matrix_at(D2_f0, r, c) += matrix_at(D2_phi, r, c);
			}
		}

		/*================================*
		 * 6. calculate the newton's step *
		 *================================*/
		matrix_scaling(0.01, D1_f0);

		matrix_inverse(D2_f0, D2_f0_inv);
		matrix_multiply(D2_f0_inv, D1_f0, newton_step);
		vector_negate(newton_step);

		VERBOSE_PRINT_MATRIX(*D1_f0);
		VERBOSE_PRINT_MATRIX(*D2_f0);
		VERBOSE_PRINT_MATRIX(*D2_f0_inv);
		VERBOSE_PRINT_MATRIX(*newton_step);

		/*=====================================*
		 * 7. update the optimization variable *
		 *=====================================*/
		for(r = 0; r < qp->x->row; r++) {
			matrix_at(qp->x, r, 0) += matrix_at(newton_step, r, 0);
		}
		VERBOSE_PRINT_MATRIX(*qp->x);

		qp->iters++;

		FLOAT resid =  vector_residual(qp->x, x_last);
		VERBOSE_PRINT("residual: %f\n", resid);
		VERBOSE_PRINT("---\n");

		/* exit if already converged */
		if(resid < qp->eps) {
			break;
		}
	}

	matrix_delete(D2_f0);
	matrix_delete(D2_f0_inv);
	matrix_delete(D1_f0);
	matrix_delete(D1_phi);
	matrix_delete(D1_fi_t);
	matrix_delete(D1_fi_D1_fi_t);
	matrix_delete(D2_phi);
	vector_delete(newton_step);
}

int qp_solve_start(qp_t *qp)
{
	if(qp->x == NULL) return QP_ERROR_NO_OPTIMIZATION_VARIABLE;
	if(qp->P == NULL) return QP_ERROR_NO_OBJECTIVE_FUNCTION;
	if(qp->A_eq == NULL && qp->b_eq != NULL) return QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT;
	if(qp->A_eq != NULL && qp->b_eq == NULL) return QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT;

	/* no constraint optimization */
	if((qp->A_eq == NULL) && (qp->lb == NULL) && (qp->ub == NULL)) {
		qp_solve_no_constraint_problem(qp);
	}

	/* equality constrained optmization */
	if((qp->A_eq != NULL) && (qp->lb == NULL) && (qp->ub == NULL)) {
		qp_solve_equality_constraint_problem(qp);
	}

	/* inequality constrained optimization */
	if((qp->A_eq == NULL) && ((qp->lb != NULL) || (qp->ub != NULL))) {
		qp_solve_inequality_constraint_problem(qp);
	} 

	/* equality-inequality constrained optimization */
	if((qp->A_eq != NULL) && ((qp->lb != NULL) || (qp->ub != NULL))) {
		qp_solve_all_constraints_problem(qp);
	}

	return QP_SUCCESS_SOLVED;
}
