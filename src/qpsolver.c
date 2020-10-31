#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "qpsolver.h"

void qp_set_default(qp_t *qp)
{
	qp->x = NULL;
	qp->P = NULL;
	qp->q = NULL;
	qp->r = NULL;
	qp->A_eq = NULL;
	qp->b_eq = NULL;
	qp->lb = NULL;
	qp->ub = NULL;
	qp->A = NULL;
	qp->b = NULL;

	qp->iters = 0;

	qp->eps = 1e-6;        //residual value to stop the optimization
	qp->a = 0.95;          //line searching parameter of the newton step
	qp->mu = 20;           //stiffness growth rate of the log barrier function
	qp->t_init = 1;        //initial stiffness of the log barrier
	qp->t_max_inc = 10000;       //maximum stiffness of the log barrier
	qp->max_iters = 1;  //maximum iteration times

	qp->line_search_num = 100;
	qp->line_search_min_step_size = 0.1;
}

bool qp_start_point_feasibility_check(qp_t *qp)
{
	float fi;
	int r, c;

	/* equality constraints */
	if((qp->A_eq != NULL) && (qp->b_eq != NULL)) {
		bool infeasible = false;
		const float tolerance = 1e-3;

		vector_t *Ax = matrix_new(qp->b_eq->row, qp->b_eq->column);
		vector_t *lin_sys_sol = matrix_new(qp->b_eq->row, qp->b_eq->column);

		matrix_multiply(qp->A_eq, qp->x, Ax);
		matrix_sub(Ax, qp->b_eq, lin_sys_sol);

		for(r = 0; r < lin_sys_sol->row; r++) {
			if(fabs(matrix_at(lin_sys_sol, r, 0)) > tolerance) {
				infeasible = true;
			}
		}

		matrix_delete(Ax);
		matrix_delete(lin_sys_sol);

		if(infeasible == true) {
			return false;
		}
	}

	/* lower bound inequality constraints */
	if(qp->lb != NULL) {
		for(r = 0; r < qp->lb->row; r++) {
			if(matrix_at(qp->x, r, 0) < matrix_at(qp->lb, r, 0)) {
				return false;
			}
		}
	}

	/* upper bound inequality constraints */
	if(qp->ub != NULL) {
		for(r = 0; r < qp->ub->row; r++) {
			if(matrix_at(qp->x, r, 0) > matrix_at(qp->ub, r, 0)) {
				return false;
			}
		}
	}

	/* affine inequality constraints */
	if((qp->A != NULL) && (qp->b != NULL)) {
		for(r = 0; r < qp->A->row; r++) {
			fi = 0;
			for(c = 0; c < qp->A->column; c++) {
				fi += matrix_at(qp->A, r, c) * matrix_at(qp->x, c, 0);
			}

			if(fi > matrix_at(qp->b, r, 0)) {
				return false;
			}
		}
	}

	return true;
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

static void qp_solve_equality_constraint_problem(qp_t *qp)
{
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

static void qp_solve_inequality_constraint_problem(qp_t *qp, bool solve_lower_bound,
        bool solve_upper_bound, bool solve_affine_inequality)
{
	const FLOAT epsilon = 1e-14; //increase numerical stability of divide by zero

	//log barrier's parameter
	float t = qp->t_init;

	//save previous optimization result
	matrix_t *x_last = matrix_new(qp->x->row, qp->x->column);
	matrix_copy(x_last, qp->x);

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
	vector_t *newton_step = matrix_new(qp->x->row, qp->x->column);

	/* i-th inenquality constraint function */
	float fi = 0;
	float div_fi = 0;
	float div_fi_squared = 0;

	int r;
	int i, j;

	/* outer loop varies the stiffness of the log barrier functions */
	while(1) {
		if(t > qp->t_max_inc) {
			break;
		}

		/* inner loop do the gradient descent */
		while(qp->iters < qp->max_iters) {
			VERBOSE_PRINT("iteration %d\n", qp->iters + 1);

			/* preseve last x for checking convergence */
			matrix_copy(x_last, qp->x);

			matrix_reset_zeros(D1_phi);
			matrix_reset_zeros(D2_phi);

#if (ENABLE_LOWER_BOUND_INEQUALITY == 1)
			/*===================================================================*
			 * calculate first and second derivative of lower bound inequalities *
			 *===================================================================*/

			if(solve_lower_bound == true) {
				for(r = 0; r < qp->lb->row; r++) {
					fi = -matrix_at(qp->x, r, 0) + matrix_at(qp->lb, r, 0);
					div_fi = 1 / (fi + epsilon);
					div_fi_squared = 1 / ((fi * fi) + epsilon);

					//first derivative of the inequality constraint function
					matrix_reset_zeros(D1_fi);
					matrix_at(D1_fi, r, 0) = -1;

					//first derivative of the log barrier function
					matrix_at(D1_phi, r, 0) += div_fi;

					//second derivative of the log barrier function
					matrix_transpose(D1_fi, D1_fi_t);
					matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
					for(i = 0; i < D2_phi->row; i++) {
						for(j = 0; j < D2_phi->column; j++) {
							matrix_at(D2_phi, i, j) +=
							    (div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j));
						}
					}
				}
			}
#endif

#if (ENABLE_UPPER_BOUND_INEQUALITY == 1)
			/*===================================================================*
			 * calculate first and second derivative of upper bound inequalities *
			 *===================================================================*/

			if(solve_upper_bound == true) {
				for(r = 0; r < qp->ub->row; r++) {
					fi = matrix_at(qp->x, r, 0) -  matrix_at(qp->ub, r, 0);
					div_fi = -1 / (fi + epsilon);
					div_fi_squared = 1 / ((fi * fi) + epsilon);

					//first derivative of the inequality constraint function
					matrix_reset_zeros(D1_fi);
					matrix_at(D1_fi, r, 0) = 1;

					//first derivative of the log barrier function
					matrix_at(D1_phi, r, 0) += div_fi;

					//second derivative of the log barrier function
					matrix_transpose(D1_fi, D1_fi_t);
					matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
					for(i = 0; i < D2_phi->row; i++) {
						for(j = 0; j < D2_phi->column; j++) {
							matrix_at(D2_phi, i, j) +=
							    (div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j));
						}
					}
				}
			}
#endif

#if (ENABLE_AFFINE_INEQUALITY == 1)
			/*==============================================================*
			 * calculate first and second derivative of affine inequalities *
			 *==============================================================*/

			if(solve_affine_inequality == true) {
				for(r = 0; r < qp->A->row; r++) {
					/* calculate constraint function value */
					fi = 0;
					for(j = 0; j < qp->A->column; j++) {
						fi += matrix_at(qp->A, r, j) * matrix_at(qp->x, j, 0);
					}
					fi -= matrix_at(qp->b, r, 0);
					div_fi = -1 / (fi + epsilon);
					div_fi_squared = 1 / ((fi * fi) + epsilon);

					/* calculate first derivative */
					for(i = 0; i < D1_phi->row; i++) {
						matrix_at(D1_phi, i, 0) +=
						    div_fi * matrix_at(qp->A, r, i);

						matrix_at(D1_fi, i, 0) = matrix_at(qp->A, r, i);
					}

					/* calculate second derivative */
					matrix_transpose(D1_fi, D1_fi_t);
					matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);

					for(i = 0; i < D2_phi->row; i++) {
						for(j = 0; j < D2_phi->column; j++) {
							matrix_at(D2_phi, i, j) +=
							    div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j);
						}
					}
				}
			}
#endif

			/*=========================================================*
			 * calculate the firt derivative of the objective function *
			 *=========================================================*/

			matrix_multiply(qp->P, qp->x, D1_f0);
			matrix_add_by(D1_f0, qp->q);
			matrix_scale_by(t, D1_f0);

			/*=========================================================*
			 * calculate the second derivate of the objective function *
			 *=========================================================*/

			matrix_copy(D2_f0, qp->P);
			matrix_scale_by(t, D2_f0);

			/*=====================================================================*
			 * combine derivatives of objective function and log barrier functions *
			 *=====================================================================*/

			matrix_add_by(D1_f0, D1_phi);
			matrix_add_by(D2_f0, D2_phi);

			/*===================================*
			 * gradient descent with newton step *
			 *===================================*/

			//calculate the newton step
			matrix_inverse(D2_f0, D2_f0_inv);
			matrix_multiply(D2_f0_inv, D1_f0, newton_step);
			matrix_scale_by(-1, newton_step);

			bool step_too_large;
			while(1) {
				step_too_large = false;

				/* update the optimization variable */
				for(r = 0; r < qp->x->row; r++) {
					matrix_add(x_last, newton_step, qp->x);

#if (ENABLE_LOWER_BOUND_INEQUALITY == 1)
					/*==================================================*
					 * check if lower bound constraints are still valid *
					 *==================================================*/
					if(solve_lower_bound == true) {
						if(matrix_at(qp->x, r, 0) < matrix_at(qp->lb, r, 0)) {
							step_too_large = true;
							break;
						}
					}
#endif

#if (ENABLE_UPPER_BOUND_INEQUALITY == 1)
					/*==================================================*
					 * check if upper bound constraints are still valid *
					 *==================================================*/
					if(solve_upper_bound == true) {
						if(matrix_at(qp->x, r, 0) > matrix_at(qp->ub, r, 0)) {
							step_too_large = true;
							break;
						}
					}
#endif
				}

#if (ENABLE_AFFINE_INEQUALITY == 1)
				/*========================================================*
				 * check if affine inequality constraints are still valid *
				 *========================================================*/
				if(solve_affine_inequality == true) {
					for(i = 0; i < qp->A->row; i++) {
						fi = 0;
						for(j = 0; j < qp->A->column; j++) {
							fi += matrix_at(qp->A, i, j) * matrix_at(qp->x, j, 0);
						}

						if(fi > matrix_at(qp->b, i, 0)) {
							step_too_large = true;
							break;
						}
					}
				}
#endif

				/*================================================================*
				 * schrink the newton step if the step is too large and break the *
				 * inequality constraints                                         *
				 *================================================================*/
				if(step_too_large == false) {
					break;
				} else {
					matrix_scale_by(qp->a, newton_step);
				}
			}

			VERBOSE_PRINT_MATRIX(*D1_phi);
			VERBOSE_PRINT_MATRIX(*D2_phi);
			VERBOSE_PRINT_MATRIX(*D1_f0);
			VERBOSE_PRINT_MATRIX(*D2_f0);
			VERBOSE_PRINT_MATRIX(*D2_f0_inv);
			VERBOSE_PRINT_MATRIX(*newton_step);
			VERBOSE_PRINT_MATRIX(*qp->x);
			VERBOSE_PRINT("t = %f\n", t);

			qp->iters++;

			FLOAT resid =  vector_residual(qp->x, x_last);
			VERBOSE_PRINT("residual: %f\n", resid);
			VERBOSE_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->eps) {
				break;
			}
		}

		t += qp->mu;
	}

	matrix_delete(x_last);
	matrix_delete(D1_f0);
	matrix_delete(D2_f0);
	matrix_delete(D2_f0_inv);
	matrix_delete(D1_fi);
	matrix_delete(D1_fi_t);
	matrix_delete(D1_fi_D1_fi_t);
	matrix_delete(D1_phi);
	matrix_delete(D2_phi);
	matrix_delete(newton_step);
}

static void qp_solve_equality_inequality_constraint_problem(qp_t *qp, bool solve_lower_bound,
        bool solve_upper_bound, bool solve_affine_inequality)
{
	const FLOAT epsilon = 1e-14; //increase numerical stability of divide by zero

	int r, c;
	int i, j;

	//log barrier's parameter
	float t = qp->t_init;

	//first derivative of the objective function
	matrix_t *D1_f0 = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the objective function
	matrix_t *D2_f0 = matrix_new(qp->P->row, qp->P->column);

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
	vector_t *newton_step = matrix_new(qp->x->row, qp->x->column);

	//newton step times step size
	vector_t *scaled_newton_step = matrix_new(qp->x->row, qp->x->column);

	/* i-th inenquality constraint function */
	float fi = 0;
	float div_fi = 0;
	float div_fi_squared = 0;

	/*========================================================*
	 * eliminate equality constraints by null space transform *
	 *========================================================*/

	matrix_t *F = matrix_zeros(qp->A_eq->column, qp->A_eq->column);
	matrix_t *F_t = matrix_new(F->column, F->column);
	matrix_t *Q, *R;

	/*==============================================================================*
	 * solve the problem with new optimization variable z in the null space of A_eq *
	 *==============================================================================*/

	matrix_t *x_hat = matrix_zeros(qp->x->row, qp->x->column);

	//new optimization variable z
	matrix_t *z_now = matrix_zeros(qp->x->row, qp->x->column);
	matrix_t *z_last = matrix_zeros(qp->x->row, qp->x->column);

	//x = Fx + x_hat
	matrix_t *x_now = matrix_zeros(qp->x->row, qp->x->column);
	matrix_t *x_last = matrix_zeros(qp->x->row, qp->x->column);

	//initialize x with z
	matrix_multiply(F, z_now, x_now);
	matrix_add_by(x_now, x_hat);

	//first derivative of the equality constraints eliminated objective function
	matrix_t *D1_f_tilde = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the equality constraints eliminated objective function
	matrix_t *D2_f_tilde = matrix_new(qp->P->row, qp->P->column);
	//inverted second derivative of the equality constraints eliminated objective function
	matrix_t *D2_f_tilde_inv = matrix_new(qp->x->row, qp->x->row);
	//D2_f0 times F
	matrix_t *D2_f0_F = matrix_new(D2_f_tilde->row, D2_f_tilde->column);

	/*=====================================*
	 * caclulate null space matrix of A_eq *
	 *=====================================*/

	/* if user input a rectangular A_eq matrix, fill up zeros to make it
	 * become a square matrix (required by QR factorization in Lapack) */
	matrix_t *A_eq_square = matrix_zeros(qp->x->row, qp->x->row);
	for(r = 0; r < qp->A_eq->row; r++) {
		for(c = 0; c < qp->A_eq->column; c++) {
			/* combine copy and tranpose in one line,
			 * transpose is for calculating the null space of A_eq */
			matrix_at(A_eq_square, c, r) = matrix_at(qp->A_eq, r, c);
		}
	}
	matrix_qr_factorization(A_eq_square, &Q, &R);
	PRINT_MATRIX(*A_eq_square);

	/* calculate the zeros row count of R matrix */
	int n_zero_cols = 0;
	for(r = (R->row - 1); r >= 0; r--) {
		FLOAT norm = 0;
		for(c = 0; c < R->column; c++) {
			norm += matrix_at(R, r, c) * matrix_at(R, r, c);
		}

		if(fabs(norm) < 1e-8) {
			n_zero_cols++;
		}
	}

	/* copy null space bias vectors from Q matrix */
	for(r = 0; r < F->row; r++) {
		for(c = 0; c < n_zero_cols; c++) {
			matrix_at(F, r, c) = matrix_at(Q, r, (Q->column - n_zero_cols + c));
		}
	}

	matrix_transpose(F, F_t);

	/*====================*
	 * optimization start *
	 *====================*/

	/* outer loop varies the stiffness of the log barrier functions */
	while(1) {
		if(t > (qp->t_init + qp->t_max_inc)) {
			break;
		}

		/* inner loop do the gradient descent */
		while(qp->iters < qp->max_iters) {
			VERBOSE_PRINT("iteration %d\n", qp->iters + 1);

			/* preseve last x and z variable */
			matrix_copy(z_last, z_now);
			matrix_copy(x_last, x_now);

			matrix_reset_zeros(D1_phi);
			matrix_reset_zeros(D2_phi);

#if (ENABLE_LOWER_BOUND_INEQUALITY == 1)
			/*===================================================================*
			 * calculate first and second derivative of lower bound inequalities *
			 *===================================================================*/

			if(solve_lower_bound == true) {
				for(r = 0; r < qp->lb->row; r++) {
					fi = -matrix_at(x_now, r, 0) + matrix_at(qp->lb, r, 0);
					div_fi = 1 / (fi + epsilon);
					div_fi_squared = 1 / ((fi * fi) + epsilon);

					//first derivative of the inequality constraint function
					matrix_reset_zeros(D1_fi);
					matrix_at(D1_fi, r, 0) = -1;

					//first derivative of the log barrier function
					matrix_at(D1_phi, r, 0) += div_fi;

					//second derivative of the log barrier function
					matrix_transpose(D1_fi, D1_fi_t);
					matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
					for(i = 0; i < D2_phi->row; i++) {
						for(j = 0; j < D2_phi->column; j++) {
							matrix_at(D2_phi, i, j) +=
							    (div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j));
						}
					}
				}
			}
#endif

#if (ENABLE_UPPER_BOUND_INEQUALITY == 1)
			/*===================================================================*
			 * calculate first and second derivative of upper bound inequalities *
			 *===================================================================*/

			if(solve_upper_bound == true) {
				for(r = 0; r < qp->ub->row; r++) {
					fi = matrix_at(x_now, r, 0) -  matrix_at(qp->ub, r, 0);
					div_fi = -1 / (fi + epsilon);
					div_fi_squared = 1 / ((fi * fi) + epsilon);

					//first derivative of the inequality constraint function
					matrix_reset_zeros(D1_fi);
					matrix_at(D1_fi, r, 0) = 1;

					//first derivative of the log barrier function
					matrix_at(D1_phi, r, 0) += div_fi;

					//second derivative of the log barrier function
					matrix_transpose(D1_fi, D1_fi_t);
					matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);
					for(i = 0; i < D2_phi->row; i++) {
						for(j = 0; j < D2_phi->column; j++) {
							matrix_at(D2_phi, i, j) +=
							    (div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j));
						}
					}
				}
			}
#endif

#if (ENABLE_AFFINE_INEQUALITY == 1)
			/*==============================================================*
			 * calculate first and second derivative of affine inequalities *
			 *==============================================================*/

			if(solve_affine_inequality == true) {
				for(r = 0; r < qp->A->row; r++) {
					/* calculate constraint function value */
					fi = 0;
					for(j = 0; j < qp->A->column; j++) {
						fi += matrix_at(qp->A, r, j) * matrix_at(x_now, j, 0);
					}
					fi -= matrix_at(qp->b, r, 0);
					div_fi = -1 / (fi + epsilon);
					div_fi_squared = 1 / ((fi * fi) + epsilon);

					/* calculate first derivative */
					for(i = 0; i < D1_phi->row; i++) {
						matrix_at(D1_phi, i, 0) +=
						    div_fi * matrix_at(qp->A, r, i);

						matrix_at(D1_fi, i, 0) = matrix_at(qp->A, r, i);
					}

					/* calculate second derivative */
					matrix_transpose(D1_fi, D1_fi_t);
					matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);

					for(i = 0; i < D2_phi->row; i++) {
						for(j = 0; j < D2_phi->column; j++) {
							matrix_at(D2_phi, i, j) +=
							    div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j);
						}
					}
				}
			}
#endif

			/*=========================================================*
			 * calculate the firt derivative of the objective function *
			 *=========================================================*/

			matrix_multiply(qp->P, x_now, D1_f0);
			matrix_add_by(D1_f0, qp->q);
			matrix_scale_by(t, D1_f0);

			/*=========================================================*
			 * calculate the second derivate of the objective function *
			 *=========================================================*/

			matrix_copy(D2_f0, qp->P);
			matrix_scale_by(t, D2_f0);

			/*=====================================================================*
			 * combine derivatives of objective function and log barrier functions *
			 *=====================================================================*/

			matrix_add_by(D1_f0, D1_phi);
			matrix_add_by(D2_f0, D2_phi);

			/*=====================================================================*
			 * multiply first and second derivatives of objective function bu null *
			 * matrix F                                                            *
			 *=====================================================================*/

			//first derivative
			matrix_multiply(F_t, D1_f0, D1_f_tilde);

			//second derivative
			matrix_multiply(D2_f0, F, D2_f0_F);
			matrix_multiply(F_t, D2_f0_F, D2_f_tilde);

			/*===================================*
			 * gradient descent with newton step *
			 *===================================*/

			//calculate the newton step
			matrix_inverse(D2_f_tilde, D2_f_tilde_inv);
			matrix_multiply(D2_f_tilde_inv, D1_f_tilde, newton_step);
			matrix_scale_by(-1, newton_step);

			/* exact line search */
			float curr_step_size;
			float best_step_size = qp->line_search_num;
			float curr_norm = 0;
			float smallest_norm = 0;

			/* initial the step size */
			matrix_add(z_last, newton_step, z_now); //scale = 1
			matrix_multiply(F, z_now, x_now);
			matrix_add_by(x_now, x_hat);
			matrix_multiply(qp->P, x_now, D1_f0);
			matrix_add_by(D1_f0, qp->q);
			matrix_scale_by(t, D1_f0);
			for(r = 0; r < D1_f0->row; r++) {
				smallest_norm += matrix_at(D1_f0, r, 0) * matrix_at(D1_f0, r, 0);
			}

			/* find the best step size iteratively */
			for(i = (qp->line_search_num - 1); i > 0; i--) {
				curr_step_size = (float)i / (float)qp->line_search_num;

				if(curr_step_size <= qp->line_search_min_step_size) {
					break;
				}

				matrix_scaling(curr_step_size, newton_step, scaled_newton_step);
				matrix_add(z_last, scaled_newton_step, z_now);

				//calculate x from z
				matrix_multiply(F, z_now, x_now);
				matrix_add_by(x_now, x_hat);

				//calculate the first derivative
				matrix_multiply(qp->P, x_now, D1_f0);
				matrix_add_by(D1_f0, qp->q);
				matrix_scale_by(t, D1_f0);

				//calculate squared norm of the first derivative
				curr_norm = 0;
				for(r = 0; r < D1_f0->row; r++) {
					curr_norm += matrix_at(D1_f0, r, 0) * matrix_at(D1_f0, r, 0);
				}

				//VERBOSE_PRINT("[exact line search] #%d: norm = %f, step_size = %f\n",
				//              i, curr_norm, curr_step_size);

				//update best step size so far
				if(curr_norm < smallest_norm) {
					smallest_norm = curr_norm;
					best_step_size = curr_step_size;

					//VERBOSE_PRINT("[exact line search]current best: #%d, step size: %f\n",
					//              i, best_step_size);
				}
			}
			//VERBOSE_PRINT("[exact line search] best step size: %f\n", best_step_size);

			/* make sure newton step won't exceed the inequality bounds */
			bool step_too_large;
			while(1) {
				step_too_large = false;

				/* update the optimization variable */
				for(r = 0; r < z_now->row; r++) {
					//update z with newton step
					matrix_scaling(best_step_size, newton_step, scaled_newton_step);
					matrix_add(z_last, scaled_newton_step, z_now);

					//calculate x from z
					matrix_multiply(F, z_now, x_now);
					matrix_add_by(x_now, x_hat);

#if (ENABLE_LOWER_BOUND_INEQUALITY == 1)
					/*==================================================*
					 * check if lower bound constraints are still valid *
					 *==================================================*/
					if(solve_lower_bound == true) {
						if(matrix_at(x_now, r, 0) < matrix_at(qp->lb, r, 0)) {
							step_too_large = true;
							break;
						}
					}
#endif

#if (ENABLE_UPPER_BOUND_INEQUALITY == 1)
					/*==================================================*
					 * check if upper bound constraints are still valid *
					 *==================================================*/
					if(solve_upper_bound == true) {
						if(matrix_at(x_now, r, 0) > matrix_at(qp->ub, r, 0)) {
							step_too_large = true;
							break;
						}
					}
#endif
				}

#if (ENABLE_AFFINE_INEQUALITY == 1)
				/*========================================================*
				 * check if affine inequality constraints are still valid *
				 *========================================================*/
				if(solve_affine_inequality == true) {
					for(i = 0; i < qp->A->row; i++) {
						fi = 0;
						for(j = 0; j < qp->A->column; j++) {
							fi += matrix_at(qp->A, i, j) * matrix_at(x_now, j, 0);
						}

						if(fi > matrix_at(qp->b, i, 0)) {
							step_too_large = true;
							break;
						}
					}
				}
#endif

				/*================================================================*
				 * schrink the newton step if the step is too large and break the *
				 * inequality constraints                                         *
				 *================================================================*/
				if(step_too_large == false) {
					break;
				} else {
					matrix_scale_by(qp->a, newton_step);
				}
			}

			VERBOSE_PRINT_MATRIX(*Q);
			VERBOSE_PRINT_MATRIX(*R);
			VERBOSE_PRINT_MATRIX(*F);
			VERBOSE_PRINT_MATRIX(*F_t);
			VERBOSE_PRINT_MATRIX(*x_hat);
			VERBOSE_PRINT_MATRIX(*D1_f0);
			VERBOSE_PRINT_MATRIX(*D2_f0);
			VERBOSE_PRINT_MATRIX(*D1_phi);
			VERBOSE_PRINT_MATRIX(*D2_phi);
			VERBOSE_PRINT_MATRIX(*D1_f_tilde);
			VERBOSE_PRINT_MATRIX(*D2_f_tilde);
			VERBOSE_PRINT_MATRIX(*D2_f_tilde_inv);
			VERBOSE_PRINT_MATRIX(*newton_step);
			VERBOSE_PRINT_MATRIX(*z_now);
			VERBOSE_PRINT_MATRIX(*x_now);
			VERBOSE_PRINT("t = %f\n", t);

			qp->iters++;

			FLOAT resid = vector_residual(x_now, x_last);
			VERBOSE_PRINT("residual: %f\n", resid);
			VERBOSE_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->eps || qp->iters == qp->max_iters) {
				matrix_copy(qp->x, x_now);
				break;
			}
		}

		t += qp->mu;
	}

	matrix_delete(D1_f0);
	matrix_delete(D2_f0);
	matrix_delete(D1_fi);
	matrix_delete(D1_fi_t);
	matrix_delete(D1_fi_D1_fi_t);
	matrix_delete(D1_phi);
	matrix_delete(D2_phi);
	matrix_delete(newton_step);
}

int qp_solve_start(qp_t *qp)
{
	if(qp->x == NULL) return QP_ERROR_NO_OPTIMIZATION_VARIABLE;
	if(qp->P == NULL) return QP_ERROR_NO_OBJECTIVE_FUNCTION;
	if(qp->A_eq == NULL && qp->b_eq != NULL) return QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT;
	if(qp->A_eq != NULL && qp->b_eq == NULL) return QP_ERROR_INCOMPLETE_EQUAILITY_CONSTRAINT;
	if(qp->A == NULL && qp->b != NULL) return QP_ERROR_INCOMPLETE_INEQUAILITY_CONSTRAINT;
	if(qp->A != NULL && qp->b == NULL) return QP_ERROR_INCOMPLETE_INEQUAILITY_CONSTRAINT;

	bool solve_equalities = (qp->A_eq != NULL) && (qp->b_eq != NULL);

	bool solve_lower_bound = qp->lb != NULL;
	bool solve_upper_bound = qp->ub != NULL;
	bool solve_affine_inequality = (qp->A != NULL) && (qp->b != NULL);
	bool solve_inequalities = solve_lower_bound | solve_upper_bound | solve_affine_inequality;

	/* no constraint optimization */
	if(!solve_equalities && !solve_inequalities) {
		VERBOSE_PRINT("identify quadratic programming problem without any "
		              "constraints\n");

		qp_solve_no_constraint_problem(qp);
	}

	/* equality constrained optmization */
	if(solve_equalities && !solve_inequalities) {
		VERBOSE_PRINT("identify qudratic programming problem with equality "
		              "constraint\n");

		qp_solve_equality_constraint_problem(qp);
	}

	/* inequality constrained optimization */
	if(!solve_equalities && solve_inequalities) {
		VERBOSE_PRINT("identify qudratic programming problem with inequality "
		              "constraint\n");

		qp_solve_inequality_constraint_problem(qp, solve_lower_bound,
		                                       solve_upper_bound, solve_affine_inequality);
	}

	/* equality-inequality constrained optimization */
	if(solve_equalities && solve_inequalities) {
		VERBOSE_PRINT("identify qudratic programming problem with equality "
		              "and inequality constraints\n");

		qp_solve_equality_inequality_constraint_problem(qp, solve_lower_bound,
		        solve_upper_bound, solve_affine_inequality);
	}

	return QP_SUCCESS_SOLVED;
}
