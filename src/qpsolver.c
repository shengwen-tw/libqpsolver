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

	qp->eps = 1e-6;         //residual value to stop the optimization
	qp->mu = 1.5;           //stiffness growth rate of the log barrier function
	qp->t_init = 0.01;      //initial stiffness of the log barrier
	qp->t_max = 50000;      //maximum stiffness of the log barrier
	qp->max_iters = 10000;  //maximum iteration times

	qp->line_search_num = 10000;
	qp->line_search_min_step_size = 0.1;
}

bool qp_start_point_feasibility_check(qp_t *qp)
{
	FLOAT fi;
	int r, c;

	/* equality constraints */
	if((qp->A_eq != NULL) && (qp->b_eq != NULL)) {
		bool infeasible = false;
		const FLOAT tolerance = 1e-3;

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

float calc_objective_func_val(FLOAT t, qp_t *qp)
{
	FLOAT f = 0;

	/* calculate transpose(x) * P * x */
	matrix_t *xPx = matrix_new(1, 1);
	matrix_t *Px = matrix_new(qp->x->row, qp->x->column);
	matrix_t *xt = matrix_new(qp->x->column, qp->x->row);
	matrix_transpose(qp->x, xt);
	matrix_multiply(qp->P, qp->x, Px);
	matrix_multiply(xt, Px, xPx);

	/* calculate tranpose(q) * r */
	matrix_t *qt_x = matrix_new(1, 1);
	matrix_t *qt = matrix_new(qp->q->column, qp->q->row);
	matrix_transpose(qp->q, qt);
	matrix_multiply(qt, qp->x, qt_x);

	f = matrix_at(xPx, 0, 0) + matrix_at(qt_x, 0, 0);

	matrix_delete(xPx);
	matrix_delete(Px);
	matrix_delete(xt);
	matrix_delete(qt_x);
	matrix_delete(qt);

	return t * f;
}

float calc_log_barrier_objective_func_val(FLOAT t, qp_t *qp, matrix_t *A_inequality)
{
	FLOAT f = 0, sum_phi = 0;

	/* calculate transpose(x) * P * x */
	matrix_t *xPx = matrix_new(1, 1);
	matrix_t *Px = matrix_new(qp->x->row, qp->x->column);
	matrix_t *xt = matrix_new(qp->x->column, qp->x->row);
	matrix_transpose(qp->x, xt);
	matrix_multiply(qp->P, qp->x, Px);
	matrix_multiply(xt, Px, xPx);

	/* calculate tranpose(q) * r */
	matrix_t *qt_x = matrix_new(1, 1);
	matrix_t *qt = matrix_new(qp->q->column, qp->q->row);
	matrix_transpose(qp->q, qt);
	matrix_multiply(qt, qp->x, qt_x);

	f = matrix_at(xPx, 0, 0) + matrix_at(qt_x, 0, 0);

	matrix_delete(xPx);
	matrix_delete(Px);
	matrix_delete(xt);
	matrix_delete(qt_x);
	matrix_delete(qt);

	return t * f + sum_phi;
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

	int r, c;
	int i, j;

	//log barrier's parameter
	FLOAT t = qp->t_init;

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
	FLOAT fi = 0;
	FLOAT div_fi = 0;
	FLOAT div_fi_squared = 0;

	/* inequality constraint matrix (lower bound + upper bound + affine inequality) */
	int lb_size = (solve_lower_bound == true) ? qp->lb->row : 0;
	int ub_size = (solve_upper_bound == true) ? qp->ub->row : 0;
	int b_size = (solve_affine_inequality == true) ? qp->b->row : 0;

	int A_inequality_row = lb_size + ub_size + b_size;
	int A_inequality_column = qp->x->row;
	matrix_t *A_inequality = matrix_zeros(A_inequality_row, A_inequality_column);

	int b_inequality_row = lb_size + ub_size + b_size;
	int b_inequality_column = 1;
	matrix_t *b_inequality = matrix_new(b_inequality_row, b_inequality_column);

	if(solve_lower_bound == true) {
		for(r = 0; r < qp->lb->row; r++) {
			matrix_at(A_inequality, r, r) = -1;
			matrix_at(b_inequality, r, 0) = -matrix_at(qp->lb, r, 0);
		}
	}

	if(solve_upper_bound == true) {
		for(r = 0; r < qp->ub->row; r++) {
			matrix_at(A_inequality, r + lb_size, r) = 1;
			matrix_at(b_inequality, r + lb_size, 0) = matrix_at(qp->ub, r, 0);
		}
	}

	if(solve_affine_inequality == true) {
		for(r = 0; r < qp->A->row; r++) {
			for(c = 0; c < qp->A->column; c++) {
				matrix_at(A_inequality, r + lb_size + ub_size, c) =
				    matrix_at(qp->A, r, c);
			}
			matrix_at(b_inequality, r + lb_size + ub_size, 0) =
			    matrix_at(qp->b, r, 0);
		}
	}

	VERBOSE_PRINT_MATRIX(*A_inequality);
	VERBOSE_PRINT_MATRIX(*b_inequality);

	/* outer loop varies the stiffness of the log barrier functions */
	while(qp->iters < qp->max_iters) {
		if(t > qp->t_max) {
			break;
		}

		/* inner loop do the gradient descent */
		while(qp->iters < qp->max_iters) {
			VERBOSE_PRINT("iteration %d\n", qp->iters + 1);

			/* preseve last x for checking convergence */
			matrix_copy(x_last, qp->x);

			matrix_reset_zeros(D1_phi);
			matrix_reset_zeros(D2_phi);

			/*=======================================================*
			 * calculate first and second derivative of log barriers *
			 *=======================================================*/

			for(r = 0; r < A_inequality->row; r++) {
				/* calculate value of the log barrier function */
				fi = 0;
				for(j = 0; j < A_inequality->column; j++) {
					fi += matrix_at(A_inequality, r, j) * matrix_at(qp->x, j, 0);
				}
				fi -= matrix_at(b_inequality, r, 0);
				div_fi = -1 / (fi + epsilon);
				div_fi_squared = 1 / ((fi * fi) + epsilon);

				/* calculate first derivative of the log barrier function */
				for(i = 0; i < D1_phi->row; i++) {
					matrix_at(D1_phi, i, 0) +=
					    div_fi * matrix_at(A_inequality, r, i);

					matrix_at(D1_fi, i, 0) = matrix_at(A_inequality, r, i);
				}

				/* calculate second derivative of the log barrier function */
				matrix_transpose(D1_fi, D1_fi_t);
				matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);

				for(i = 0; i < D2_phi->row; i++) {
					for(j = 0; j < D2_phi->column; j++) {
						matrix_at(D2_phi, i, j) +=
						    div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j);
					}
				}
			}

			/*=================================================================*
			 * calculate first and second derivative of the objective function *
			 *=================================================================*/

			//calculate the firt derivative of the objective function
			matrix_multiply(qp->P, qp->x, D1_f0);
			matrix_add_by(D1_f0, qp->q);
			matrix_scale_by(t, D1_f0);

			//calculate the second derivate of the objective function
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

			//update the optimization variable
			matrix_add(x_last, newton_step, qp->x);

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

		t *= qp->mu;
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
	//FIXME
	solve_lower_bound = false;
	solve_upper_bound = false;
	solve_affine_inequality = false;

	const FLOAT epsilon = 1e-14; //increase numerical stability of divide by zero

	int r, c;
	int i, j;

	//log barrier's parameter
	FLOAT t = qp->t_init;

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

	/*=================================*
	 * gradient descent step variables *
	 *=================================*/

	//newton step's vector
	vector_t *newton_step_obj = matrix_new(qp->x->row, qp->x->column);
	//newton step times step size
	vector_t *scaled_newton_step_obj = matrix_new(qp->x->row, qp->x->column);

	//newton step's vector
	vector_t *newton_step_barrier = matrix_new(qp->x->row, qp->x->column);
	//newton step times step size
	vector_t *scaled_newton_step_barrier = matrix_new(qp->x->row, qp->x->column);

	/*==============================================================================*
	 * inequality constraint matrix (lower bound + upper bound + affine inequality) *
	 *==============================================================================*/
	int lb_size = (solve_lower_bound == true) ? qp->lb->row : 0;
	int ub_size = (solve_upper_bound == true) ? qp->ub->row : 0;
	int b_size = (solve_affine_inequality == true) ? qp->b->row : 0;

	int A_inequality_row = lb_size + ub_size + b_size;
	int A_inequality_column = qp->x->row;
	matrix_t *A_inequality = matrix_zeros(A_inequality_row, A_inequality_column);

	int b_inequality_row = lb_size + ub_size + b_size;
	int b_inequality_column = 1;
	matrix_t *b_inequality = matrix_new(b_inequality_row, b_inequality_column);

	if(solve_lower_bound == true) {
		for(r = 0; r < qp->lb->row; r++) {
			matrix_at(A_inequality, r, r) = -1;
			matrix_at(b_inequality, r, 0) = -matrix_at(qp->lb, r, 0);
		}
	}

	if(solve_upper_bound == true) {
		for(r = 0; r < qp->ub->row; r++) {
			matrix_at(A_inequality, r + lb_size, r) = 1;
			matrix_at(b_inequality, r + lb_size, 0) = matrix_at(qp->ub, r, 0);
		}
	}

	if(solve_affine_inequality == true) {
		for(r = 0; r < qp->A->row; r++) {
			for(c = 0; c < qp->A->column; c++) {
				matrix_at(A_inequality, r + lb_size + ub_size, c) =
				    matrix_at(qp->A, r, c);
			}
			matrix_at(b_inequality, r + lb_size + ub_size, 0) =
			    matrix_at(qp->b, r, 0);
		}
	}

	/* i-th inenquality constraint function */
	FLOAT fi = 0;
	FLOAT div_fi = 0;
	FLOAT div_fi_squared = 0;

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
	matrix_t *x_last = matrix_zeros(qp->x->row, qp->x->column);

	//initialize x with z
	matrix_multiply(F, z_now, qp->x);
	matrix_add_by(qp->x, x_hat);

	//first derivative of the equality constraints eliminated objective function
	matrix_t *D1_f_tilde = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the equality constraints eliminated objective function
	matrix_t *D2_f_tilde = matrix_new(qp->P->row, qp->P->column);
	//D2_f0 times F
	matrix_t *D2_f0_F = matrix_new(D2_f_tilde->row, D2_f_tilde->column);
	//inverted second derivative of the equality constraints eliminated objective function
	matrix_t *D2_f_tilde_inv = matrix_new(qp->x->row, qp->x->row);

	//first derivative of the equality constraints eliminated objective function
	matrix_t *D1_phi_tilde = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the equality constraints eliminated objective function
	matrix_t *D2_phi_tilde = matrix_new(qp->P->row, qp->P->column);
	//D2_f0 times F
	matrix_t *D2_phi_F = matrix_new(D2_phi_tilde->row, D2_phi_tilde->column);
	//inverted second derivative of the summation of the log barrier functions
	matrix_t *D2_phi_tilde_inv = matrix_new(qp->x->row, qp->x->row);

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
	while(qp->iters < qp->max_iters) {
		if(t > (qp->t_init + qp->t_max)) {
			break;
		}

		/* inner loop do the gradient descent */
		while(qp->iters < qp->max_iters) {
			VERBOSE_PRINT("iteration %d\n", qp->iters + 1);

			/* preseve last x and z variable */
			matrix_copy(z_last, z_now);
			matrix_copy(x_last, qp->x);

			matrix_reset_zeros(D1_phi);
			matrix_reset_zeros(D2_phi);

			/*=======================================================*
			 * calculate first and second derivative of log barriers *
			 *=======================================================*/

			for(r = 0; r < A_inequality->row; r++) {
				/* calculate value of the log barrier function */
				fi = 0;
				for(j = 0; j < A_inequality->column; j++) {
					fi += matrix_at(A_inequality, r, j) * matrix_at(qp->x, j, 0);
				}
				fi -= matrix_at(b_inequality, r, 0);
				div_fi = -1 / (fi + epsilon);
				div_fi_squared = 1 / ((fi * fi) + epsilon);

				/* calculate first derivative of the log barrier function */
				for(i = 0; i < D1_phi->row; i++) {
					matrix_at(D1_phi, i, 0) +=
					    div_fi * matrix_at(A_inequality, r, i);

					matrix_at(D1_fi, i, 0) = matrix_at(A_inequality, r, i);
				}

				/* calculate second derivative of the log barrier function */
				matrix_transpose(D1_fi, D1_fi_t);
				matrix_multiply(D1_fi, D1_fi_t, D1_fi_D1_fi_t);

				for(i = 0; i < D2_phi->row; i++) {
					for(j = 0; j < D2_phi->column; j++) {
						matrix_at(D2_phi, i, j) +=
						    div_fi_squared * matrix_at(D1_fi_D1_fi_t, i, j);
					}
				}
			}

			/*====================================================*
			 * calculate newton step of the log barrier functions *
			 *====================================================*/

			//null space transform of first derivative
			matrix_multiply(F_t, D1_phi, D1_phi_tilde);

			//null space transform of second derivative
			matrix_multiply(D2_phi, F, D2_phi_F);
			matrix_multiply(F_t, D2_phi_F, D2_phi_tilde);

			//calculate newton step
			matrix_inverse(D2_phi_tilde, D2_phi_tilde_inv);
			matrix_multiply(D2_phi_tilde_inv, D1_phi_tilde, newton_step_barrier);
			matrix_scale_by(-1, newton_step_barrier);

			//exact line search
			//TODO

			/*==================================================*
			 * calculate newton step of the objective functions *
			 *==================================================*/

			//first derivative of the objective function
			matrix_multiply(qp->P, qp->x, D1_f0);
			matrix_add_by(D1_f0, qp->q);
			matrix_scale_by(t, D1_f0);

			//second derivative of the objective function
			matrix_copy(D2_f0, qp->P);
			matrix_scale_by(t, D2_f0);

			//null space transform of the first derivative
			matrix_multiply(F_t, D1_f0, D1_f_tilde);

			//null space transform of the second derivative
			matrix_multiply(D2_f0, F, D2_f0_F);
			matrix_multiply(F_t, D2_f0_F, D2_f_tilde);

			//calculate newton step
			matrix_inverse(D2_f_tilde, D2_f_tilde_inv);
			matrix_multiply(D2_f_tilde_inv, D1_f_tilde, newton_step_obj);
			matrix_scale_by(-1, newton_step_obj);

			/* exact line search */
			FLOAT curr_step_size;
			FLOAT best_step_size = 1.0f;
			FLOAT curr_cost = 0;
			FLOAT best_cost = 0;

			/* initialize the step size */
			matrix_add(z_last, newton_step_obj, z_now); //scale = 1
			matrix_multiply(F, z_now, qp->x);
			matrix_add_by(qp->x, x_hat);
			curr_cost = calc_objective_func_val(t, qp);

			/* find the best step size iteratively */
			for(i = (qp->line_search_num - 1); i > 0; i--) {
				curr_step_size = (FLOAT)i / (FLOAT)qp->line_search_num;

				if(curr_step_size <= qp->line_search_min_step_size) {
					break;
				}

				matrix_scaling(curr_step_size, newton_step_obj, scaled_newton_step_obj);
				matrix_add(z_last, scaled_newton_step_obj, z_now);

				//calculate x from z
				matrix_multiply(F, z_now, qp->x);
				matrix_add_by(qp->x, x_hat);

				//calculate (x.' * P * x) + (q.' * r)
				curr_cost = calc_objective_func_val(t, qp);

				//VERBOSE_PRINT("[exact line search] #%d: cost = %f, step_size = %f\n",
				//              i, curr_cost, curr_step_size);

				//update best step size so far
				if(curr_cost < best_cost) {
					best_cost = curr_cost;
					best_step_size = curr_step_size;

					//VERBOSE_PRINT("[exact line search]current best: #%d, step size: %f\n",
					//              i, best_step_size);
				}
			}
			//VERBOSE_PRINT("[exact line search] best step size: %f\n", best_step_size);

			/*==================================*
			 * update the optimization variable *
			 *==================================*/

			//update z with newton step
			matrix_scaling(1, newton_step_obj, scaled_newton_step_obj);
			matrix_add(z_last, scaled_newton_step_obj, z_now);

			matrix_scaling(1, newton_step_barrier, scaled_newton_step_barrier);
			matrix_add_by(z_now, scaled_newton_step_barrier);

			//calculate x from z
			matrix_multiply(F, z_now, qp->x);
			matrix_add_by(qp->x, x_hat);

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
			VERBOSE_PRINT_MATRIX(*newton_step_obj);
			VERBOSE_PRINT_MATRIX(*newton_step_barrier);
			VERBOSE_PRINT_MATRIX(*z_now);
			VERBOSE_PRINT_MATRIX(*qp->x);
			VERBOSE_PRINT("t = %f\n", t);

			qp->iters++;

			FLOAT resid = vector_residual(qp->x, x_last);
			VERBOSE_PRINT("residual: %f\n", resid);
			VERBOSE_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->eps || qp->iters == qp->max_iters) {
				break;
			}
		}

		t *= qp->mu;
	}

	matrix_delete(D1_f0);
	matrix_delete(D2_f0);
	matrix_delete(D1_fi);
	matrix_delete(D1_fi_t);
	matrix_delete(D1_fi_D1_fi_t);
	matrix_delete(D1_phi);
	matrix_delete(D2_phi);
	matrix_delete(newton_step_obj);
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
