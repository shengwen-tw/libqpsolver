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

	/* parameters of phase1 (feasibilty) solver */
	qp->phase1.max_iters = 10000;
	qp->phase1.iters = 0;
	qp->phase1.eps = 1e-5;
	qp->phase1.s_margin = 10;
	qp->phase1.beta = -0.1;
	qp->phase1.t_init = 0.05;
	qp->phase1.t_max = 1000;
	qp->phase1.mu = 2.0;
	qp->phase1.backtracking_alpha = 0.25;
	qp->phase1.backtracking_beta = 0.9;
	//qp->phase1.step_size = 0.1;

	/* parameters of phase2 (quadratic programming) solver */
	qp->phase2.iters = 0;
	qp->phase2.eps = 1e-6;         //residual value to stop the gradient descent inner loop
	qp->phase2.mu = 1.5;           //stiffness growth rate of the log barrier function
	qp->phase2.t_init = 0.01;      //initial stiffness of the log barrier
	qp->phase2.t_max = 50000;      //maximum stiffness of the log barrier
	qp->phase2.max_iters = 10000;  //maximum iteration times
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
	VERBOSE_PRINT("[solver] problem type: unconstrained QP\n");

	/* the closed form solution is given by setting the first derivative equal
	 * to zero, i.e: Px = -q */

	/* construct -q vector */
	int r;
	for(r = 0; r < qp->q->row; r++) {
		qp->q->data[r] *= -1;
	}

	/* solve Px = -q */
	solve_linear_system(qp->P, qp->x, qp->q);
}

static void qp_solve_equality_constraint_problem(qp_t *qp)
{
	VERBOSE_PRINT("[solver] problem type: equality constrainted QP\n");

	/* the closed form solution of the problem can be obtained by solving the *
	 * KKT system, i.e: [P  A.'][ x*] = [-q]                                  *
	 *                  [A   0 ][nu*]   [ b]                                  *
	 * where x* is the optimal solution and nu* is the optimal dual variable  */

	int r, c;

	/* construct the [-q; b] matrix */
	int qb_vec_row = qp->q->row + qp->b_eq->row;
	int qb_vec_column = 1;
	matrix_t *qb_vec = matrix_new(qb_vec_row, qb_vec_column);

	//copy -q
	for(r = 0; r < qp->q->row; r++) {
		//copy -q
		matrix_at(qb_vec, r, 0) = -matrix_at(qp->q, r, 0);
	}
	//copy b
	for(r = 0; r < qp->b_eq->row; r++) {
		//copy b
		matrix_at(qb_vec, r + qp->q->row, 0) =
		    matrix_at(qp->b_eq, r, 0);
	}

	/* construct the KKT matrix */
	int kkt_row = qp->P->row + qp->A_eq->row;
	int kkt_column = qp->P->column + qp->A_eq->row;
	matrix_t * KKT = matrix_zeros(kkt_row, kkt_column);

	//copy P
	for(r = 0; r < qp->P->row; r++) {
		for(c = 0; c < qp->P->column; c++) {
			matrix_at(KKT, r, c) = matrix_at(qp->P, r, c);
		}
	}
	//copy A.'
	for(r = 0; r < qp->A_eq->column; r++) {
		for(c = 0; c < qp->A_eq->row; c++) {
			matrix_at(KKT, r, (c + qp->A_eq->column)) =
			    matrix_at(qp->A_eq, c, r);
		}
	}
	//copy A
	for(r = 0; r < qp->A_eq->row; r++) {
		for(c = 0; c < qp->A_eq->column; c++) {
			matrix_at(KKT, (r + qp->P->row), c) =
			    matrix_at(qp->A_eq, r, c);
		}
	}
	//set zero matrix
	for(r = 0; r < qp->A_eq->row; r++) {
		for(c = 0; c < qp->A_eq->row; c++) {
			matrix_at(KKT, (r + qp->P->row), (c + qp->A_eq->column)) = 0;
		}
	}

	/* construct kkt solution vector */
	int kkt_sol_row = qp->x->row + qp->b_eq->row;
	int kkt_sol_column = 1;
	matrix_t *kkt_sol = matrix_new(kkt_sol_row, kkt_sol_column);

	/* solve the KKT system */
	solve_linear_system(KKT, kkt_sol, qb_vec);

	/* copy the optimal solution back to x */
	for(r = 0; r < qp->x->row; r++) {
		matrix_at(qp->x, r, 0) = matrix_at(kkt_sol, r, 0);
	}

	DEBUG_PRINT_MATRIX(*qb_vec);
	DEBUG_PRINT_MATRIX(*KKT);
	DEBUG_PRINT_MATRIX(*kkt_sol);

	/* free mallocs */
	matrix_delete(KKT);
	matrix_delete(KKT);
	matrix_delete(qb_vec);
}

static FLOAT qp_phase1_cost_function(FLOAT t, matrix_t *x_prime,
                                     matrix_t *A_inequality, matrix_t* b_inequality,
                                     FLOAT beta, FLOAT s_min)
{
	int r, c;

	FLOAT f, fi;
	FLOAT div_by_t = 1.0f / t;

	/* cost of slack variable s */
	f = matrix_at(x_prime, x_prime->row - 1, 0);

	/* cost of log barrier functions for x */
	for(r = 0; r < (A_inequality->row - 1); r++) {
		/* calculate value of the inequality functions */
		fi = 0;
		for(c = 0; c < (A_inequality->column - 1); c++) {
			fi += matrix_at(A_inequality, r, c) * matrix_at(x_prime, c, 0);
		}
		fi = fi - matrix_at(b_inequality, r, 0) -
		     matrix_at(x_prime, x_prime->row-1, 0);

		f -= div_by_t * log(-fi);
	}

	/* cost of log barrier function for s */
	fi = -matrix_at(x_prime, x_prime->row - 1, 0) + beta;
	f -= div_by_t * log(-fi);

	/* cost of log barrier function for limiting s not cross over s */
	fi = -matrix_at(x_prime, x_prime->row - 1, 0) + s_min;
	f -= div_by_t * log(-fi);

	return f;
}

static int qp_inequality_constraint_phase1(qp_t *qp, bool solve_lower_bound,
        bool solve_upper_bound, bool solve_affine_inequality)
{
	VERBOSE_PRINT("[solver] infeasible start point, phase1 start\n");

	const FLOAT epsilon = 1e-14; //increase numerical stability of divide by zero

	int r, c;
	int i, j;

	//log barrier's parameter
	FLOAT t = qp->phase1.t_init;
	FLOAT div_by_t;

	//new optimization vector = original vector + slack variable s
	matrix_t *x_prime = matrix_new(qp->x->row + 1, qp->x->column);
	matrix_t *x_prime_last = matrix_new(qp->x->row + 1, qp->x->column);

	for(r = 0; r < qp->x->row; r++) {
		matrix_at(x_prime, r, 0) = matrix_at(qp->x, r, 0);
		matrix_at(x_prime_last, r, 0) = matrix_at(qp->x, r, 0);
	}

	//first derivative of the objective function
	matrix_t *D1_f0 = matrix_zeros(qp->x->row + 1, qp->x->column);

	//first derivative of the sumation of log barrier functions
	matrix_t *D1_phi_x = matrix_new(qp->x->row + 1, qp->x->column);
	matrix_t *D1_phi_s = matrix_zeros(qp->x->row + 1, qp->x->column);

	/* i-th inenquality constraint function */
	FLOAT fi = 0;
	FLOAT div_fi = 0;

	//gradient descent step
	vector_t *descent_step = matrix_new(qp->x->row + 1, qp->x->column);

	/* inequality constraint matrix (lower bound + upper bound + affine inequality) */
	int lb_size = (solve_lower_bound == true) ? qp->lb->row : 0;
	int ub_size = (solve_upper_bound == true) ? qp->ub->row : 0;
	int b_size = (solve_affine_inequality == true) ? qp->b->row : 0;

	int A_inequality_row = lb_size + ub_size + b_size + 1;
	int A_inequality_column = qp->x->row + 1;
	matrix_t *A_inequality = matrix_zeros(A_inequality_row, A_inequality_column);

	int b_inequality_row = lb_size + ub_size + b_size + 1;
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

	/*============================================*
	 * initialize the inequality slack variable s *
	 *============================================*/
	FLOAT fi_max = 0;

	/* initial value of fi_max */
	for(j = 0; j < A_inequality->column - 1; j++) {
		fi_max += matrix_at(A_inequality, 0, j) * matrix_at(qp->x, j, 0);
	}
	fi_max -= matrix_at(b_inequality, 0, 0);

	/* search for the largest fi_max */
	for(r = 1; r < b_inequality->row - 1; r++) {
		/* calculate value of the log barrier function */
		fi = 0;
		for(c = 0; c < A_inequality->column; c++) {
			fi += matrix_at(A_inequality, r, c) * matrix_at(qp->x, c, 0);
		}
		fi -= matrix_at(b_inequality, r, 0);

		if(fi > fi_max) {
			fi_max = fi;
		}
	}

	//initialize slack variable s
	matrix_at(x_prime, qp->x->row, 0) = fi_max + qp->phase1.s_margin;

	//lower bound of s can be for now (avoiding s cross over the inequality constraints of x)
	FLOAT s_min_now = fi;

	/*====================*
	 * optimization start *
	 *====================*/

	while(qp->phase1.iters < qp->phase1.max_iters) {
		if(t > (qp->phase1.t_init + qp->phase1.t_max)) {
			break;
		}

		if(matrix_at(x_prime, x_prime->row - 1, 0) < 0) {
			break;
		}

		div_by_t = 1 / t;
		while(qp->phase1.iters < qp->phase1.max_iters) {
			DEBUG_PRINT("iteration %d\n", qp->phase1.iters + 1);

			matrix_copy(x_prime_last, x_prime);

			matrix_reset_zeros(D1_f0);
			matrix_reset_zeros(D1_phi_x);

			/* first derivative of the objective function */
			matrix_at(D1_f0, x_prime->row-1, 0) = 1;

			/* first derivarive of x's inequality constraints */
			for(r = 0; r < (A_inequality->row - 1); r++) {
				/* calculate value of the inequality functions */
				fi = 0;
				for(j = 0; j < (A_inequality->column - 1); j++) {
					fi += matrix_at(A_inequality, r, j) * matrix_at(x_prime, j, 0);
				}
				fi = fi - matrix_at(b_inequality, r, 0) -
				     matrix_at(x_prime, x_prime->row-1, 0);
				div_fi = -1 / (fi + epsilon);

				/* calculate first derivative of the log barrier function */
				for(i = 0; i < D1_phi_x->row - 1; i++) {
					matrix_at(D1_phi_x, i, 0) +=
					    div_by_t * div_fi * matrix_at(A_inequality, r, i);
				}
			}

			/*===============================================*
			 * first derivative of s's inequality constraint *
			 *===============================================*/

			//stop minimization when s is negative enough
			fi = matrix_at(x_prime, x_prime->row - 1, 0) - qp->phase1.beta;
			div_fi = -1.0 / (fi + epsilon);
			matrix_at(D1_phi_s, qp->x->row, 0) = div_by_t * div_fi;

			//avoiding s cross over the inequality constraints of x
			fi = matrix_at(x_prime, x_prime->row - 1, 0) - s_min_now;
			div_fi = -1.0 / (fi + epsilon);
			matrix_at(D1_phi_s, qp->x->row, 0) += div_by_t * div_fi;

			/* combine first derivative of objective function and log barriers */
			matrix_add_by(D1_f0, D1_phi_x);
			matrix_add_by(D1_f0, D1_phi_s);

			/*================================================*
			 * gradient descent with backtracking line search *
			 *================================================*/

			/* initialization of backtracking line search */
			FLOAT backtracking_t = 1.0f;
			FLOAT bt_cost_now;

			//f(x)
			FLOAT bt_cost_origin = qp_phase1_cost_function(t, x_prime, A_inequality,
			                       b_inequality, qp->phase1.beta, s_min_now);
			//f(x) + (a * t * D1_f(x).' * D1_f(x))
			FLOAT bt_cost_alpha_line = 0;
			//(a * t * D1_f(x).' * D1_f(x))
			FLOAT alpha_line_change = 0;

			/* backtracking line search loop */
			while(1) {
				alpha_line_change = 0;
				for(r = 0; r < x_prime->row; r++) {
					alpha_line_change += matrix_at(D1_f0, r, 0) * matrix_at(D1_f0, r, 0);
				}
				alpha_line_change *= qp->phase1.backtracking_alpha * backtracking_t;
				bt_cost_alpha_line = bt_cost_origin - alpha_line_change;

				matrix_scaling(-backtracking_t, D1_f0, descent_step);
				matrix_add(x_prime_last, descent_step, x_prime);
				bt_cost_now =
				    qp_phase1_cost_function(t, x_prime, A_inequality,
				                            b_inequality, qp->phase1.beta, s_min_now);

				if(bt_cost_now <= bt_cost_alpha_line) {
					break;
				}

				backtracking_t *= qp->phase1.backtracking_beta;
			}

			/*===================================================*
			 * gradient descent without backtracking line search *
			 *===================================================*/
			//matrix_scaling(-qp->phase1.step_size, D1_f0, descent_step);
			//matrix_add(x_prime_last, descent_step, x_prime);

			/*=====================*
			 * find new s boundary *
			 *=====================*/

			/* initial value of fi_max */
			fi_max = 0;
			for(j = 0; j < A_inequality->column - 1; j++) {
				fi_max += matrix_at(A_inequality, 0, j) * matrix_at(x_prime, j, 0);
			}
			fi_max -= matrix_at(b_inequality, 0, 0);

			/* search for the largest fi_max */
			for(r = 1; r < b_inequality->row - 1; r++) {
				/* calculate value of the log barrier function */
				fi = 0;
				for(c = 0; c < A_inequality->column; c++) {
					fi += matrix_at(A_inequality, r, c) * matrix_at(x_prime, c, 0);
				}
				fi -= matrix_at(b_inequality, r, 0);

				if(fi > fi_max) {
					fi_max = fi;
				}
			}
			s_min_now = fi;

			DEBUG_PRINT_VAR(s_min_now);

			qp->phase1.iters++;

			FLOAT resid = vector_residual(x_prime, x_prime_last);

			DEBUG_PRINT_MATRIX(*D1_f0);
			DEBUG_PRINT_MATRIX(*D1_phi_x);
			DEBUG_PRINT_MATRIX(*D1_phi_s);
			DEBUG_PRINT_VAR(t);
			DEBUG_PRINT_MATRIX(*x_prime);
			DEBUG_PRINT_VAR(resid);
			DEBUG_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->phase1.eps || qp->phase1.iters == qp->phase1.max_iters) {
				break;
			}
		}
		t *= qp->phase1.mu;
	}

	for(r = 0; r < qp->x->row; r++) {
		matrix_at(qp->x, r, 0) = matrix_at(x_prime, r, 0);
	}

	int ret_val;

	//if s <= 0 then the problem is feasible
	if(matrix_at(x_prime, x_prime->row-1, 0) < 0) {
		ret_val = QP_PHASE1_FEASIBLE;
	} else {
		ret_val = QP_PHASE1_INFEASIBLE;
	}

	matrix_delete(x_prime);
	matrix_delete(x_prime_last);
	matrix_delete(D1_f0);
	matrix_delete(D1_phi_x);
	matrix_delete(D1_phi_s);
	matrix_delete(descent_step);
	matrix_delete(A_inequality);
	matrix_delete(b_inequality);

	return ret_val;
}

static void qp_solve_inequality_constraint_problem(qp_t *qp, bool solve_lower_bound,
        bool solve_upper_bound, bool solve_affine_inequality)
{
	VERBOSE_PRINT("[solver] problem type: inequality constrained QP\n");

#if (ENABLE_INFEASIBLE_START != 0)
	bool feasible = qp_start_point_feasibility_check(qp);

	if(feasible == false) {
		int phase1 = qp_inequality_constraint_phase1(qp, solve_lower_bound,
		             solve_upper_bound, solve_affine_inequality);

		if(phase1 == QP_PHASE1_FEASIBLE) {
			VERBOSE_PRINT("[solver] phase1 end, feasible start point found\n");
			VERBOSE_PRINT("[solver] phase2 start\n");
		} else {
			VERBOSE_PRINT("[solver] phase1 end: infeasible problem\n"
			              "[solver] abort\n");
			return;
		}
	}
#endif

	const FLOAT epsilon = 1e-14; //increase numerical stability of divide by zero

	int r, c;
	int i, j;

	//log barrier's parameter
	FLOAT t = qp->phase2.t_init;

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

	/* outer loop varies the stiffness of the log barrier functions */
	while(qp->phase2.iters < qp->phase2.max_iters) {
		if(t > qp->phase2.t_max) {
			break;
		}

		/* inner loop do the gradient descent */
		while(qp->phase2.iters < qp->phase2.max_iters) {
			DEBUG_PRINT("iteration %d\n", qp->phase2.iters + 1);

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

			qp->phase2.iters++;
			FLOAT resid = vector_residual(qp->x, x_last);

			DEBUG_PRINT_MATRIX(*D1_phi);
			DEBUG_PRINT_MATRIX(*D2_phi);
			DEBUG_PRINT_MATRIX(*D1_f0);
			DEBUG_PRINT_MATRIX(*D2_f0);
			DEBUG_PRINT_MATRIX(*D2_f0_inv);
			DEBUG_PRINT_MATRIX(*newton_step);
			DEBUG_PRINT_MATRIX(*qp->x);
			DEBUG_PRINT_VAR(t);
			DEBUG_PRINT_VAR(resid);
			DEBUG_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->phase2.eps) {
				break;
			}
		}

		t *= qp->phase2.mu;
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

static int qp_equality_inequality_constraint_phase1(qp_t *qp, bool solve_lower_bound,
        bool solve_upper_bound, bool solve_affine_inequality)
{
	VERBOSE_PRINT("[solver] infeasible start point, phase1 start\n");

	const FLOAT epsilon = 1e-14; //increase numerical stability of divide by zero

	int r, c;
	int i, j;

	//log barrier's parameter
	FLOAT t = qp->phase1.t_init;
	FLOAT div_by_t;

	//new optimization vector = original vector + slack variable s
	matrix_t *z_prime = matrix_zeros(qp->x->row + 1, qp->x->column);
	matrix_t *z_prime_last = matrix_new(qp->x->row + 1, qp->x->column);

	//vectos for evaulating the objective function during the backtracking line search
	matrix_t *z = matrix_new(qp->x->row, qp->x->column);
	matrix_t *Fz = matrix_new(qp->x->row, qp->x->column);
	matrix_t *x_prime = matrix_new(qp->x->row + 1, qp->x->column);

	//first derivative of the objective function
	matrix_t *D1_f0 = matrix_zeros(qp->x->row + 1, qp->x->column);

	//first derivative of the sumation of log barrier functions
	matrix_t *D1_phi_z = matrix_new(qp->x->row, qp->x->column);
	matrix_t *D1_phi_s = matrix_zeros(qp->x->row + 1, qp->x->column);

	matrix_t *D1_phi_tilde_z = matrix_new(qp->x->row, qp->x->column);

	/* i-th inenquality constraint function */
	FLOAT fi = 0;
	FLOAT div_fi = 0;

	//gradient descent step
	vector_t *descent_step = matrix_new(qp->x->row + 1, qp->x->column);

	/* inequality constraint matrix (lower bound + upper bound + affine inequality) */
	int lb_size = (solve_lower_bound == true) ? qp->lb->row : 0;
	int ub_size = (solve_upper_bound == true) ? qp->ub->row : 0;
	int b_size = (solve_affine_inequality == true) ? qp->b->row : 0;

	int A_inequality_row = lb_size + ub_size + b_size + 1;
	int A_inequality_column = qp->x->row + 1;
	matrix_t *A_inequality = matrix_zeros(A_inequality_row, A_inequality_column);

	int b_inequality_row = lb_size + ub_size + b_size + 1;
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

	/*========================================================*
	 * eliminate equality constraints by null space transform *
	 *========================================================*/

	matrix_t *F = matrix_zeros(qp->A_eq->column, qp->A_eq->column);
	matrix_t *F_t = matrix_new(F->column, F->column);
	matrix_t *Q, *R;

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

	/*============================================*
	 * initialize the inequality slack variable s *
	 *============================================*/
	FLOAT fi_max = 0;

	/* initial value of fi_max */
	for(j = 0; j < A_inequality->column - 1; j++) {
		fi_max += matrix_at(A_inequality, 0, j) * matrix_at(qp->x, j, 0);
	}
	fi_max -= matrix_at(b_inequality, 0, 0);

	/* search for the largest fi_max */
	for(r = 1; r < b_inequality->row - 1; r++) {
		/* calculate value of the log barrier function */
		fi = 0;
		for(c = 0; c < A_inequality->column; c++) {
			fi += matrix_at(A_inequality, r, c) * matrix_at(qp->x, c, 0);
		}
		fi -= matrix_at(b_inequality, r, 0);

		if(fi > fi_max) {
			fi_max = fi;
		}
	}

	//initialize slack variable s
	matrix_at(z_prime, z_prime->row - 1, 0) = fi_max + qp->phase1.s_margin;

	//lower bound of s can be for now (avoiding s cross over the inequality constraints of x)
	FLOAT s_min_now = fi;

	/*====================*
	 * optimization start *
	 *====================*/

	while(qp->phase1.iters < qp->phase1.max_iters) {
		if(t > (qp->phase1.t_init + qp->phase1.t_max)) {
			break;
		}

		if(matrix_at(z_prime, z_prime->row - 1, 0) < 0) {
			break;
		}

		div_by_t = 1 / t;
		while(qp->phase1.iters < qp->phase1.max_iters) {
			DEBUG_PRINT("iteration %d\n", qp->phase1.iters + 1);

			matrix_copy(z_prime_last, z_prime);

			matrix_reset_zeros(D1_f0);
			matrix_reset_zeros(D1_phi_z);

			/* first derivative of the objective function */
			matrix_at(D1_f0, z_prime->row-1, 0) = 1;

			//copy z from z_prime (XXX: this code could be optmized)
			for(r = 0; r < z->row; r++) {
				matrix_at(z, r, 0) = matrix_at(z_prime, r, 0);
			}

			//calculate x = Fz
			matrix_multiply(F, z, Fz);

			/* first derivarive of z's inequality constraints */
			for(r = 0; r < (A_inequality->row - 1); r++) {
				/* calculate value of the inequality functions */
				fi = 0;
				for(j = 0; j < (A_inequality->column - 1); j++) {
					fi += matrix_at(A_inequality, r, j) * matrix_at(Fz, j, 0);
				}
				fi = fi - matrix_at(b_inequality, r, 0) -
				     matrix_at(z_prime, z_prime->row-1, 0);
				div_fi = -1 / (fi + epsilon);

				/* calculate first derivative of the log barrier function */
				for(i = 0; i < D1_phi_z->row - 1; i++) {
					matrix_at(D1_phi_z, i, 0) +=
					    div_by_t * div_fi * matrix_at(A_inequality, r, i);
				}
			}
			matrix_multiply(F_t, D1_phi_z, D1_phi_tilde_z);

			/*===============================================*
			 * first derivative of s's inequality constraint *
			 *===============================================*/

			//stop minimization when s is negative enough
			fi = matrix_at(z_prime, z_prime->row - 1, 0) - qp->phase1.beta;
			div_fi = -1.0 / (fi + epsilon);
			matrix_at(D1_phi_s, qp->x->row, 0) = div_by_t * div_fi;

			//avoiding s cross over the inequality constraints of x
			fi = matrix_at(z_prime, z_prime->row - 1, 0) - s_min_now;
			div_fi = -1.0 / (fi + epsilon);
			matrix_at(D1_phi_s, qp->x->row, 0) += div_by_t * div_fi;

			/* combine first derivative of objective function and log barriers */
			//first derivative from z's inequality
			for(r = 0; r < D1_phi_tilde_z->row; r++) {
				matrix_at(D1_f0, r, 0) += matrix_at(D1_phi_tilde_z, r, 0);
			}
			//first derivative from s's inequality
			matrix_add_by(D1_f0, D1_phi_s);

			/*================================================*
			 * gradient descent with backtracking line search *
			 *================================================*/

			/* initialization of backtracking line search */
			FLOAT backtracking_t = 1.0f;
			FLOAT bt_cost_now;

			//copy z from z_prime (XXX: this code could be optmized)
			for(r = 0; r < z->row; r++) {
				matrix_at(z, r, 0) = matrix_at(z_prime, r, 0);
			}

			//calculate x = Fz
			matrix_multiply(F, z, Fz);

			//x_prime = [Fz; s] (XXX: this code could be optmized)
			for(r = 0; r < x_prime->row - 1; r++) {
				matrix_at(x_prime, r, 0) = matrix_at(Fz, r, 0);
			}
			matrix_at(x_prime, x_prime->row - 1, 0) =
			    matrix_at(z_prime, z_prime->row - 1, 0);

			//f(x)
			FLOAT bt_cost_origin = qp_phase1_cost_function(t, x_prime, A_inequality,
			                       b_inequality, qp->phase1.beta, s_min_now);
			//f(x) + (a * t * D1_f(x).' * D1_f(x))
			FLOAT bt_cost_alpha_line = 0;
			//(a * t * D1_f(x).' * D1_f(x))
			FLOAT alpha_line_change = 0;

			/* backtracking line search loop */
			while(1) {
				/* calculate backtraking stop line */
				alpha_line_change = 0;
				for(r = 0; r < z_prime->row; r++) {
					alpha_line_change += matrix_at(D1_f0, r, 0) * matrix_at(D1_f0, r, 0);
				}
				alpha_line_change *= qp->phase1.backtracking_alpha * backtracking_t;
				bt_cost_alpha_line = bt_cost_origin - alpha_line_change;

				/* calculate f(x + t * delta_x) */
				matrix_scaling(-backtracking_t, D1_f0, descent_step);
				matrix_add(z_prime_last, descent_step, z_prime);

				//copy z from z_prime (XXX: this code could be optmized)
				for(r = 0; r < z->row; r++) {
					matrix_at(z, r, 0) = matrix_at(z_prime, r, 0);
				}

				//x = Fz
				matrix_multiply(F, z, Fz);

				//x_prime = [Fz; s] (XXX: this code could be optmized)
				for(r = 0; r < x_prime->row - 1; r++) {
					matrix_at(x_prime, r, 0) = matrix_at(Fz, r, 0);
				}
				matrix_at(x_prime, x_prime->row - 1, 0) =
				    matrix_at(z_prime, z_prime->row - 1, 0);

				//f(x + t * delta_x)
				bt_cost_now =
				    qp_phase1_cost_function(t, x_prime, A_inequality,
				                            b_inequality, qp->phase1.beta, s_min_now);

				if(bt_cost_now <= bt_cost_alpha_line) {
					break;
				}

				backtracking_t *= qp->phase1.backtracking_beta;
			}

			/*===================================================*
			 * gradient descent without backtracking line search *
			 *===================================================*/
			//matrix_scaling(-qp->phase1.step_size, D1_f0, descent_step);
			//matrix_add(z_prime_last, descent_step, z_prime);

			/*=====================*
			 * find new s boundary *
			 *=====================*/

			/* initial value of fi_max */
			fi_max = 0;
			for(j = 0; j < A_inequality->column - 1; j++) {
				fi_max += matrix_at(A_inequality, 0, j) * matrix_at(Fz, j, 0);
			}
			fi_max -= matrix_at(b_inequality, 0, 0);

			/* search for the largest fi_max */
			for(r = 1; r < b_inequality->row - 1; r++) {
				/* calculate value of the log barrier function */
				fi = 0;
				for(c = 0; c < A_inequality->column; c++) {
					fi += matrix_at(A_inequality, r, c) * matrix_at(Fz, c, 0);
				}
				fi -= matrix_at(b_inequality, r, 0);

				if(fi > fi_max) {
					fi_max = fi;
				}
			}
			s_min_now = fi;

			DEBUG_PRINT_VAR(s_min_now);

			qp->phase1.iters++;

			FLOAT resid = vector_residual(z_prime, z_prime_last);

			DEBUG_PRINT_MATRIX(*D1_f0);
			DEBUG_PRINT_MATRIX(*D1_phi_z);
			DEBUG_PRINT_MATRIX(*D1_phi_s);
			DEBUG_PRINT_VAR(t);
			DEBUG_PRINT_MATRIX(*z_prime);
			DEBUG_PRINT_VAR(resid);
			DEBUG_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->phase1.eps || qp->phase1.iters == qp->phase1.max_iters) {
				break;
			}
		}
		t *= qp->phase1.mu;
	}

	//return the optimal x vector for feasibility optmization
	matrix_multiply(F, z_prime, qp->x);

	int ret_val;

	//if s <= 0 then the problem is feasible
	if(matrix_at(z_prime, z_prime->row-1, 0) < 0) {
		ret_val = QP_PHASE1_FEASIBLE;
	} else {
		ret_val = QP_PHASE1_INFEASIBLE;
	}

	matrix_delete(z_prime);
	matrix_delete(z_prime_last);
	matrix_delete(D1_f0);
	matrix_delete(D1_phi_z);
	matrix_delete(D1_phi_s);
	matrix_delete(descent_step);
	matrix_delete(A_inequality);
	matrix_delete(b_inequality);

	return ret_val;
}

static void qp_solve_equality_inequality_constraint_problem(qp_t *qp, bool solve_lower_bound,
        bool solve_upper_bound, bool solve_affine_inequality)
{
	VERBOSE_PRINT("[solver] problem type: equality and inequality constrained QP\n");

#if (ENABLE_INFEASIBLE_START != 0)
	bool feasible = qp_start_point_feasibility_check(qp);

	if(feasible == false) {
		int phase1 = qp_equality_inequality_constraint_phase1(qp, solve_lower_bound,
		             solve_upper_bound, solve_affine_inequality);

		if(phase1 == QP_PHASE1_FEASIBLE) {
			VERBOSE_PRINT("[solver] phase1 end, feasible start point found\n");
			VERBOSE_PRINT("[solver] phase2 start\n");
		} else {
			VERBOSE_PRINT("[solver] phase1 end: infeasible problem\n"
			              "[solver] abort\n");
			return;
		}
	}
#endif

	int r, c;
	int i, j;

	/* increase numerical stability */
	//avoid divide by zero
	const FLOAT epsilon = 1e-14;
	//avoid inversion of singular matrix
	matrix_t *epsilon_matrix = matrix_new(qp->x->row, qp->x->row);
	for(r = 0; r < epsilon_matrix->row; r++) {
		for(c = 0; c < epsilon_matrix->column; c++) {
			matrix_at(epsilon_matrix, r, c) = epsilon;
		}
	}

	//log barrier's parameter
	FLOAT t = qp->phase2.t_init;

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

	//newton step
	vector_t *newton_step = matrix_new(qp->x->row, qp->x->column);

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

	//new optimization variable z
	matrix_t *z_now = matrix_zeros(qp->x->row, qp->x->column);
	matrix_t *z_last = matrix_zeros(qp->x->row, qp->x->column);

	//x = Fx
	matrix_t *x_last = matrix_zeros(qp->x->row, qp->x->column);

	//initialize x with z
	matrix_multiply(F, z_now, qp->x);

	//first derivative of the equality constraints eliminated objective function
	matrix_t *D1_f_tilde = matrix_new(qp->x->row, qp->x->column);
	//second derivative of the equality constraints eliminated objective function
	matrix_t *D2_f_tilde = matrix_new(qp->P->row, qp->P->column);
	//D2_f0 times F
	matrix_t *D2_f0_F = matrix_new(D2_f_tilde->row, D2_f_tilde->column);
	//inverted second derivative of the equality constraints eliminated objective function
	matrix_t *D2_f_tilde_inv = matrix_new(qp->x->row, qp->x->row);

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
	while(qp->phase2.iters < qp->phase2.max_iters) {
		if(t > (qp->phase2.t_init + qp->phase2.t_max)) {
			break;
		}

		/* inner loop do the gradient descent */
		while(qp->phase2.iters < qp->phase2.max_iters) {
			DEBUG_PRINT("iteration %d\n", qp->phase2.iters + 1);

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

			//combine derivatives of objective function and log barrier functions
			matrix_add_by(D1_f0, D1_phi);
			matrix_add_by(D2_f0, D2_phi);

			//null space transform of the first derivative
			matrix_multiply(F_t, D1_f0, D1_f_tilde);

			//null space transform of the second derivative
			matrix_multiply(D2_f0, F, D2_f0_F);
			matrix_multiply(F_t, D2_f0_F, D2_f_tilde);

			//calculate newton step
			matrix_add_by(D2_f_tilde, epsilon_matrix); //increase numerical stability
			matrix_inverse(D2_f_tilde, D2_f_tilde_inv);
			matrix_multiply(D2_f_tilde_inv, D1_f_tilde, newton_step);
			matrix_scale_by(-1, newton_step);

			/*==================================*
			 * update the optimization variable *
			 *==================================*/

			//update z with newton step
			matrix_add(z_last, newton_step, z_now);

			//calculate x from z
			matrix_multiply(F, z_now, qp->x);

			qp->phase2.iters++;

			FLOAT resid = vector_residual(qp->x, x_last);

			DEBUG_PRINT_MATRIX(*Q);
			DEBUG_PRINT_MATRIX(*R);
			DEBUG_PRINT_MATRIX(*F);
			DEBUG_PRINT_MATRIX(*F_t);
			DEBUG_PRINT_MATRIX(*D1_f0);
			DEBUG_PRINT_MATRIX(*D2_f0);
			DEBUG_PRINT_MATRIX(*D1_phi);
			DEBUG_PRINT_MATRIX(*D2_phi);
			DEBUG_PRINT_MATRIX(*D1_f_tilde);
			DEBUG_PRINT_MATRIX(*D2_f_tilde);
			DEBUG_PRINT_MATRIX(*D2_f_tilde_inv);
			DEBUG_PRINT_MATRIX(*newton_step);
			DEBUG_PRINT_MATRIX(*z_now);
			DEBUG_PRINT_MATRIX(*qp->x);
			DEBUG_PRINT_VAR(t);
			DEBUG_PRINT_VAR(resid);
			DEBUG_PRINT("---\n");

			/* exit if already converged */
			if(resid < qp->phase2.eps || qp->phase2.iters == qp->phase2.max_iters) {
				break;
			}
		}

		t *= qp->phase2.mu;
	}

	matrix_delete(D1_f0);
	matrix_delete(D2_f0);
	matrix_delete(D1_fi);
	matrix_delete(D1_fi_t);
	matrix_delete(D1_fi_D1_fi_t);
	matrix_delete(D1_phi);
	matrix_delete(D2_phi);
	matrix_delete(newton_step);
	matrix_delete(A_inequality);
	matrix_delete(b_inequality);
	matrix_delete(F);
	matrix_delete(F_t);
	matrix_delete(Q);
	matrix_delete(R);
	matrix_delete(z_now);
	matrix_delete(z_last);
	matrix_delete(x_last);
	matrix_delete(D1_f_tilde);
	matrix_delete(D2_f_tilde);
	matrix_delete(D2_f0_F);
	matrix_delete(D2_f_tilde_inv);
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
		qp_solve_no_constraint_problem(qp);
	}

	/* equality constrained optmization */
	if(solve_equalities && !solve_inequalities) {
		qp_solve_equality_constraint_problem(qp);
	}

	/* inequality constrained optimization */
	if(!solve_equalities && solve_inequalities) {
		qp_solve_inequality_constraint_problem(qp, solve_lower_bound,
		                                       solve_upper_bound, solve_affine_inequality);
	}

	/* equality-inequality constrained optimization */
	if(solve_equalities && solve_inequalities) {
		qp_solve_equality_inequality_constraint_problem(qp, solve_lower_bound,
		        solve_upper_bound, solve_affine_inequality);
	}

	return QP_SUCCESS_SOLVED;
}
