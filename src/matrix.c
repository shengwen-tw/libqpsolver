#include <string.h>
#include <math.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

int matrix_rank(matrix_t* mat)
{
	matrix_t *mat_copy = matrix_new(mat->row, mat->column);
	matrix_copy(mat_copy, mat);

	int rank = 0;
	double epsilon = 1e-6;

	matrix_t *S = matrix_new(mat_copy->column, 1);
	matrix_t *superb = matrix_new(mat_copy->column, 1);

	LAPACKE_dgesvd(LAPACK_ROW_MAJOR,
	               'N',
	               'N',
	               mat_copy->row,
	               mat_copy->column,
	               mat_copy->data,
	               mat_copy->column,
	               S->data,
	               NULL,
	               mat_copy->row,
	               NULL,
	               mat_copy->column,
	               superb->data
	              );

	int r;
	for(r = 0; r < S->row; r++) {
		if(matrix_at(S, r, 0) > epsilon) {
			rank++;
		}
	}

	matrix_delete(mat_copy);
	matrix_delete(S);
	matrix_delete(superb);

	return rank;
}

void solve_linear_system(matrix_t *A, matrix_t *X, matrix_t *B)
{
	memcpy(X->data, B->data, sizeof(double) * X->row * X->column);

	int *pivots = (int *)malloc(sizeof(int) * A->row);
	LAPACKE_dgesv(LAPACK_ROW_MAJOR, A->row, X->column,
	              A->data, A->column, pivots, X->data, X->column);
	free(pivots);
}

void matrix_qr_factorization(matrix_t *A, matrix_t **Q_ret, matrix_t **R_ret)
{
	//check: https://www.netlib.org/lapack/lug/node40.html

	matrix_t *Q = matrix_new(A->row, A->row);    //m-by-m
	matrix_t *R = matrix_new(A->row, A->column); //m-by-n

	int rank;
	if(A->row < A->column) {
		rank = A->row;
	} else {
		rank = A->column;
	}

	int m = A->row, n = A->column;
	double *tau = (double *)malloc(sizeof(double) * rank);

	matrix_copy(R, A);
	LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, m, n, R->data, n, tau);

	matrix_t *Q_packed = matrix_new(A->row, A->column);
	matrix_copy(Q_packed, R);
	LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m, rank, rank, Q_packed->data, n, tau);

	//Q matrix
	int r, c;
	for(r = 0; r < Q->row; r++) {
		for(c = 0; c < Q->column; c++) {
			matrix_at(Q, r, c) = matrix_at(Q_packed, r, c);
		}
	}
	//R matrix
	int diag_start = 0;
	for(c = 0; c < n; c++) {
		for(r = diag_start + 1; r < m; r++) {
			matrix_at(R, r, c) = 0;
		}
		diag_start++;
	}

	*Q_ret = Q;
	*R_ret = R;
}

void matrix_construct(matrix_t *mat, int r, int c, double *data)
{
	mat->row = r;
	mat->column = c;
	mat->data = data;
}

matrix_t* matrix_new(int r, int c)
{
	matrix_t *mat = (matrix_t *)malloc(sizeof(matrix_t));

	mat->row = r;
	mat->column = c;
	mat->data = (double *)malloc(sizeof(double) * r * c);

	return mat;
}

matrix_t* matrix_zeros(int r, int c)
{
	matrix_t *mat = (matrix_t *)malloc(sizeof(matrix_t));

	mat->row = r;
	mat->column = c;
	mat->data = (double *)calloc(r * c, sizeof(double));

	return mat;
}

void matrix_delete(matrix_t *mat)
{
	if(mat->data != NULL) free(mat->data);
	if(mat != NULL) free(mat);
}

void matrix_reset_zeros(matrix_t *mat)
{
	size_t size = mat->row * mat->column * sizeof(double);
	memset(mat->data, 0, size);
}

void matrix_inverse(matrix_t *mat, matrix_t *mat_inv)
{
	memcpy(mat_inv->data, mat->data, sizeof(double) * mat_inv->row * mat_inv->column);

	int *pivots = (int *)malloc(sizeof(int) * mat->row);

	//solve matrix inversion by LU decomposition
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, mat_inv->row, mat_inv->column, mat_inv->data,
	               mat_inv->column, pivots);

	LAPACKE_dgetri(LAPACK_ROW_MAJOR, mat_inv->row,  mat_inv->data,
	               mat_inv->column, pivots);
	free(pivots);
}

void matrix_copy(matrix_t *dest, matrix_t *src)
{
	size_t size = dest->row * dest->column * sizeof(double);
	memcpy(dest->data, src->data, size);
}

void matrix_add(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result)
{
	int r, c;
	for(r = 0; r < mat1->row; r++) {
		for(c = 0; c < mat1->column; c++) {
			matrix_at(mat_result, r, c) =
			    matrix_at(mat1, r, c) + matrix_at(mat2, r, c);
		}
	}
}

void matrix_add_by(matrix_t *lhs, matrix_t *rhs)
{
	int r, c;
	for(r = 0; r < lhs->row; r++) {
		for(c = 0; c < lhs->column; c++) {
			matrix_at(lhs, r, c) += matrix_at(rhs, r, c);
		}
	}
}

void matrix_sub(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result)
{
	int r, c;
	for(r = 0; r < mat1->row; r++) {
		for(c = 0; c < mat1->column; c++) {
			matrix_at(mat_result, r, c) =
			    matrix_at(mat1, r, c) - matrix_at(mat2, r, c);
		}
	}
}

void matrix_multiply(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mat_result->row,
	            mat_result->column, mat2->row, 1, mat1->data, mat2->row,
	            mat2->data, mat_result->column, 0, mat_result->data,
	            mat_result->column);
}

void matrix_transpose(matrix_t *mat, matrix_t *trans_mat)
{
	trans_mat->row = mat->column;
	trans_mat->column = mat->row;

	int r, c;
	for(r = 0; r < mat->row; r++) {
		for(c = 0; c < mat->column; c++) {
			matrix_at(trans_mat, c, r) = matrix_at(mat, r, c);
		}
	}
}

void matrix_scaling(double scaler, matrix_t *in, matrix_t *out)
{
	int r, c;
	for(r = 0; r < out->row; r++) {
		for(c = 0; c < out->column; c++) {
			matrix_at(out, r, c) = scaler * matrix_at(in, r, c);
		}
	}
}

void matrix_scale_by(double scaler, matrix_t *mat)
{
	int r, c;
	for(r = 0; r < mat->row; r++) {
		for(c = 0; c < mat->column; c++) {
			matrix_at(mat, r, c) *= scaler;
		}
	}
}

double vector_residual(vector_t *vec1, vector_t *vec2)
{
	double sum_of_squared = 0;
	double diff;

	int r;
	for(r = 0; r < vec1->row; r++) {
		diff = matrix_at(vec1, r, 0) - matrix_at(vec2, r, 0);
		sum_of_squared += diff * diff;
	}

	return sqrt(sum_of_squared);
}

double vector_norm(vector_t *vec)
{
	double sum_of_squared = 0;

	int r;
	for(r = 0; r < vec->row; r++) {
		sum_of_squared += (matrix_at(vec, r, 0) * matrix_at(vec, r, 0));
	}

	return sqrt(sum_of_squared);
}

void print_matrix(char *prompt, matrix_t *mat)
{
	printf("%s (%dx%d) = \n", prompt, mat->row, mat->column);

	int r, c;
	for(r = 0; r < mat->row; r++) {
		printf("   ");
		for(c = 0; c < mat->column; c++) {
			printf("%f  ", mat->data[r * mat->column + c]);
		}
		printf("\n");
	}
}

void print_var(char *prompt, double var)
{
	printf("%s (1x1) = \n   %f\n", prompt, var);
}
