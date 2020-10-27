#include <string.h>
#include <math.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "libqpsolver.h"
#include "matrix.h"

void solve_linear_system(matrix_t *A, matrix_t *X, matrix_t *B)
{
	memcpy(X->data, B->data, sizeof(FLOAT) * X->row * X->column);

	int *pivots = (int *)malloc(sizeof(int) * A->row);
	GESV(LAPACK_ROW_MAJOR, A->row, X->column,
	     A->data, A->column, pivots, X->data, X->column);
	free(pivots);
}

void matrix_qr_factorization(matrix_t *A, matrix_t *Q, matrix_t *R)
{
	//check: https://www.netlib.org/lapack/lug/node40.html

	matrix_t *Q_packed = matrix_new(A->row, A->column);

	int rank;
	if(A->row < A->column) {
		rank = A->row;
	} else {
		rank = A->column;
	}

	int m = A->row, n = A->column;
	FLOAT *tau = (FLOAT *)malloc(sizeof(FLOAT) * rank);

	matrix_copy(R, A);
	GEQRF(LAPACK_ROW_MAJOR, m, n, R->data, n, tau);

	matrix_copy(Q_packed, R);
	ORGQR(LAPACK_ROW_MAJOR, m, rank, rank, Q_packed->data, n, tau);

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
}

void matrix_construct(matrix_t *mat, int r, int c, FLOAT *data)
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
	mat->data = (FLOAT *)malloc(sizeof(FLOAT) * r * c);

	return mat;
}

matrix_t* matrix_zeros(int r, int c)
{
	matrix_t *mat = (matrix_t *)malloc(sizeof(matrix_t));

	mat->row = r;
	mat->column = c;
	mat->data = (FLOAT *)calloc(r * c, sizeof(FLOAT));

	return mat;
}

void matrix_delete(matrix_t *mat)
{
	if(mat->data != NULL) free(mat->data);
}

void matrix_reset_zeros(matrix_t *mat)
{
	int r, c;
	for(r = 0; r < mat->row; r++) {
		for(c = 0; c < mat->column; c++) {
			matrix_at(mat, r, c) = 0;
		}
	}
}

void matrix_inverse(matrix_t *mat, matrix_t *mat_inv)
{
	memcpy(mat_inv->data, mat->data, sizeof(FLOAT) * mat_inv->row * mat_inv->column);

	int *pivots = (int *)malloc(sizeof(int) * mat->row);

	//solve matrix inversion by LU decomposition
	GETRF(LAPACK_ROW_MAJOR, mat_inv->row, mat_inv->column, mat_inv->data,
	      mat_inv->column, pivots);

	GETRI(LAPACK_ROW_MAJOR, mat_inv->row,  mat_inv->data,
	      mat_inv->column, pivots);
	free(pivots);
}

void matrix_copy(matrix_t *dest, matrix_t *src)
{
	int r, c;
	for(r = 0; r < dest->row; r++) {
		for(c = 0; c < dest->column; c++) {
			matrix_at(dest, r, c) = matrix_at(src, r, c);
		}
	}
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
	GEMM(CblasRowMajor, CblasNoTrans, CblasNoTrans, mat_result->row,
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

void matrix_scaling(float scaler, matrix_t *mat)
{
	int r, c;
	for(r = 0; r < mat->row; r++) {
		for(c = 0; c < mat->column; c++) {
			matrix_at(mat, r, c) *= scaler;
		}
	}
}

float vector_residual(vector_t *vec1, vector_t *vec2)
{
	float sum_of_squared = 0;
	float diff;

	int r;
	for(r = 0; r < vec1->row; r++) {
		diff = matrix_at(vec1, r, 0) - matrix_at(vec2, r, 0);
		sum_of_squared += diff * diff;
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
