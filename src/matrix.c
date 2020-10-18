#include <string.h>
#include <math.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "matrix.h"

void solve_linear_system(matrix_t *A, matrix_t *X, matrix_t *B, int *pivots)
{
	memcpy(X->data, B->data, sizeof(FLOAT) * X->row * X->column);

	LAPACKE_sgesv(LAPACK_ROW_MAJOR, A->row, X->column,
		      A->data, A->column, pivots, X->data, X->column);
}

void matrix_inverse(matrix_t *mat, matrix_t *mat_inv, int *pivots)
{
	memcpy(mat_inv->data, mat->data, sizeof(FLOAT) * mat_inv->row * mat_inv->column);

	//solve matrix inversion by LU decomposition
	LAPACKE_sgetrf(LAPACK_ROW_MAJOR, mat_inv->row, mat_inv->column, mat_inv->data,
		       mat_inv->column, pivots);

	LAPACKE_sgetri(LAPACK_ROW_MAJOR, mat_inv->row,  mat_inv->data,
		       mat_inv->column, pivots);
}

void matrix_multiply(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, mat_result->row,
		    mat_result->column, mat2->row, 1, mat1->data, mat2->row,
		    mat2->data, mat_result->column, 0, mat_result->data,
		    mat_result->column);
}

void matrix_scaling(float scaler, matrix_t *mat)
{
	int r, c;
	for(r = 0; r < mat->row; r++) {
		for(c = 0; c < mat->column; c++) {
			MATRIX_DATA(mat, r, c) *= scaler;
		}
	}
}

void vector_scaling(float scaler, vector_t *vec)
{
	int r;
	for(r = 0; r < vec->row; r++) {
		MATRIX_DATA(vec, r, 0) *= scaler;
	}
}

void vector_copy(vector_t *dest, vector_t *src)
{
	int r;
	for(r = 0; r < dest->row; r++) {
		MATRIX_DATA(dest, r, 0) = MATRIX_DATA(src, r, 0);
	}
}

void vector_negate(vector_t *vec)
{
	int r;
	for(r = 0; r < vec->row; r++) {
		MATRIX_DATA(vec, r, 0) *= -1;
	}
}

float vector_residual(vector_t *vec1, vector_t *vec2)
{
	float sum_of_squared = 0;
	float diff;

	int r;
	for(r = 0; r < vec1->row; r++) {
		diff = MATRIX_DATA(vec1, r, 0) - MATRIX_DATA(vec2, r, 0);
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
