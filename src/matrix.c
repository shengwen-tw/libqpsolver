#include <string.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "matrix.h"

void solver_linear_system(matrix_t *A, matrix_t *X, matrix_t *B, int *pivots)
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

void print_matrix(char *prompt, matrix_t *mat)
{
	printf("%s = \n", prompt);

	int r, c;
	for(r = 0; r < mat->row; r++) {
		printf("   ");
		for(c = 0; c < mat->column; c++) {
			printf("%f  ", mat->data[r * mat->column + c]);
		}
		printf("\n");
	}
}
