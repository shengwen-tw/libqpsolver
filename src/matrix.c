#include <mkl.h>
#include <mkl_lapacke.h>
#include "matrix.h"

void solver_linear_system(matrix_t *A, matrix_t *X, matrix_t *B, int *pivots)
{
	/* in sgesv, the B matrix serves the purpose of both input and output.
	 * since we want to store the data in X matrix, therefore do the copy
	 * first than pass X as argument */
	*X = *B;

	LAPACKE_sgesv(LAPACK_ROW_MAJOR, A->row, X->column,
		      A->data, A->column, pivots, X->data, X->column);
}

void matrix_inverse(matrix_t *mat, matrix_t *mat_inv, int *pivots)
{
	*mat_inv = *mat;

	LAPACKE_sgetrf(LAPACK_ROW_MAJOR, mat_inv->row, mat_inv->column, mat_inv->data,
		       mat_inv->column, pivots);

	LAPACKE_sgetri(LAPACK_ROW_MAJOR, mat_inv->row,  mat_inv->data,
		       mat_inv->column, pivots);
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
