#include <stdio.h>
#include <mkl.h>
#include <mkl_lapacke.h>

void print_matrix(float *mat, int r, int c)
{
	int i, j;
	for(i = 0; i < r; i++) {
		for(j = 0; j < c; j++) {
			printf("%f ", mat[i * r + j]);
		}
		printf("\n");
	}
}

int main(void)
{
	/* solve Ax = b */
	float A[3*3] = {3, 5, 2,
			2, 1, 3,
			4, 3, 2};
	int A_size = 3;

	float b[2*3] = {57, 23,
	 		22, 12,
			41, 84};
	int b_col = 3;

	int ipiv[3] = {0};

	printf("A matrix:\n");
	print_matrix(A, 3, 3);

	printf("\nb matrix:\n");
	print_matrix(b, 3, 2);

	LAPACKE_sgesv(LAPACK_COL_MAJOR,
			      A_size,
			      b_col,
			      A,
			      A_size,
			      ipiv,
			      b,
			      b_col);

	printf("\nsolve Ax=b\n");
	print_matrix(b, 3, 2);

	return 0;
}
