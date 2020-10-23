#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>
#include "libqpsolver.h"

#define ELEMENTS(...) (FLOAT []){__VA_ARGS__}

#define matrix_at(mat_ptr, r, c) (mat_ptr)->data[(r * (mat_ptr)->column) + c]
#define vector_at(vec_ptr, r, c) (vec_ptr)->data[(r * (vec_ptr)->column) + c]

#define PRINT_MATRIX(mat) print_matrix(#mat, &mat)
#define PRINT_VECTOR(vec) print_matrix(#vec, &vec)

typedef struct {
	FLOAT *data;
	int row;
	int column;
} matrix_t;

typedef matrix_t vector_t;

void solve_linear_system(matrix_t *A, matrix_t *X, matrix_t *B);

void matrix_construct(matrix_t *mat, int r, int c, FLOAT *data);
matrix_t* matrix_new(int r, int c);
matrix_t* matrix_zeros(int r, int c);
void matrix_delete(matrix_t *mat);
void matrix_reset_zeros(matrix_t *mat);
void matrix_inverse(matrix_t *mat, matrix_t *mat_inv);
void matrix_copy(matrix_t *dest, matrix_t *src);
void matrix_add(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result);
void matrix_multiply(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result);
void matrix_scaling(float scaler, matrix_t *mat);
void matrix_transpose(matrix_t *mat, matrix_t *trans_mat);

void vector_construct(vector_t *vec, int r, int c, FLOAT *data);
vector_t* vector_new(int r, int c);
vector_t* vector_zeros(int r, int c);
void vector_reset_zeros(vector_t *vec);
void vector_delete(vector_t *vec);
void vector_scaling(float scaler, vector_t *vec);
void vector_copy(vector_t *dest, vector_t *src);
void vector_negate(vector_t *vec);
float vector_residual(vector_t *vec1, vector_t *vec2);

void print_matrix(char *prompt, matrix_t *mat);

#endif
