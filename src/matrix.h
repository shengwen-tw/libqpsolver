#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>
#include "libqpsolver.h"

/* matrix initialization */
#define MATRIX_ARGS(...) __VA_ARGS__
#define DECLARE_MATRIX(name, _row, _column, ...) \
	FLOAT name##_data[_row * _column] = {MATRIX_ARGS __VA_ARGS__}; \
	matrix_t name = { \
		.data = name##_data, \
		.row = _row, \
		.column = _column \
	}

/* malloc new matrix */
#define MALLOC_MATRIX(name, _row, _column) \
	FLOAT *name##_data = (FLOAT *)malloc(sizeof(FLOAT) * _row * _column); \
	matrix_t name = { \
		.data = name##_data, \
		.row = _row, \
		.column = _column \
	};

/* calloc new matrix */
#define CALLOC_MATRIX(name, _row, _column) \
	FLOAT *name##_data = (FLOAT *)calloc(_row * _column, sizeof(FLOAT)); \
	matrix_t name = { \
		.data = name##_data, \
		.row = _row, \
		.column = _column \
	};

#define DELETE_MATRIX(name) \
	free(name.data)

/* access matrix element with row and column index */
#define MATRIX_DATA(mat_ptr, r, c) \
	(mat_ptr)->data[(r * (mat_ptr)->column) + c]

/* matrix debug print wrapper*/
#define PRINT_MATRIX(mat) \
	print_matrix(#mat, &mat)

#define DECLARE_VECTOR DECLARE_MATRIX
#define MALLOC_VECTOR MALLOC_MATRIX
#define CALLOC_VECTOR CALLOC_MATRIX
#define DELETE_VECTOR DELETE_MATRIX
#define VECTOR_DATA MATRIX_DATA
#define PRINT_VECTOR PRINT_MATRIX

typedef struct {
	FLOAT *data;
	int row;
	int column;
} matrix_t;

typedef matrix_t vector_t;

void solve_linear_system(matrix_t *A, matrix_t *X, matrix_t *B);

matrix_t* matrix_new(int r, int c);
matrix_t* matrix_zeros(int r, int c);
void matrix_delete(matrix_t *mat);
void matrix_inverse(matrix_t *mat, matrix_t *mat_inv);
void matrix_copy(matrix_t *dest, matrix_t *src);
void matrix_add(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result);
void matrix_multiply(matrix_t *mat1, matrix_t *mat2, matrix_t *mat_result);
void matrix_scaling(float scaler, matrix_t *mat);
void matrix_transpose(matrix_t *mat, matrix_t *trans_mat);

void vector_scaling(float scaler, vector_t *vec);
void vector_copy(vector_t *dest, vector_t *src);
void vector_negate(vector_t *vec);
float vector_residual(vector_t *vec1, vector_t *vec2);

void print_matrix(char *prompt, matrix_t *mat);

#endif
