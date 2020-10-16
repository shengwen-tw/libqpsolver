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

/* matrix debug print wrapper*/
#define PRINT_MATRIX(mat) \
	print_matrix(#mat, &mat)

typedef struct {
	FLOAT *data;
	int row;
	int column;
} matrix_t;

typedef struct {
	FLOAT *data;
	int size;
} vector_t;

void solver_linear_system(matrix_t *A, matrix_t *X, matrix_t *B, int *pivot);

void print_matrix(char *prompt, matrix_t *mat);

#endif
