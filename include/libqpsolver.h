#ifndef __LIBQPSOLVER_H__
#define __LIBQPSOLVER_H__

/*======================================*
 * enable internal debug message output *
 *======================================*/
#define VERBOSE_MESSAGE 1
#define DEBUG_MESSAGE   0

/*====================================*
 * floating-point data type precision *
 *====================================*/
#define USE_FLOAT  0
#define USE_DOUBLE 1
#define FLOAT_PRECISION USE_DOUBLE

/*===========================================*
 * enable infeasible start (phase1) function *
 *===========================================*/
#define ENABLE_INFEASIBLE_START 1

/*======================================================*
 * enable matrix/vector dimension check (for debugging) *
 *======================================================*/
#define ENABLE_DIMENSION_CHECK  0

#endif
