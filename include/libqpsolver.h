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

/*=======================================================*
 * select lapack/cblas helper function according to the  *
 * precision of the floating-point data type             *
 *=======================================================*/
#if (FLOAT_PRECISION == USE_FLOAT)
#define FLOAT float
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define FLOAT double
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define GESV LAPACKE_sgesv
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define GESV LAPACKE_dgesv
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define GETRF LAPACKE_sgetrf
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define GETRF LAPACKE_dgetrf
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define GETRI LAPACKE_sgetri
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define GETRI LAPACKE_dgetri
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define GEMM cblas_sgemm
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define GEMM cblas_dgemm
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define GEQRF LAPACKE_sgeqrf
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define GEQRF LAPACKE_dgeqrf
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define ORGQR LAPACKE_sorgqr
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define ORGQR LAPACKE_dorgqr
#endif

#if (FLOAT_PRECISION == USE_FLOAT)
#define GESVD LAPACKE_sgesvd
#elif (FLOAT_PRECISION == USE_DOUBLE)
#define GESVD LAPACKE_dgesvd
#endif

#endif
