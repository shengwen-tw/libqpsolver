# libqpsolver

A quadratic programming solver library written in C

<img src="https://github.com/shengwen-tw/libqpsolver/blob/master/materials/qp.png?raw=true" width="18%" height="18%">

## Algorithm

The libqpsolver use **Primal Interior Point Method** to solve quadratic programming problems based on Chapter 9 to 11 of [Convex Optimization](https://web.stanford.edu/~boyd/cvxbook/) by Stephen Boyd and Lieven Vandenberghe

## Features

1. Support solving quadratic programming problems with or without equality, inequality, lower bound and upper bound constraints
2. Based on Primal Interrior Point Method
3. Support infeasible start with basic phase1 method (Not yet been fully tested)
4. Based on Intel Math Kernel Libray (Planed to support OpenBLAS in the future)

## Prerequisite

[Intel Math Kernel Library (MKL)](https://software.intel.com/content/www/us/en/develop/tools/performance-libraries.html)

Installation:

```
./mkl_install.sh
```

Make sure the MKL_PATH variable in ``config.mk`` is as same as your installation, the default directory path is set as:

```
MKL_PATH=/opt/intel/mkl
```

## Build and Run

```
git clone https://github.com/shengwen-tw/libqpsolver.git
cd libqpsolver
make example
./example/qp_example
```

## Unit test

The unit test launchs the solver to solve randomly generated QP problems several times (~10000) and compare the results with [CVXOPT](http://cvxopt.org/), it may take several minutes to complete the process.

```
make unit_test
```
