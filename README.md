# libqpsolver

A quadratic programming solving library written in C

## Features

1. Solving quadratic programming problem with or without equality, inequality, lower bound and upper bound constraints
2. Primal based interrior point method (log barrier method)
3. Support infeasible start with basic phase1 method

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
