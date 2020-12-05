# libqpsolver

A quadratic programming solver library written in C

<img src="https://github.com/shengwen-tw/libqpsolver/blob/master/images/qp.png?raw=true" width="18%" height="18%">

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

**Dependencies**

```
pip3 install setuptools numpy pybind11 progressbar

sudo apt install libsuitesparse-dev
git clone https://github.com/cvxopt/cvxopt.git
cd ./cvxopt
sudo python3 setup.py install
```

**Run unit test of the linear algrebra functions**

The program test the basic linear algebra functions underlying of the libqpsolver with randomly generated problems and compare the results with the [Numpy](https://numpy.org/). It is designed to run for several times (~250) and may take minutes to complete.

```
make unit_test_matrix
```

**Run unit test of the quadratic programming solver**

The unit test launchs the solver to solve randomly generated QP problems and compare the results with [CVXOPT](http://cvxopt.org/). It is designed to run for several times (~13000) and may take minutes to complete.
```
make unit_test_qp
```
