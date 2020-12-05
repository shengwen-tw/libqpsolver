CC=gcc

MKL_PATH=/opt/intel/mkl
MKL_LDFLAGS=$(MKL_PATH)/lib/intel64/libmkl_intel_lp64.so \
	$(MKL_PATH)/lib/intel64/libmkl_sequential.so \
	$(MKL_PATH)/lib/intel64/libmkl_core.so

CFLAGS=
CFLAGS+=-O3 -Wall -g

LDFLAGS=-lm

CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I../../include

SRC=../../src/matrix.c

test_lapack: test_lapack.c $(SRC)
	$(CC) $(CFLAGS) $^ $(MKL_LDFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf test_lapack

.PHONY: all clean
