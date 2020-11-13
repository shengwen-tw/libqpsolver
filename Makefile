EXECUTABLE=./example/qp_solver
CC=gcc

MKL_PATH=/opt/intel/mkl
MKL_LDFLAGS=$(MKL_PATH)/lib/intel64/libmkl_intel_lp64.so \
	$(MKL_PATH)/lib/intel64/libmkl_sequential.so \
	$(MKL_PATH)/lib/intel64/libmkl_core.so

CFLAGS=
CFLAGS+=-O3 -g -Wall -Wno-unused-label
CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I./include

LDFLAGS=-lm

SRC=example/qp_test.c \
	src/matrix.c \
	src/qpsolver.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	@echo "LD" $@
	@$(CC) $(CFLAGS) $(OBJS) $(MKL_LDFLAGS) $(LDFLAGS) -o $@

-include $(DEPEND)

%.o: %.c
	@echo "CC" $@
	@$(CC) $(CFLAGS) -MMD -MP -c $< $(LDFLAGS) -o $@

astyle:
	astyle --style=linux --indent=tab=4 --recursive "*.c,*.h"

clean:
	rm -rf $(EXECUTABLE)
	rm -rf $(OBJS)
	rm -rf $(DEPEND)

.PHONY: all clean astyle
