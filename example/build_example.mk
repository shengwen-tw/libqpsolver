-include ../config.mk

EXECUTABLE=qp_solver

LIB_QPSOLVER_INC=../include
LIB_QPSOLVER=../src/libqpsolver.a

CFLAGS=-O2 -g -Wall -Wno-unused-label

CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I$(LIB_QPSOLVER_INC)

LDFLAGS=-L../src -lqpsolver $(MKL_LDFLAGS) -lm

SRC=qp_example.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

$(EXECUTABLE): $(OBJS) $(LIB_QPSOLVER)
	@echo "LD" $@
	@$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $@

$(LIB_QPSOLVER):
	@echo "Build libqpsolver.a"
	@$(MAKE) -C ../src -f libqpsolver.mk

-include $(DEPEND)
%.o: %.c
	@echo "CC" $@
	@$(CC)  $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(EXECUTABLE)
	rm -rf $(OBJS)
	rm -rf $(DEPEND)

.PHONY: clean
