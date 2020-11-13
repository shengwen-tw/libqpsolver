-include config.mk

EXECUTABLE=qp_solver

CFLAGS=-O3 -g -Wall -Wno-unused-label
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
