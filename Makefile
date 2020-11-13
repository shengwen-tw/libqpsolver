-include config.mk

EXECUTABLE=qp_solver

CFLAGS=-O2 -g -Wall -Wno-unused-label
CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I./include

LDFLAGS=-Lsrc -lqpsolver $(MKL_LDFLAGS) -lm

SRC=example/qp_test.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) src/libqpsolver.a
	@echo "LD" $@
	@$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $@

src/libqpsolver.a:
	@echo "Build libqpsolver.a"
	@$(MAKE) -C ./src -f libqpsolver.mk

-include $(DEPEND)
%.o: %.c
	@echo "CC" $@
	@$(CC) -I./include $(CFLAGS) -MMD -MP -c $< -o $@

astyle:
	astyle --style=linux --indent=tab=4 --recursive "*.c,*.h"

clean:
	rm -rf $(EXECUTABLE)
	rm -rf $(OBJS)
	rm -rf $(DEPEND)

.PHONY: all clean astyle
