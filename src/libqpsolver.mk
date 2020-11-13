-include ../config.mk

CFLAGS=-O3 -g -Wall -Wno-unused-label
CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I../include

LDFLAGS=-lm

SRC=qpsolver.c \
	matrix.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

libqpsolver.a: $(OBJS)
	@echo "AR" $@
	@$(AR) -rcs $@ $(OBJS)

-include $(DEPEND)

%.o: %.c
	@echo "CC" $@
	@$(CC) $(CFLAGS) -MMD -MP -c $< $(LDFLAGS) -o $@
