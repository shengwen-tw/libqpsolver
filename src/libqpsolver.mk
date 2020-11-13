-include ../config.mk

CFLAGS=-O3 -g -Wall -Wno-unused-label
CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I../include

LDFLAGS=-lm

SRC=matrix.c \
	qpsolver.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

libqpsolver.a: $(OBJS)
	@echo "AR" $@
	@$(AR) rcs $@ $^

-include $(DEPEND)

$(OBJS): $(SRC)
	@echo "CC" $@
	@$(CC) $(CFLAGS) -MMD -MP -c $< $(LDFLAGS) -o $@
