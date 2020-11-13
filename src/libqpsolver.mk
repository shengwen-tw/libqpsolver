-include ../config.mk

LIB_QPSOLVER=libqpsolver.a
LIB_QPSOLVER_PATH=../include

CFLAGS=-O3 -g -Wall -Wno-unused-label

CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I$(LIB_QPSOLVER_PATH)

SRC=qpsolver.c \
	matrix.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

$(LIB_QPSOLVER): $(OBJS)
	@echo "AR" $@
	@$(AR) -rcs $@ $(OBJS)

-include $(DEPEND)

%.o: %.c
	@echo "CC" $@
	@$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(LIB_QPSOLVER)
	rm -rf $(OBJS)
	rm -rf $(DEPEND)
