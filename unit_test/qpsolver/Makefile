-include ../../config.mk

CXX=g++

TARGET=libqpsolver
LIB_QPSOLVER = $(TARGET)$(shell python3-config --extension-suffix)

LIB_QPSOLVER_PATH=../../include

CFLAGS+=-std=c++11 -O3 -g -Wall -Wno-unused-label -shared -fPIC
CFLAGS+=`python3 -m pybind11 --includes`

CFLAGS+=-I$(MKL_PATH)/include
CFLAGS+=-I$(LIB_QPSOLVER_PATH)

LDFLAGS=$(MKL_LDFLAGS) -lm -lstdc++

SRC=../../src/qpsolver.c \
	../../src/matrix.c \
	pybind.c

OBJS=$(SRC:.c=.o)
DEPEND=$(SRC:.c=.d)

$(LIB_QPSOLVER): $(OBJS)
	@echo "LD" $@
	@$(CXX) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $@

-include $(DEPEND)

%.o: %.c
	@echo "CXX" $@
	@$(CXX) $(CFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(LIB_QPSOLVER)
	rm -rf $(OBJS)
	rm -rf $(DEPEND)
