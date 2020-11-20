-include ../../config.mk

CXX=g++

TARGET=libqpsolver
LIB_QPSOLVER = $(TARGET)$(shell python3-config --extension-suffix)

LIB_QPSOLVER_PATH=../..

CXXFLAGS+=-std=c++17 -O3 -g -Wall -Wno-unused-label -shared -fPIC
CXXFLAGS+=`python3 -m pybind11 --includes`

CXXFLAGS+=-I$(MKL_PATH)/include
CXXFLAGS+=-I$(LIB_QPSOLVER_PATH)/include

LDFLAGS=$(MKL_LDFLAGS) -lm

CSRC=qpsolver.c \
	matrix.c
CXXSRC=pybind.cpp

COBJS=$(CSRC:.c=.o)
CXXOBJS=$(CXXSRC:.cpp=.o)

CDEPEND=$(CSRC:.c=.d)
CXXDEPEND=$(CXXSRC:.cpp=.d)

$(LIB_QPSOLVER): $(COBJS) $(CXXOBJS)
	@echo "LD" $@
	@$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

-include $(CDEPEND)
-include $(CPPDEPEND)

%.o: $(LIB_QPSOLVER_PATH)/src/%.c
	@echo "CXX" $@
	@$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

%.o: %.cpp
	@echo "CXX" $@
	@$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(LIB_QPSOLVER)
	rm -rf $(COBJS)
	rm -rf $(CXXOBJS)
	rm -rf $(CDEPEND)
	rm -rf $(CXXDEPEND)
