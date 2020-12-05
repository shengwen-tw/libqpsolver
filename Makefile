-include config.mk

example:
	#build example program of libqpsolver
	@echo "Build example program"
	@$(MAKE) -C ./example -f build_example.mk

libqpsolver:
	#build library as static object
	@echo "Build libqpsolver.a"
	@$(MAKE) -C ./src -f libqpsolver.mk

unit_test_qp:
	#test of linear algebra functions
	@$(MAKE) -C ./unit_test/qpsolver -f Makefile
	@./unit_test/qpsolver/run_test.sh

unit_test_matrix:
	#test of quadratic programming solver
	@$(MAKE) -C ./unit_test/matrix -f Makefile
	@./unit_test/matrix/run_test.sh

astyle:
	astyle --style=linux --indent=tab=4 --recursive "*.c,*.h"

clean:
	@$(MAKE) -C ./src -f libqpsolver.mk clean
	@$(MAKE) -C ./example -f build_example.mk clean
	@$(MAKE) -C ./unit_test/qpsolver -f Makefile clean

.PHONY: example libqpsolver unit_test astyle clean
