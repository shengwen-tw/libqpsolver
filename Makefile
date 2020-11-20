-include config.mk

example:
	@echo "Build example program"
	@$(MAKE) -C ./example -f build_example.mk

libqpsolver:
	@echo "Build libqpsolver.a"
	@$(MAKE) -C ./src -f libqpsolver.mk

unit_test:
	@$(MAKE) -C ./unit_test/qpsolver -f Makefile
	@python3 ./unit_test/qpsolver/qp_test.py

astyle:
	astyle --style=linux --indent=tab=4 --recursive "*.c,*.h"

clean:
	@$(MAKE) -C ./src -f libqpsolver.mk clean
	@$(MAKE) -C ./example -f build_example.mk clean
	@$(MAKE) -C ./unit_test/qpsolver -f Makefile clean

.PHONY: example libqpsolver unit_test astyle clean
