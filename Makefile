-include config.mk

example:
	@echo "Build example program"
	@$(MAKE) -C ./example -f build_example.mk

libqpsolver:
	@echo "Build libqpsolver.a"
	@$(MAKE) -C ./src -f libqpsolver.mk

astyle:
	astyle --style=linux --indent=tab=4 --recursive "*.c,*.h"

clean:
	@$(MAKE) -C ./src -f libqpsolver.mk clean
	@$(MAKE) -C ./example -f build_example.mk clean

.PHONY: example libqpsolver all clean astyle
