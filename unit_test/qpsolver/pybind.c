#include <pybind11/pybind11.h>
#include "libqpsolver.h"
#include "qpsolver.h"

PYBIND11_MODULE(libqpsolver, qp)
{
	qp.doc() = "libqpsolver, a quadratic programming library written in C";
	qp.def("solver_set_default", &qp_set_default, "");
	qp.def("set_optimization_variable", &qp_solve_set_optimization_variable, "");
	qp.def("set_cost_function", &qp_solve_set_cost_function, "");
	qp.def("set_equality_constraint", &qp_solve_set_equality_constraints, "");
	qp.def("set_upper_bound", &qp_solve_set_upper_bound_inequality_constraints, "");
	qp.def("set_lower_bound", &qp_solve_set_lower_bound_inequality_constraints, "");
	qp.def("set_inequality_constraint", &qp_solve_set_affine_inequality_constraints, "");
	qp.def("solve", &qp_solve_start, "");
}
