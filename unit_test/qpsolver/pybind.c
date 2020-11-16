#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "libqpsolver.h"
#include "qpsolver.h"

using namespace std;
namespace py = pybind11;

using pyarray = std::optional<py::array_t<FLOAT>>;

bool unit_test_debug_print = true;

void test_print_matrix(string prompt, matrix_t *mat)
{
	if(unit_test_debug_print == false) {
		return;
	}

	printf("%s (%dx%d) = \n", prompt.c_str(), mat->row, mat->column);

	int r, c;
	for(r = 0; r < mat->row; r++) {
		printf("   ");
		for(c = 0; c < mat->column; c++) {
			printf("%f  ", mat->data[r * mat->column + c]);
		}
		printf("\n");
	}
}

bool convert_np_array_to_matrix(matrix_t **mat, pyarray &py_arr)
{
	/* use C++17 std::optional to identify if python object exists (!= None) */
	if(py_arr.has_value() == false) {
		return false;
	}

	*mat = matrix_new((int)py_arr.value().shape(0), (int)py_arr.value().shape(1));

	auto data = py_arr.value().unchecked<2>();
	for (py::ssize_t i = 0; i < data.shape(0); i++) {
		for (py::ssize_t j = 0; j < data.shape(1); j++) {
			matrix_at(*mat, i, j) = data(i, j);
		}
	}

	return true;
}

void my_quadprog(pyarray P_numpy,    pyarray q_numpy,
                 pyarray A_numpy,    pyarray b_numpy,
                 pyarray A_eq_numpy, pyarray b_eq_numpy,
                 pyarray lb_numpy,   pyarray ub_numpy)
{
	matrix_t *x, *P, *q, *A, *b, *A_eq, *b_eq, *lb, *ub;

	bool has_P = convert_np_array_to_matrix(&P, P_numpy);
	bool has_q = convert_np_array_to_matrix(&q, q_numpy);
	bool has_A = convert_np_array_to_matrix(&A, A_numpy);
	bool has_b = convert_np_array_to_matrix(&b, b_numpy);
	bool has_A_eq = convert_np_array_to_matrix(&A_eq, A_eq_numpy);
	bool has_b_eq = convert_np_array_to_matrix(&b_eq, b_eq_numpy);
	bool has_lb = convert_np_array_to_matrix(&lb, lb_numpy);
	bool has_ub = convert_np_array_to_matrix(&ub, ub_numpy);

	qp_t qp;
	qp_set_default(&qp);

	if(has_P == true) {
		x = matrix_zeros(P_numpy.value().shape(1), 1);
		qp_solve_set_optimization_variable(&qp, x);
		test_print_matrix("x", x);

		test_print_matrix("P", P);

		if(has_q == true) {
			qp_solve_set_cost_function(&qp, P, q, NULL);
			test_print_matrix("q", q);
		} else {
			qp_solve_set_cost_function(&qp, P, NULL, NULL);
		}
	}

	if((has_A == true) && (has_b == true)) {
		qp_solve_set_affine_inequality_constraints(&qp, A, b);
		test_print_matrix("A", A);
		test_print_matrix("b", b);
	}

	if((has_A_eq == true) && (has_b_eq == true)) {
		qp_solve_set_equality_constraints(&qp, A_eq, b_eq);
		test_print_matrix("A_eq", A_eq);
		test_print_matrix("b_eq", b_eq);
	}

	if(has_lb == true) {
		qp_solve_set_lower_bound_inequality_constraints(&qp, lb);
		test_print_matrix("lb", lb);
	}

	if(has_ub == true) {
		qp_solve_set_upper_bound_inequality_constraints(&qp, ub);
		test_print_matrix("ub", ub);
	}

	qp_solve_start(&qp);
}

PYBIND11_MODULE(libqpsolver, qp)
{
	qp.doc() = "libqpsolver, a quadratic programming library written in C";
	qp.def("my_quadprog", &my_quadprog, "");
}