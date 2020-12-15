#include <iostream>
#include <sys/time.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "libqpsolver.h"
#include "qpsolver.h"

using namespace std;
namespace py = pybind11;

using pyarray = std::optional<py::array_t<double>>;

bool unit_test_debug_print = false;

double time(void)
{
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, NULL);

	if (sec < 0) sec = tv.tv_sec;
	return (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

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

py::array_t<double> convert_np_array_to_matrix(matrix_t *mat)
{
	unsigned long n_row = mat->row;
	unsigned long n_column = mat->column;

	auto buf_info =
	        py::buffer_info(mat->data,
	                        sizeof(double),
	                        py::format_descriptor<double>::format(),
	                        2,
	                        std::vector<size_t> {n_row, n_column},
	                        std::vector<size_t> {n_column * sizeof(double), sizeof(double)}
	                       );
	return py::array_t<double>(buf_info);
}

py::array_t<double> my_quadprog(pyarray P_numpy,    pyarray q_numpy,
                                pyarray A_numpy,    pyarray b_numpy,
                                pyarray A_eq_numpy, pyarray b_eq_numpy,
                                pyarray lb_numpy,   pyarray ub_numpy)
{
	matrix_t *x = NULL;
	matrix_t *P = NULL;
	matrix_t *q = NULL;
	matrix_t *A = NULL;
	matrix_t *b = NULL;
	matrix_t *A_eq = NULL;
	matrix_t *b_eq = NULL;
	matrix_t *lb = NULL;
	matrix_t *ub = NULL;

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
			qp_solve_set_cost_function(&qp, P, q);
			test_print_matrix("q", q);
		} else {
			qp_solve_set_cost_function(&qp, P, NULL);
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

	double start_time = time();
	int qp_ret = qp_solve_start(&qp);
	double end_time = time();

	if(unit_test_debug_print == true) {
		if(qp_ret == QP_SUCCESS_SOLVED) {
			printf("the optimal solution of the problem is:\n");
			test_print_matrix("x", x);
		} else {
			printf("unable to solve the problem!\n");
		}
	}

	if(unit_test_debug_print == true) {
		printf("run time: %lf seconds\n"
		       "phase1 stage took %d iterations\n"
		       "phase2 stage took %d iterations\n",
		       end_time - start_time, qp.phase1.iters, qp.phase2.iters + 1);
	}

	return convert_np_array_to_matrix(x);
}

PYBIND11_MODULE(libqpsolver, mod)
{
	mod.doc() = "libqpsolver, a quadratic programming library written in C";
	mod.def("quadprog", &my_quadprog, "");
}
