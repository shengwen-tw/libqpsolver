#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "libqpsolver.h"
#include "qpsolver.h"

using namespace std;
namespace py = pybind11;

using pyarray = std::optional<py::array_t<FLOAT>>;

void print_py_array(string prompt, pyarray &py_arr)
{
    /* use C++17 std::optional to identify if python object exists (!= None) */
    if(py_arr.has_value() == true) {
        printf("%s (%dx%d) =\n",
               prompt.c_str(),
               (int)py_arr.value().shape(0),
               (int)py_arr.value().shape(1));

        auto data = py_arr.value().unchecked<2>();
        for (py::ssize_t r = 0; r < data.shape(0); r++) {
            printf("   ");
            for (py::ssize_t c = 0; c < data.shape(1); c++) {
                printf("%f ", data(r, c));
            }
            printf("\n");
        }
    }
}

void my_quadprog(pyarray P,    pyarray q,
                 pyarray A,    pyarray b,
                 pyarray A_eq, pyarray b_eq,
                 pyarray lb,   pyarray ub)
{
    print_py_array("P", P);
    print_py_array("q", q);
    print_py_array("A", A);
    print_py_array("b", b);
    print_py_array("A_eq", A_eq);
    print_py_array("b_eq", b_eq);
    print_py_array("lb", lb);
    print_py_array("ub", ub);

    qp_t qp;
    qp_set_default(&qp);
}

PYBIND11_MODULE(libqpsolver, qp)
{
	qp.doc() = "libqpsolver, a quadratic programming library written in C";
	qp.def("my_quadprog", &my_quadprog, "");
}
