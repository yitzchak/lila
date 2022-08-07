#include <clasp/clasp.h>
#include <clasp/core/translators.h>
#include <lla/lla.h>

PACKAGE_USE("COMMON-LISP");
NAMESPACE_PACKAGE_ASSOCIATION(lla, lla_pkg, "LLA");

namespace lla {

typedef vector<float> v1;
typedef vector<double> v2;
typedef matrix<float> m1;
typedef matrix<double> m2;

CL_EXPOSE void lla_startup() {
  clbind::package_ pkg(lla_pkg);

  clbind::class_<v1>(pkg, "vector-float")
      .def_constructor("make-vector-float", clbind::constructor<unsigned long, float>())
      .def("dimension", &v1::dimension)
      .def("dot", &v1::dot)
      .def("inner", &v1::inner);
  clbind::class_<v2>(pkg, "vector-double")
      .def_constructor("make-vector-double", clbind::constructor<unsigned long, double>())
      .def("dimension", &v2::dimension)
      .def("dot", &v2::dot)
      .def("inner", &v2::inner);

  pkg.def(
      "vector-float",
      +[](core::Vaslist_sp args) {
        v1 res = v1(args->total_nargs());
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res[i] = core::clasp_to_float(args->next_arg());
        }
        return res;
      },
      "(core:&va-rest args)"_ll);
  pkg.def(
      "vector-double",
      +[](core::Vaslist_sp args) {
        v2 res = v2(args->total_nargs());
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res[i] = core::clasp_to_double(args->next_arg());
        }
        return res;
      },
      "(core:&va-rest args)"_ll);

  pkg.def(
      "vref@v1", +[](const v1 &x, size_t i) { return x[i]; }, clbind::noAutoExport());
  pkg.def(
      "vref@v2", +[](const v2 &x, size_t i) { return x[i]; }, clbind::noAutoExport());

  pkg.def(
      "setf-vref@v1", +[](float value, v1 &x, size_t i) { return x[i] = value; }, clbind::noAutoExport());
  pkg.def(
      "setf-vref@v2", +[](double value, v2 &x, size_t i) { return x[i] = value; }, clbind::noAutoExport());

  clbind::class_<m1>(pkg, "matrix-float")
      .def_constructor("make-matrix-float", clbind::constructor<unsigned long, unsigned long, float>())
      .def("column-dimension", &m1::column_dimension)
      .def("row-dimension", &m1::row_dimension)
      .def("transpose", &m1::transpose);
  clbind::class_<m2>(pkg, "matrix-double")
      .def_constructor("make-matrix-double", clbind::constructor<unsigned long, unsigned long, double>())
      .def("column-dimension", &m2::column_dimension)
      .def("row-dimension", &m2::row_dimension)
      .def("transpose", &m2::transpose);

  pkg.def(
      "mref@m1", +[](const m1 &x, size_t i, size_t j) { return x(i, j); }, clbind::noAutoExport());
  pkg.def(
      "mref@m2", +[](const m2 &x, size_t i, size_t j) { return x(i, j); }, clbind::noAutoExport());

  pkg.def(
      "setf-mref@m1", +[](float value, m1 &x, size_t i, size_t j) { return x(i, j) = value; }, clbind::noAutoExport());
  pkg.def(
      "setf-mref@m2", +[](double value, m2 &x, size_t i, size_t j) { return x(i, j) = value; }, clbind::noAutoExport());
}

} // namespace lla
