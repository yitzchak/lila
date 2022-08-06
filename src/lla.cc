#include <clasp/clasp.h>
#include <clasp/core/translators.h>
#include <lla/lla.h>

PACKAGE_USE("COMMON-LISP");
NAMESPACE_PACKAGE_ASSOCIATION(lla, lla_pkg, "LLA");

namespace lla {

typedef vector<float> v1;
typedef vector<double> v2;

CL_EXPOSE void lla_startup() {
  clbind::package_ pkg(lla_pkg);

  clbind::class_<v1>(pkg, "vector-float")
      .def_constructor("make-vector-float", clbind::constructor<unsigned long, bool>())
      .def("dimension", &v1::dimension)
      .def("rowp", &v1::row)
      .def("transpose", &v1::transpose);
  clbind::class_<v2>(pkg, "vector-double")
      .def_constructor("make-vector-double", clbind::constructor<unsigned long, bool>())
      .def("dimension", &v2::dimension)
      .def("rowp", &v2::row)
      .def("transpose", &v2::transpose);

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
}

} // namespace lla
