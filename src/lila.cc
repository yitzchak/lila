#include <clasp/clasp.h>
#include <clasp/core/translators.h>
#include <lila/lila.h>

namespace translate {

template <class T> struct to_object<std::complex<T>> {
  static core::T_sp convert(const std::complex<T> &v) { return core::Complex_O::create(v.real(), v.imag()); }
};

template <> struct from_object<std::complex<float>, std::true_type> {
  typedef std::complex<float> DeclareType;
  DeclareType _v;

  from_object(core::T_sp o) {
    if (gctools::IsA<core::Complex_sp>(o)) {
      _v = std::complex<double>(core::clasp_to_float(gctools::As<core::Complex_sp>(o)->real()),
                                core::clasp_to_float(gctools::As<core::Complex_sp>(o)->imaginary()));
    } else {
      _v = core::clasp_to_float(o);
    }
  }
};

template <> struct from_object<std::complex<double>, std::true_type> {
  typedef std::complex<double> DeclareType;
  DeclareType _v;

  from_object(core::T_sp o) {
    if (gctools::IsA<core::Complex_sp>(o)) {
      _v = std::complex<double>(core::clasp_to_double(gctools::As<core::Complex_sp>(o)->real()),
                                core::clasp_to_double(gctools::As<core::Complex_sp>(o)->imaginary()));
    } else {
      _v = core::clasp_to_double(o);
    }
  }
};

} // namespace translate

PACKAGE_USE("COMMON-LISP");
NAMESPACE_PACKAGE_ASSOCIATION(lila, lila_pkg, "LILA");

namespace lila {

void vector_coerce(core::T_sp obj, bool &double_vec, bool &complex_vec, std::size_t &dimension) {
  core::WrappedPointer_sp p = gctools::As<core::WrappedPointer_sp>(obj);

  r1v *r1v_p = p->castOrNull<r1v>();
  if (r1v_p) {
    dimension = std::max(dimension, r1v_p->dimension());
    return;
  }

  r2v *r2v_p = p->castOrNull<r2v>();
  if (r2v_p) {
    dimension = std::max(dimension, r2v_p->dimension());
    double_vec = true;
    return;
  }

  c1v *c1v_p = p->castOrNull<c1v>();
  if (c1v_p) {
    dimension = std::max(dimension, c1v_p->dimension());
    complex_vec = true;
    return;
  }

  c2v *c2v_p = p->castOrNull<c2v>();
  if (c2v_p) {
    dimension = std::max(dimension, c2v_p->dimension());
    double_vec = true;
    complex_vec = true;
    return;
  }
}

void vector_coerce(core::Vaslist_sp args, bool &double_vec, bool &complex_vec, std::size_t &dimension) {
  for (int i = args->remaining_nargs() - 1; i > -1; --i) {
    vector_coerce(args->next_arg_indexed(i), double_vec, complex_vec, dimension);
  }
}

CL_EXPOSE void lila_startup() {
  clbind::package_ pkg(lila_pkg);

  clbind::class_<r1v>(pkg, "real-single-vector")
      .def_constructor("make-real-single-vector", clbind::constructor<unsigned long, r1>())
      .def("dimension@r1v", &r1v::dimension, clbind::noAutoExport())
      .def("l1-norm@r1v", &r1v::l1_norm, clbind::noAutoExport())
      .def("l2-norm@r1v", &r1v::l2_norm, clbind::noAutoExport())
      .def("l2-norm-sqr@r1v", &r1v::l2_norm_sqr, clbind::noAutoExport());
  clbind::class_<r2v>(pkg, "real-double-vector")
      .def_constructor("make-real-double-vector", clbind::constructor<unsigned long, r2>())
      .def("dimension@r2v", &r2v::dimension, clbind::noAutoExport())
      .def("l1-norm@r2v", &r2v::l1_norm, clbind::noAutoExport())
      .def("l2-norm@r2v", &r2v::l2_norm, clbind::noAutoExport())
      .def("l2-norm-sqr@r2v", &r2v::l2_norm_sqr, clbind::noAutoExport());
  clbind::class_<c1v>(pkg, "complex-single-vector")
      .def_constructor("make-complex-single-vector", clbind::constructor<unsigned long, c1>())
      .def("dimension@c1v", &c1v::dimension, clbind::noAutoExport())
      .def("l1-norm@c1v", &c1v::l1_norm, clbind::noAutoExport())
      .def("l2-norm@c1v", &c1v::l2_norm, clbind::noAutoExport())
      .def("l2-norm-sqr@c1v", &c1v::l2_norm_sqr, clbind::noAutoExport());
  clbind::class_<c2v>(pkg, "complex-double-vector")
      .def_constructor("make-complex-double-vector", clbind::constructor<unsigned long, c2>())
      .def("dimension@c2v", &c2v::dimension, clbind::noAutoExport())
      .def("l1-norm@c2v", &c2v::l1_norm, clbind::noAutoExport())
      .def("l2-norm@c2v", &c2v::l2_norm, clbind::noAutoExport())
      .def("l2-norm-sqr@c2v", &c2v::l2_norm_sqr, clbind::noAutoExport());

  pkg.def(
      "real-single-vector",
      +[](core::Vaslist_sp args) {
        r1v res = r1v(args->total_nargs());
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res[i] = core::clasp_to_float(args->next_arg());
        }
        return res;
      },
      "(core:&va-rest args)"_ll);
  pkg.def(
      "real-double-vector",
      +[](core::Vaslist_sp args) {
        r2v res = r2v(args->total_nargs());
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res[i] = core::clasp_to_double(args->next_arg());
        }
        return res;
      },
      "(core:&va-rest args)"_ll);
  pkg.def(
      "complex-single-vector",
      +[](core::Vaslist_sp args) {
        c1v res = c1v(args->total_nargs());
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res[i] = translate::from_object<c1>(args->next_arg())._v;
        }
        return res;
      },
      "(core:&va-rest args)"_ll);
  pkg.def(
      "complex-double-vector",
      +[](core::Vaslist_sp args) {
        c2v res = c2v(args->total_nargs());
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res[i] = translate::from_object<c2>(args->next_arg())._v;
        }
        return res;
      },
      "(core:&va-rest args)"_ll);

  pkg.def(
      "vref@r1v", +[](const r1v &x, size_t i) { return x[i]; }, clbind::noAutoExport());
  pkg.def(
      "vref@r2v", +[](const r2v &x, size_t i) { return x[i]; }, clbind::noAutoExport());
  pkg.def(
      "vref@c1v", +[](const c1v &x, size_t i) { return x[i]; }, clbind::noAutoExport());
  pkg.def(
      "vref@c2v", +[](const c2v &x, size_t i) { return x[i]; }, clbind::noAutoExport());

  pkg.def(
      "setf-vref@r1v", +[](float value, r1v &x, size_t i) { return x[i] = value; }, clbind::noAutoExport());
  pkg.def(
      "setf-vref@r2v", +[](double value, r2v &x, size_t i) { return x[i] = value; }, clbind::noAutoExport());
  pkg.def(
      "setf-vref@c1v", +[](float value, c1v &x, size_t i) { return x[i] = value; }, clbind::noAutoExport());
  pkg.def(
      "setf-vref@c2v", +[](double value, c2v &x, size_t i) { return x[i] = value; }, clbind::noAutoExport());

  pkg.def(
      "v+",
      +[](core::Vaslist_sp args) {
        bool double_vec = false, complex_vec = false;
        std::size_t dimension = 0;
        vector_coerce(args, double_vec, complex_vec, dimension);
        if (double_vec && complex_vec) {
          c2v res = c2v(dimension);
          for (size_t i = 0; args->remaining_nargs() > 0; i++) {
            core::WrappedPointer_sp wp_sp = gctools::As<core::WrappedPointer_sp>(args->next_arg());
            c2v *c2v_p = wp_sp->castOrNull<c2v>();
            c1v *c1v_p = wp_sp->castOrNull<c1v>();
            r2v *r2v_p = wp_sp->castOrNull<r2v>();
            if (c2v_p) {
              res += *c2v_p;
            } else if (c1v_p) {
              res += *c1v_p;
            } else if (r2v_p) {
              res += *r2v_p;
            } else {
              res += *wp_sp->cast<r1v>();
            }
          }
          return translate::to_object<c2v>::convert(res);
        }
        if (double_vec) {
          r2v res = r2v(dimension);
          for (size_t i = 0; args->remaining_nargs() > 0; i++) {
            core::WrappedPointer_sp wp_sp = gctools::As<core::WrappedPointer_sp>(args->next_arg());
            r2v *r2v_p = wp_sp->castOrNull<r2v>();
            if (r2v_p) {
              res += *r2v_p;
            } else {
              res += *wp_sp->cast<r1v>();
            }
          }
          return translate::to_object<r2v>::convert(res);
        }
        if (complex_vec) {
          c1v res = c1v(dimension);
          for (size_t i = 0; args->remaining_nargs() > 0; i++) {
            core::WrappedPointer_sp wp_sp = gctools::As<core::WrappedPointer_sp>(args->next_arg());
            c1v *c1v_p = wp_sp->castOrNull<c1v>();
            if (c1v_p) {
              res += *c1v_p;
            } else {
              res += *wp_sp->cast<r1v>();
            }
          }
          return translate::to_object<c1v>::convert(res);
        }
        r1v res = r1v(dimension);
        for (size_t i = 0; args->remaining_nargs() > 0; i++) {
          res += *gctools::As<core::WrappedPointer_sp>(args->next_arg())->cast<r1v>();
        }
        return translate::to_object<r1v>::convert(res);
      },
      "(core:&va-rest args)"_ll);

  clbind::class_<r1m>(pkg, "matrix-float")
      .def_constructor("make-matrix-float", clbind::constructor<unsigned long, unsigned long, float>())
      .def("column-dimension", &r1m::column_dimension)
      .def("row-dimension", &r1m::row_dimension)
      .def("transpose", &r1m::transpose);
  clbind::class_<r2m>(pkg, "matrix-double")
      .def_constructor("make-matrix-double", clbind::constructor<unsigned long, unsigned long, double>())
      .def("column-dimension", &r2m::column_dimension)
      .def("row-dimension", &r2m::row_dimension)
      .def("transpose", &r2m::transpose);

  pkg.def(
      "mref@r1m", +[](const r1m &x, size_t i, size_t j) { return x(i, j); }, clbind::noAutoExport());
  pkg.def(
      "mref@r2m", +[](const r2m &x, size_t i, size_t j) { return x(i, j); }, clbind::noAutoExport());

  pkg.def(
      "setf-mref@r1m", +[](float value, r1m &x, size_t i, size_t j) { return x(i, j) = value; }, clbind::noAutoExport());
  pkg.def(
      "setf-mref@r2m", +[](double value, r2m &x, size_t i, size_t j) { return x(i, j) = value; }, clbind::noAutoExport());
}

} // namespace lila
