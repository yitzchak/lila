#include <clasp/clasp.h>
#include <clasp/core/translators.h>
#include <lila/lila.h>

namespace lila {

void vector_desc(core::T_sp obj, bool &double_vec, bool &complex_vec, std::size_t &dimension) {
  if (gctools::IsA<RealSingleVector_sp>(obj)) {
    dimension = std::max(dimension, gctools::As<RealSingleVector_sp>(obj)->dimension());
  } else if (gctools::IsA<RealDoubleVector_sp>(obj)) {
    dimension = std::max(dimension, gctools::As<RealDoubleVector_sp>(obj)->dimension());
    double_vec = true;
  } else if (gctools::IsA<ComplexSingleVector_sp>(obj)) {
    dimension = std::max(dimension, gctools::As<ComplexSingleVector_sp>(obj)->dimension());
    complex_vec = true;
  } else if (gctools::IsA<ComplexDoubleVector_sp>(obj)) {
    dimension = std::max(dimension, gctools::As<ComplexDoubleVector_sp>(obj)->dimension());
    double_vec = true;
    complex_vec = true;
  } else if (gctools::IsA<core::Number_sp>(obj)) {
    switch (core::clasp_t_of(gctools::As<core::Number_sp>(obj))) {
    case core::number_DoubleFloat:
    case core::number_LongFloat:
      double_vec = true;
      break;
    case core::number_Complex:
      complex_vec = true;
      break;
    default:
      break;
    }
  }
}

void vector_desc(core::Vaslist_sp args, bool &double_vec, bool &complex_vec, std::size_t &dimension) {
  for (int i = args->remaining_nargs() - 1; i > -1; --i) {
    vector_desc(args->next_arg_indexed(i), double_vec, complex_vec, dimension);
  }
}

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("real-single-vector");
CL_DEFUN RealSingleVector_sp RealSingleVector_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<RealSingleVector_O>::allocate_with_default_constructor();
  res->_Value.resize(args->total_nargs());
  for (size_t i = 0; args->remaining_nargs() > 0; i++) {
    res->_Value[i] = core::clasp_to_float(args->next_arg());
  }
  return res;
}

CL_LAMBDA(vector)
CL_LISPIFY_NAME("dimension")
CL_DEFMETHOD std::size_t RealSingleVector_O::dimension() const { return _Value.dimension(); }

CL_LAMBDA(vector index)
CL_LISPIFY_NAME("vref")
CL_DEFMETHOD float RealSingleVector_O::vref(std::size_t index) const { return _Value[index]; }

CL_LAMBDA(vector new-value index)
CL_LISPIFY_NAME("setf_vref")
CL_DEFMETHOD float RealSingleVector_O::setf_vref(float new_value, std::size_t index) { return _Value[index] = new_value; }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l1_norm")
CL_DEFMETHOD r1 RealSingleVector_O::l1_norm() const { return _Value.l1_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm")
CL_DEFMETHOD r1 RealSingleVector_O::l2_norm() const { return _Value.l2_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm_sqr")
CL_DEFMETHOD r1 RealSingleVector_O::l2_norm_sqr() const { return _Value.l2_norm_sqr(); }

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("real-double-vector");
CL_DEFUN RealDoubleVector_sp RealDoubleVector_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<RealDoubleVector_O>::allocate_with_default_constructor();
  res->_Value.resize(args->total_nargs());
  for (size_t i = 0; args->remaining_nargs() > 0; i++) {
    res->_Value[i] = core::clasp_to_float(args->next_arg());
  }
  return res;
}

CL_LAMBDA(vector)
CL_LISPIFY_NAME("dimension")
CL_DEFMETHOD std::size_t RealDoubleVector_O::dimension() const { return _Value.dimension(); }

CL_LAMBDA(vector index)
CL_LISPIFY_NAME("vref")
CL_DEFMETHOD double RealDoubleVector_O::vref(std::size_t index) const { return _Value[index]; }

CL_LAMBDA(vector new-value index)
CL_LISPIFY_NAME("setf_vref")
CL_DEFMETHOD double RealDoubleVector_O::setf_vref(double new_value, std::size_t index) { return _Value[index] = new_value; }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l1_norm")
CL_DEFMETHOD r2 RealDoubleVector_O::l1_norm() const { return _Value.l1_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm")
CL_DEFMETHOD r2 RealDoubleVector_O::l2_norm() const { return _Value.l2_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm_sqr")
CL_DEFMETHOD r2 RealDoubleVector_O::l2_norm_sqr() const { return _Value.l2_norm_sqr(); }

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("complex-single-vector");
CL_DEFUN ComplexSingleVector_sp ComplexSingleVector_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<ComplexSingleVector_O>::allocate_with_default_constructor();
  res->_Value.resize(args->total_nargs());
  for (size_t i = 0; args->remaining_nargs() > 0; i++) {
    res->_Value[i] = translate::from_object<c1>(args->next_arg())._v;
  }
  return res;
}

CL_LAMBDA(vector)
CL_LISPIFY_NAME("dimension")
CL_DEFMETHOD std::size_t ComplexSingleVector_O::dimension() const { return _Value.dimension(); }

CL_LAMBDA(vector index)
CL_LISPIFY_NAME("vref")
CL_DEFMETHOD c1 ComplexSingleVector_O::vref(std::size_t index) const { return _Value[index]; }

CL_LAMBDA(vector new-value index)
CL_LISPIFY_NAME("setf_vref")
CL_DEFMETHOD c1 ComplexSingleVector_O::setf_vref(c1 new_value, std::size_t index) { return _Value[index] = new_value; }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l1_norm")
CL_DEFMETHOD c1 ComplexSingleVector_O::l1_norm() const { return _Value.l1_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm")
CL_DEFMETHOD c1 ComplexSingleVector_O::l2_norm() const { return _Value.l2_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm_sqr")
CL_DEFMETHOD c1 ComplexSingleVector_O::l2_norm_sqr() const { return _Value.l2_norm_sqr(); }

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("complex-double-vector");
CL_DEFUN ComplexDoubleVector_sp ComplexDoubleVector_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<ComplexDoubleVector_O>::allocate_with_default_constructor();
  res->_Value.resize(args->total_nargs());
  for (size_t i = 0; args->remaining_nargs() > 0; i++) {
    res->_Value[i] = translate::from_object<c2>(args->next_arg())._v;
  }
  return res;
}

CL_LAMBDA(vector)
CL_LISPIFY_NAME("dimension")
CL_DEFMETHOD std::size_t ComplexDoubleVector_O::dimension() const { return _Value.dimension(); }

CL_LAMBDA(vector index)
CL_LISPIFY_NAME("vref")
CL_DEFMETHOD c2 ComplexDoubleVector_O::vref(std::size_t index) const { return _Value[index]; }

CL_LAMBDA(vector new-value index)
CL_LISPIFY_NAME("setf_vref")
CL_DEFMETHOD c2 ComplexDoubleVector_O::setf_vref(c2 new_value, std::size_t index) { return _Value[index] = new_value; }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l1_norm")
CL_DEFMETHOD c2 ComplexDoubleVector_O::l1_norm() const { return _Value.l1_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm")
CL_DEFMETHOD c2 ComplexDoubleVector_O::l2_norm() const { return _Value.l2_norm(); }

CL_LAMBDA(vector)
CL_LISPIFY_NAME("l2_norm_sqr")
CL_DEFMETHOD c2 ComplexDoubleVector_O::l2_norm_sqr() const { return _Value.l2_norm_sqr(); }

// CL_LAMBDA(x y)
CL_LISPIFY_NAME("dot");
CL_DEFUN core::T_sp lila__dot(core::T_sp x, core::T_sp y) {
  if (gctools::IsA<RealSingleVector_sp>(x)) {
    if (gctools::IsA<RealSingleVector_sp>(y)) {
      return core::clasp_make_single_float(
          dot(gctools::As<RealSingleVector_sp>(x)->_Value, gctools::As<RealSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<RealDoubleVector_sp>(y)) {
      return core::clasp_make_double_float(
          dot((r2v)gctools::As<RealSingleVector_sp>(x)->_Value, gctools::As<RealDoubleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexSingleVector_sp>(y)) {
      return translate::to_object<c1>::convert(
          dot((c1v)gctools::As<RealSingleVector_sp>(x)->_Value, gctools::As<ComplexSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexDoubleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot((c2v)gctools::As<RealSingleVector_sp>(x)->_Value, gctools::As<ComplexDoubleVector_sp>(y)->_Value));
    }
  } else if (gctools::IsA<RealDoubleVector_sp>(x)) {
    if (gctools::IsA<RealSingleVector_sp>(y)) {
      return core::clasp_make_double_float(
          dot(gctools::As<RealDoubleVector_sp>(x)->_Value, (r2v)gctools::As<RealSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<RealDoubleVector_sp>(y)) {
      return core::clasp_make_double_float(
          dot(gctools::As<RealDoubleVector_sp>(x)->_Value, gctools::As<RealDoubleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexSingleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot((c2v)gctools::As<RealDoubleVector_sp>(x)->_Value, (c2v)gctools::As<ComplexSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexDoubleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot((c2v)gctools::As<RealDoubleVector_sp>(x)->_Value, gctools::As<ComplexDoubleVector_sp>(y)->_Value));
    }
  } else if (gctools::IsA<ComplexSingleVector_sp>(x)) {
    if (gctools::IsA<RealSingleVector_sp>(y)) {
      return translate::to_object<c1>::convert(
          dot(gctools::As<ComplexSingleVector_sp>(x)->_Value, (c1v)gctools::As<RealSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<RealDoubleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot((c2v)gctools::As<ComplexSingleVector_sp>(x)->_Value, (c2v)gctools::As<RealDoubleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexSingleVector_sp>(y)) {
      return translate::to_object<c1>::convert(
          dot(gctools::As<ComplexSingleVector_sp>(x)->_Value, gctools::As<ComplexSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexDoubleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot((c2v)gctools::As<ComplexSingleVector_sp>(x)->_Value, gctools::As<ComplexDoubleVector_sp>(y)->_Value));
    }
  } else if (gctools::IsA<ComplexDoubleVector_sp>(x)) {
    if (gctools::IsA<RealSingleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot(gctools::As<ComplexDoubleVector_sp>(x)->_Value, (c2v)gctools::As<RealSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<RealDoubleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot(gctools::As<ComplexDoubleVector_sp>(x)->_Value, (c2v)gctools::As<RealDoubleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexSingleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot(gctools::As<ComplexDoubleVector_sp>(x)->_Value, (c2v)gctools::As<ComplexSingleVector_sp>(y)->_Value));
    }
    if (gctools::IsA<ComplexDoubleVector_sp>(y)) {
      return translate::to_object<c2>::convert(
          dot((c2v)gctools::As<ComplexDoubleVector_sp>(x)->_Value, gctools::As<ComplexDoubleVector_sp>(y)->_Value));
    }
  }

  return nil<core::T_O>();
}

} // Namespace lila
