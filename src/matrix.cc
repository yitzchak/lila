#include <clasp/clasp.h>
#include <lila/lila.h>

namespace lila {

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("real-single-matrix");
CL_DEFUN RealSingleMatrix_sp RealSingleMatrix_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<RealSingleMatrix_O>::allocate_with_default_constructor();
  auto row_dimension = args->total_nargs();
  if (row_dimension > 0) {
    auto column_dimension = gctools::As<core::Cons_sp>(args->next_arg_indexed(0))->length();
    res->_Value.resize(row_dimension, column_dimension);
    for (size_t row = 0; args->remaining_nargs() > 0; row++) {
      for (auto [column, cons] = std::tuple{0, args->next_arg()}; column < column_dimension; column++, cons = core::oCdr(cons)) {
        res->_Value(row, column) = core::clasp_to_float(oCar(cons));
      }
    }
  }
  return res;
}

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("row-dimension")
CL_DEFMETHOD std::size_t RealSingleMatrix_O::row_dimension() const { return _Value.row_dimension(); }

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("column-dimension")
CL_DEFMETHOD std::size_t RealSingleMatrix_O::column_dimension() const { return _Value.column_dimension(); }

CL_LAMBDA(matrix row column)
CL_LISPIFY_NAME("mref")
CL_DEFMETHOD r1 RealSingleMatrix_O::mref(std::size_t row, std::size_t column) const { return _Value(row, column); }

CL_LAMBDA(matrix new - value row column)
CL_LISPIFY_NAME("setf_mref")
CL_DEFMETHOD r1 RealSingleMatrix_O::setf_mref(r1 new_value, std::size_t row, std::size_t column) {
  return _Value(row, column) = new_value;
}

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("real-double-matrix");
CL_DEFUN RealDoubleMatrix_sp RealDoubleMatrix_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<RealDoubleMatrix_O>::allocate_with_default_constructor();
  auto row_dimension = args->total_nargs();
  if (row_dimension > 0) {
    auto column_dimension = gctools::As<core::Cons_sp>(args->next_arg_indexed(0))->length();
    res->_Value.resize(row_dimension, column_dimension);
    for (size_t row = 0; args->remaining_nargs() > 0; row++) {
      for (auto [column, cons] = std::tuple{0, args->next_arg()}; column < column_dimension; column++, cons = core::oCdr(cons)) {
        res->_Value(row, column) = core::clasp_to_double(oCar(cons));
      }
    }
  }
  return res;
}

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("row-dimension")
CL_DEFMETHOD std::size_t RealDoubleMatrix_O::row_dimension() const { return _Value.row_dimension(); }

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("column-dimension")
CL_DEFMETHOD std::size_t RealDoubleMatrix_O::column_dimension() const { return _Value.column_dimension(); }

CL_LAMBDA(matrix row column)
CL_LISPIFY_NAME("mref")
CL_DEFMETHOD r2 RealDoubleMatrix_O::mref(std::size_t row, std::size_t column) const { return _Value(row, column); }

CL_LAMBDA(matrix new - value row column)
CL_LISPIFY_NAME("setf_mref")
CL_DEFMETHOD r2 RealDoubleMatrix_O::setf_mref(r2 new_value, std::size_t row, std::size_t column) {
  return _Value(row, column) = new_value;
}

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("complex-single-matrix");
CL_DEFUN ComplexSingleMatrix_sp ComplexSingleMatrix_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<ComplexSingleMatrix_O>::allocate_with_default_constructor();
  auto row_dimension = args->total_nargs();
  if (row_dimension > 0) {
    auto column_dimension = gctools::As<core::Cons_sp>(args->next_arg_indexed(0))->length();
    res->_Value.resize(row_dimension, column_dimension);
    for (size_t row = 0; args->remaining_nargs() > 0; row++) {
      for (auto [column, cons] = std::tuple{0, args->next_arg()}; column < column_dimension; column++, cons = core::oCdr(cons)) {
        res->_Value(row, column) = translate::from_object<c1>(oCar(cons))._v;
      }
    }
  }
  return res;
}

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("row-dimension")
CL_DEFMETHOD std::size_t ComplexSingleMatrix_O::row_dimension() const { return _Value.row_dimension(); }

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("column-dimension")
CL_DEFMETHOD std::size_t ComplexSingleMatrix_O::column_dimension() const { return _Value.column_dimension(); }

CL_LAMBDA(matrix row column)
CL_LISPIFY_NAME("mref")
CL_DEFMETHOD c1 ComplexSingleMatrix_O::mref(std::size_t row, std::size_t column) const { return _Value(row, column); }

CL_LAMBDA(matrix new - value row column)
CL_LISPIFY_NAME("setf_mref")
CL_DEFMETHOD c1 ComplexSingleMatrix_O::setf_mref(c1 new_value, std::size_t row, std::size_t column) {
  return _Value(row, column) = new_value;
}

CL_LAMBDA(core:&va-rest args)
CL_LISPIFY_NAME("complex-double-matrix");
CL_DEFUN ComplexDoubleMatrix_sp ComplexDoubleMatrix_O::make(core::Vaslist_sp args) {
  auto res = gctools::GC<ComplexDoubleMatrix_O>::allocate_with_default_constructor();
  auto row_dimension = args->total_nargs();
  if (row_dimension > 0) {
    auto column_dimension = gctools::As<core::Cons_sp>(args->next_arg_indexed(0))->length();
    res->_Value.resize(row_dimension, column_dimension);
    for (size_t row = 0; args->remaining_nargs() > 0; row++) {
      for (auto [column, cons] = std::tuple{0, args->next_arg()}; column < column_dimension; column++, cons = core::oCdr(cons)) {
        res->_Value(row, column) = translate::from_object<c2>(oCar(cons))._v;
      }
    }
  }
  return res;
}

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("row-dimension")
CL_DEFMETHOD std::size_t ComplexDoubleMatrix_O::row_dimension() const { return _Value.row_dimension(); }

CL_LAMBDA(matrix)
CL_LISPIFY_NAME("column-dimension")
CL_DEFMETHOD std::size_t ComplexDoubleMatrix_O::column_dimension() const { return _Value.column_dimension(); }

CL_LAMBDA(matrix row column)
CL_LISPIFY_NAME("mref")
CL_DEFMETHOD c2 ComplexDoubleMatrix_O::mref(std::size_t row, std::size_t column) const { return _Value(row, column); }

CL_LAMBDA(matrix new - value row column)
CL_LISPIFY_NAME("setf_mref")
CL_DEFMETHOD c2 ComplexDoubleMatrix_O::setf_mref(c2 new_value, std::size_t row, std::size_t column) {
  return _Value(row, column) = new_value;
}

}
