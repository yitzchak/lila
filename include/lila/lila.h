#ifndef LILA_H
#define LILA_H

#ifdef _TARGET_OS_DARWIN
#include <Accelerate/Accelerate.h>
typedef CBLAS_ORDER CBLAS_LAYOUT;
#else
#include <cblas.h>
#endif

PACKAGE_USE("COMMON-LISP");
NAMESPACE_PACKAGE_ASSOCIATION(lila, lila_pkg, "LILA");
// SYMBOL_SHADOW_EXPORT_SC_(lila_pkg, vector);

namespace lila {

typedef enum { transpose_none = CblasNoTrans, transpose_normal = CblasTrans, transpose_hermitian = CblasConjTrans } transpose_type;

template <class T> class matrix;

template <class T> class vector {
  using V = vector<T>;
  using M = matrix<T>;

  std::size_t _dimension;
  gctools::Vec0<T> _data;

  friend class matrix<T>;

public:
  vector(std::size_t dimension = 0, const T &value = 0) : _dimension(dimension) { _data.resize(_dimension, value); }

  vector(const V &x) : _dimension(x._dimension), _data(x._data) {}

  // template <class U> operator vector<U>() const { return vector<U>(*this); }

  inline std::size_t dimension() const { return _dimension; }
  inline T *data() { return _data.data(); }
  inline const T *data() const { return _data.data(); }

  template <class U> vector(const vector<U> &x) : _dimension(x.dimension()) {
    _data.resize(dimension());
    for (int i = 0; i < dimension(); i++)
      _data[i] = x[i];
  }

  inline T operator[](std::size_t i) const { return _data[i]; }
  inline T &operator[](std::size_t i) { return _data[i]; }

  inline V &resize(std::size_t dimension = 0, const T &value = 0) {
    _data.resize(_dimension = dimension, value);
    return *this;
  }

  // V &update(T a, const V &x);
  template <class U> V &update(U a, const vector<U> &x);
  V &scale(T a);
  V &copy(const V &x);
  V &product(T alpha, const M &A, transpose_type type, const V &x, T beta);

  inline V &operator*=(T a) { return scale(a); }
  inline V &operator+=(const V &x) { return update(1, x); }
  inline V &operator-=(const V &x) { return update(-1, x); }

  T l1_norm() const;
  T l2_norm() const;
  T l2_norm_sqr() const;

  template <class M, class N> inline V &cross(const vector<M> &x, const vector<N> &y) {
    _data.resize(3);
    _dimension = 3;
    _data[0] = static_cast<T>(x[1]) * static_cast<T>(y[2]) - static_cast<T>(x[2]) * static_cast<T>(y[1]);
    _data[1] = static_cast<T>(x[2]) * static_cast<T>(y[0]) - static_cast<T>(x[0]) * static_cast<T>(y[2]);
    _data[2] = static_cast<T>(x[0]) * static_cast<T>(y[1]) - static_cast<T>(x[1]) * static_cast<T>(y[0]);
    return *this;
  }
};

template <class T> class matrix {
  using V = vector<T>;
  using M = matrix<T>;

  std::size_t _row_dimension, _column_dimension, _size;
  bool _column_major = false;
  gctools::Vec0<T> _data;

  friend class vector<T>;

  inline std::size_t index(std::size_t row, std::size_t column) const {
    return _column_major ? column * _row_dimension + row : row * _column_dimension + column;
  }

  CBLAS_LAYOUT layout() const { return _column_major ? CblasColMajor : CblasRowMajor; }

public:
  inline T *data() { return _data.data(); }
  inline const T *data() const { return _data.data(); }
  inline std::size_t column_dimension() const { return _column_dimension; }
  inline std::size_t row_dimension() const { return _row_dimension; }
  inline bool column_major() const { return _column_major; }

  matrix(std::size_t row_dimension = 0, std::size_t column_dimension = 0, const T &value = 0)
      : _row_dimension(row_dimension), _column_dimension(column_dimension) {
    _data.resize(_row_dimension * _column_dimension, value);
  }

  template <class U>
  matrix(const matrix<U> &other)
      : _column_major(other.column_major()), _row_dimension(other.row_dimension()), _column_dimension(other.column_dimension()) {
    _data.resize(_row_dimension * _column_dimension);
    for (int i = 0; i < _row_dimension; i++)
      for (int j = 0; j < _column_dimension; j++)
        _data(i, j) = other(i, j);
  }

  inline M &resize(std::size_t row_dimension = 0, std::size_t column_dimension = 0, const T &value = 0) {
    _row_dimension = row_dimension;
    _column_dimension = column_dimension;
    _data.resize(row_dimension * column_dimension, value);
    return *this;
  }
  inline T operator()(std::size_t row, std::size_t column) const { return data()[index(row, column)]; }
  inline T &operator()(std::size_t row, std::size_t column) { return data()[index(row, column)]; }

  inline void transpose() {
    _column_major = !_column_major;
    std::swap(_row_dimension, _column_dimension);
  }
};

typedef float r1;
typedef double r2;
typedef std::complex<float> c1;
typedef std::complex<double> c2;
typedef vector<r1> r1v;
typedef vector<r2> r2v;
typedef vector<c1> c1v;
typedef vector<c2> c2v;
typedef matrix<r1> r1m;
typedef matrix<r2> r2m;
typedef matrix<c1> c1m;
typedef matrix<c2> c2m;

template <> template <> inline r1v &r1v::update(r1 a, const r1v &x) {
  cblas_saxpy(std::min(x.dimension(), dimension()), a, x.data(), 1, data(), 1);
  return *this;
}

template <> template <> inline r2v &r2v::update(r2 a, const r2v &x) {
  cblas_daxpy(std::min(x.dimension(), dimension()), a, x.data(), 1, data(), 1);
  return *this;
}

template <> template <> inline c1v &c1v::update(c1 a, const c1v &x) {
  cblas_caxpy(std::min(x.dimension(), dimension()), &a, x.data(), 1, data(), 1);
  return *this;
}

template <> template <> inline c1v &c1v::update(r1 a, const r1v &x) {
  cblas_saxpy(std::min(x.dimension(), dimension()), a, x.data(), 1, reinterpret_cast<r1 *>(data()), 2);
  return *this;
}

template <> template <> inline c2v &c2v::update(c2 a, const c2v &x) {
  cblas_zaxpy(std::min(x.dimension(), dimension()), &a, x.data(), 1, data(), 1);
  return *this;
}

template <> template <> inline c2v &c2v::update(r2 a, const r2v &x) {
  cblas_daxpy(std::min(x.dimension(), dimension()), a, x.data(), 1, reinterpret_cast<r2 *>(data()), 2);
  return *this;
}

template <> inline r1v &r1v::scale(r1 x) {
  cblas_sscal(dimension(), x, data(), 1);
  return *this;
}

template <> inline r2v &r2v::scale(r2 x) {
  cblas_dscal(dimension(), x, data(), 1);
  return *this;
}

template <> inline c1v &c1v::scale(c1 x) {
  cblas_cscal(dimension(), &x, data(), 1);
  return *this;
}

template <> inline c2v &c2v::scale(c2 x) {
  cblas_zscal(dimension(), &x, data(), 1);
  return *this;
}

template <> inline r1v &r1v::copy(const r1v &x) {
  cblas_scopy(std::min(x.dimension(), dimension()), x.data(), 1, data(), 1);
  if (dimension() > x.dimension())
    std::fill(data() + x.dimension(), data() + dimension(), 0.0f);
  return *this;
}

template <> inline r2v &r2v::copy(const r2v &x) {
  cblas_dcopy(std::min(x.dimension(), dimension()), x.data(), 1, data(), 1);
  if (dimension() > x.dimension())
    std::fill(data() + x.dimension(), data() + dimension(), 0.0);
  return *this;
}

template <> inline c1v &c1v::copy(const c1v &x) {
  cblas_ccopy(std::min(x.dimension(), dimension()), x.data(), 1, data(), 1);
  if (dimension() > x.dimension())
    std::fill(data() + x.dimension(), data() + dimension(), 0.0f);
  return *this;
}

template <> inline c2v &c2v::copy(const c2v &x) {
  cblas_zcopy(std::min(x.dimension(), dimension()), x.data(), 1, data(), 1);
  if (dimension() > x.dimension())
    std::fill(data() + x.dimension(), data() + dimension(), 0.0);
  return *this;
}

template <> inline r1v &r1v::product(r1 alpha, const r1m &A, transpose_type type, const r1v &x, r1 beta) {
  cblas_sgemv(A.layout(), static_cast<CBLAS_TRANSPOSE>(type), A.row_dimension(), A.column_dimension(), alpha, A.data(),
              A.row_dimension(), x.data(), 1, beta, data(), 1);
  return *this;
}

template <> inline r2v &r2v::product(r2 alpha, const r2m &A, transpose_type type, const r2v &x, r2 beta) {
  cblas_dgemv(A.layout(), static_cast<CBLAS_TRANSPOSE>(type), A.row_dimension(), A.column_dimension(), alpha, A.data(),
              A.row_dimension(), x.data(), 1, beta, data(), 1);
  return *this;
}

template <> inline c1v &c1v::product(c1 alpha, const c1m &A, transpose_type type, const c1v &x, c1 beta) {
  cblas_cgemv(A.layout(), static_cast<CBLAS_TRANSPOSE>(type), A.row_dimension(), A.column_dimension(), &alpha, A.data(),
              A.row_dimension(), x.data(), 1, &beta, data(), 1);
  return *this;
}

template <> inline c2v &c2v::product(c2 alpha, const c2m &A, transpose_type type, const c2v &x, c2 beta) {
  cblas_zgemv(A.layout(), static_cast<CBLAS_TRANSPOSE>(type), A.row_dimension(), A.column_dimension(), &alpha, A.data(),
              A.row_dimension(), x.data(), 1, &beta, data(), 1);
  return *this;
}

template <> inline r1 r1v::l1_norm() const { return cblas_sasum(dimension(), data(), 1); }

template <> inline r2 r2v::l1_norm() const { return cblas_dasum(dimension(), data(), 1); }

template <> inline c1 c1v::l1_norm() const { return cblas_scasum(dimension(), data(), 1); }

template <> inline c2 c2v::l1_norm() const { return cblas_dzasum(dimension(), data(), 1); }

template <> inline r1 r1v::l2_norm() const { return cblas_snrm2(dimension(), data(), 1); }

template <> inline r2 r2v::l2_norm() const { return cblas_dnrm2(dimension(), data(), 1); }

template <> inline c1 c1v::l2_norm() const { return cblas_scnrm2(dimension(), data(), 1); }

template <> inline c2 c2v::l2_norm() const { return cblas_dznrm2(dimension(), data(), 1); }

template <> inline r1 r1v::l2_norm_sqr() const { return cblas_sdot(dimension(), data(), 1, data(), 1); }

template <> inline r2 r2v::l2_norm_sqr() const { return cblas_ddot(dimension(), data(), 1, data(), 1); }

template <> inline c1 c1v::l2_norm_sqr() const {
  c1 result;
  cblas_cdotc_sub(dimension(), data(), 1, data(), 1, &result);
  return result;
}

template <> inline c2 c2v::l2_norm_sqr() const {
  c2 result;
  cblas_zdotc_sub(dimension(), data(), 1, data(), 1, &result);
  return result;
}

template <class T> inline vector<T> &cross(const vector<T> &x, const vector<T> &y) {
  vector<T> result(3);
  result.cross(x, y);
  return result;
}

template <class T> T triple(const vector<T> &x, const vector<T> &y, const vector<T> &z) {
  return x[0] * y[1] * z[2] - x[1] * y[0] * z[2] + x[2] * y[0] * z[1] - x[0] * y[2] * z[1] + +x[1] * y[2] * z[0] -
         x[2] * y[1] * z[0];
}

inline r1 dot(const r1v &x, const r1v &y) { return cblas_sdot(std::min(x.dimension(), y.dimension()), x.data(), 1, y.data(), 1); }

inline r2 dot(const r2v &x, const r2v &y) { return cblas_ddot(std::min(x.dimension(), y.dimension()), x.data(), 1, y.data(), 1); }

inline c1 dot(const c1v &x, const c1v &y) {
  c1 result;
  cblas_cdotu_sub(std::min(y.dimension(), x.dimension()), y.data(), 1, x.data(), 1, &result);
  return result;
}

inline c2 dot(const c2v &x, const c2v &y) {
  c2 result;
  cblas_zdotu_sub(std::min(y.dimension(), x.dimension()), y.data(), 1, x.data(), 1, &result);
  return result;
}

inline r1 dotc(const r1v &x, const r1v &y) { return cblas_sdot(std::min(x.dimension(), y.dimension()), x.data(), 1, y.data(), 1); }

inline r2 dotc(const r2v &x, const r2v &y) { return cblas_ddot(std::min(x.dimension(), y.dimension()), x.data(), 1, y.data(), 1); }

inline c1 dotc(const c1v &x, const c1v &y) {
  c1 result;
  cblas_cdotc_sub(std::min(y.dimension(), x.dimension()), y.data(), 1, x.data(), 1, &result);
  return result;
}

inline c2 dotc(const c2v &x, const c2v &y) {
  c2 result;
  cblas_zdotc_sub(std::min(y.dimension(), x.dimension()), y.data(), 1, x.data(), 1, &result);
  return result;
}

SMART(Vector)

class Vector_O : public core::CxxObject_O {
  LISP_ABSTRACT_CLASS(lila, lila_pkg, Vector_O, "lvector", core::CxxObject_O);

public:
  Vector_O() {}
};

SMART(RealSingleVector);

class RealSingleVector_O : public Vector_O {
  LISP_CLASS(lila, lila_pkg, RealSingleVector_O, "REAL-SINGLE-VECTOR", Vector_O);
  r1v _Value;

public:
  RealSingleVector_O() {}

  static RealSingleVector_sp make(core::Vaslist_sp args);

  r1v &value() { return _Value; }
  const r1v &value() const { return _Value; }

  std::size_t dimension() const;

  float vref(std::size_t index) const;
  float setf_vref(float new_value, std::size_t index);

  r1 l1_norm() const;
  r1 l2_norm() const;
  r1 l2_norm_sqr() const;
};

SMART(RealDoubleVector);

class RealDoubleVector_O : public Vector_O {
  LISP_CLASS(lila, lila_pkg, RealDoubleVector_O, "REAL-DOUBLE-VECTOR", Vector_O);
  r2v _Value;

public:
  RealDoubleVector_O() {}

  static RealDoubleVector_sp make(core::Vaslist_sp args);

  std::size_t dimension() const;

  double vref(std::size_t index) const;
  double setf_vref(double new_value, std::size_t index);

  r2 l1_norm() const;
  r2 l2_norm() const;
  r2 l2_norm_sqr() const;
};

SMART(ComplexSingleVector);

class ComplexSingleVector_O : public Vector_O {
  LISP_CLASS(lila, lila_pkg, ComplexSingleVector_O, "COMPLEX-SINGLE-VECTOR", Vector_O);
  c1v _Value;

public:
  ComplexSingleVector_O() {}

  static ComplexSingleVector_sp make(core::Vaslist_sp args);

  std::size_t dimension() const;

  c1 vref(std::size_t index) const;
  c1 setf_vref(c1 new_value, std::size_t index);

  c1 l1_norm() const;
  c1 l2_norm() const;
  c1 l2_norm_sqr() const;
};

SMART(ComplexDoubleVector);

class ComplexDoubleVector_O : public Vector_O {
  LISP_CLASS(lila, lila_pkg, ComplexDoubleVector_O, "COMPLEX-DOUBLE-VECTOR", Vector_O);
  c2v _Value;

public:
  ComplexDoubleVector_O() {}

  static ComplexDoubleVector_sp make(core::Vaslist_sp args);

  std::size_t dimension() const;

  c2 vref(std::size_t index) const;
  c2 setf_vref(c2 new_value, std::size_t index);

  c2 l1_norm() const;
  c2 l2_norm() const;
  c2 l2_norm_sqr() const;
};

SMART(Matrix)

class Matrix_O : public core::CxxObject_O {
  LISP_ABSTRACT_CLASS(lila, lila_pkg, Matrix_O, "matrix", core::CxxObject_O);

public:
  Matrix_O() {}
};

SMART(RealSingleMatrix);

class RealSingleMatrix_O : public Matrix_O {
  LISP_CLASS(lila, lila_pkg, RealSingleMatrix_O, "REAL-SINGLE-MATRIX", Matrix_O);
  r1m _Value;

public:
  RealSingleMatrix_O() {}

  static RealSingleMatrix_sp make(core::Vaslist_sp args);

  r1m &value() { return _Value; }
  const r1m &value() const { return _Value; }

  std::size_t row_dimension() const;
  std::size_t column_dimension() const;

  r1 mref(std::size_t row, std::size_t column) const;
  r1 setf_mref(r1 new_value, std::size_t row, std::size_t column);
};

SMART(RealDoubleMatrix);

class RealDoubleMatrix_O : public Matrix_O {
  LISP_CLASS(lila, lila_pkg, RealDoubleMatrix_O, "REAL-DOUBLE-MATRIX", Matrix_O);
  r2m _Value;

public:
  RealDoubleMatrix_O() {}

  static RealDoubleMatrix_sp make(core::Vaslist_sp args);

  std::size_t row_dimension() const;
  std::size_t column_dimension() const;

  r2 mref(std::size_t row, std::size_t column) const;
  r2 setf_mref(r2 new_value, std::size_t row, std::size_t column);
};

SMART(ComplexSingleMatrix);

class ComplexSingleMatrix_O : public Matrix_O {
  LISP_CLASS(lila, lila_pkg, ComplexSingleMatrix_O, "COMPLEX-SINGLE-MATRIX", Matrix_O);
  c1m _Value;

public:
  ComplexSingleMatrix_O() {}

  static ComplexSingleMatrix_sp make(core::Vaslist_sp args);

  std::size_t row_dimension() const;
  std::size_t column_dimension() const;

  c1 mref(std::size_t row, std::size_t column) const;
  c1 setf_mref(c1 new_value, std::size_t row, std::size_t column);
};

SMART(ComplexDoubleMatrix);

class ComplexDoubleMatrix_O : public Matrix_O {
  LISP_CLASS(lila, lila_pkg, ComplexDoubleMatrix_O, "COMPLEX-DOUBLE-MATRIX", Matrix_O);
  c2m _Value;

public:
  ComplexDoubleMatrix_O() {}

  static ComplexDoubleMatrix_sp make(core::Vaslist_sp args);

  std::size_t row_dimension() const;
  std::size_t column_dimension() const;

  c2 mref(std::size_t row, std::size_t column) const;
  c2 setf_mref(c2 new_value, std::size_t row, std::size_t column);
};

} // namespace lila

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

#endif
