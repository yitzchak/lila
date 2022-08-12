#ifndef LILA_H
#define LILA_H

#ifdef _TARGET_OS_DARWIN
#include <Accelerate/Accelerate.h>
typedef CBLAS_ORDER CBLAS_LAYOUT;
#else
#include <cblas.h>
#endif

namespace lila {

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

  template <class U> vector(const vector<U> &x) : _dimension(x.dimension()) {
    _data.resize(_dimension);
    for (int i = 0; i < _dimension; i++)
      _data[i] = x[i];
  }

  std::size_t dimension() const { return _dimension; }

  T operator[](std::size_t i) const { return _data[i]; }
  T &operator[](std::size_t i) { return _data[i]; }

  V &resize(std::size_t dimension = 0, const T &value = 0) {
    _data.resize(_dimension = dimension, value);
    return *this;
  }

  V &update(T a, const V &x);
  V &scale(T a);
  V &copy(const V &x);

  V &operator*=(T a) { return scale(a); }
  V &operator+=(const V &x) { return update(1, x); }
  V &operator-=(const V &x) { return update(-1, x); }

  T l1_norm() const;
  T l2_norm() const;
  T l2_norm_sqr() const;

  T dot(const V &x) const;
};

template <class T> class matrix {
  using V = vector<T>;
  using M = matrix<T>;

  static std::allocator<T> alloc;

  std::size_t _row_dimension, _column_dimension, _size;
  bool _column_major = false;
  T *_data;

  friend class vector<T>;

  std::size_t index(std::size_t row, std::size_t column) const {
    return _column_major ? column * _row_dimension + row : row * _column_dimension + column;
  }

  CBLAS_LAYOUT layout() const { return _column_major ? CblasColMajor : CblasRowMajor; }

public:
  matrix(std::size_t row_dimension, std::size_t column_dimension, const T &value = 0)
      : _row_dimension(row_dimension), _column_dimension(column_dimension) {
    _size = _row_dimension * _column_dimension;
    _data = alloc.allocate(_size);
    std::fill(_data, _data + _size, value);
  }

  matrix(const M &other)
      : _column_major(other._column_major), _row_dimension(other._row_dimension), _column_dimension(other._column_dimension),
        _size(other._size) {}

  ~matrix() { alloc.deallocate(_data, _size); }

  std::size_t column_dimension() const { return _column_dimension; }
  std::size_t row_dimension() const { return _row_dimension; }

  T operator()(std::size_t row, std::size_t column) const { return _data[index(row, column)]; }
  T &operator()(std::size_t row, std::size_t column) { return _data[index(row, column)]; }

  void transpose() {
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

template <class T> std::allocator<T> matrix<T>::alloc = std::allocator<T>();

template <> r1v &r1v::update(r1 a, const r1v &x) {
  cblas_saxpy(std::min(x._dimension, _dimension), a, x._data.data(), 1, _data.data(), 1);
  return *this;
}

template <> r2v &r2v::update(r2 a, const r2v &x) {
  cblas_daxpy(std::min(x._dimension, _dimension), a, x._data.data(), 1, _data.data(), 1);
  return *this;
}

template <> c1v &c1v::update(c1 a, const c1v &x) {
  cblas_caxpy(std::min(x._dimension, _dimension), &a, x._data.data(), 1, _data.data(), 1);
  return *this;
}

template <> c2v &c2v::update(c2 a, const c2v &x) {
  cblas_zaxpy(std::min(x._dimension, _dimension), &a, x._data.data(), 1, _data.data(), 1);
  return *this;
}

template <> r1v &r1v::scale(r1 x) {
  cblas_sscal(_dimension, x, _data.data(), 1);
  return *this;
}

template <> r2v &r2v::scale(r2 x) {
  cblas_dscal(_dimension, x, _data.data(), 1);
  return *this;
}

template <> c1v &c1v::scale(c1 x) {
  cblas_cscal(_dimension, &x, _data.data(), 1);
  return *this;
}

template <> c2v &c2v::scale(c2 x) {
  cblas_zscal(_dimension, &x, _data.data(), 1);
  return *this;
}

template <> r1v &r1v::copy(const r1v &x) {
  cblas_scopy(std::min(x._dimension, _dimension), x._data.data(), 1, _data.data(), 1);
  if (_dimension > x._dimension)
    std::fill(_data.data() + x._dimension, _data.data() + _dimension, 0.0f);
  return *this;
}

template <> r2v &r2v::copy(const r2v &x) {
  cblas_dcopy(std::min(x._dimension, _dimension), x._data.data(), 1, _data.data(), 1);
  if (_dimension > x._dimension)
    std::fill(_data.data() + x._dimension, _data.data() + _dimension, 0.0);
  return *this;
}

template <> c1v &c1v::copy(const c1v &x) {
  cblas_ccopy(std::min(x._dimension, _dimension), x._data.data(), 1, _data.data(), 1);
  if (_dimension > x._dimension)
    std::fill(_data.data() + x._dimension, _data.data() + _dimension, 0.0f);
  return *this;
}

template <> c2v &c2v::copy(const c2v &x) {
  cblas_zcopy(std::min(x._dimension, _dimension), x._data.data(), 1, _data.data(), 1);
  if (_dimension > x._dimension)
    std::fill(_data.data() + x._dimension, _data.data() + _dimension, 0.0);
  return *this;
}

template <> r1 r1v::l1_norm() const { return cblas_sasum(_dimension, _data.data(), 1); }

template <> r2 r2v::l1_norm() const { return cblas_dasum(_dimension, _data.data(), 1); }

template <> c1 c1v::l1_norm() const { return cblas_scasum(_dimension, _data.data(), 1); }

template <> c2 c2v::l1_norm() const { return cblas_dzasum(_dimension, _data.data(), 1); }

template <> r1 r1v::l2_norm() const { return cblas_snrm2(_dimension, _data.data(), 1); }

template <> r2 r2v::l2_norm() const { return cblas_dnrm2(_dimension, _data.data(), 1); }

template <> c1 c1v::l2_norm() const { return cblas_scnrm2(_dimension, _data.data(), 1); }

template <> c2 c2v::l2_norm() const { return cblas_dznrm2(_dimension, _data.data(), 1); }

template <> r1 r1v::l2_norm_sqr() const { return cblas_sdot(_dimension, _data.data(), 1, _data.data(), 1); }

template <> r2 r2v::l2_norm_sqr() const { return cblas_ddot(_dimension, _data.data(), 1, _data.data(), 1); }

template <> c1 c1v::l2_norm_sqr() const {
  c1 result;
  cblas_cdotc_sub(_dimension, _data.data(), 1, _data.data(), 1, &result);
  return result;
}

template <> c2 c2v::l2_norm_sqr() const {
  c2 result;
  cblas_zdotc_sub(_dimension, _data.data(), 1, _data.data(), 1, &result);
  return result;
}

template <> r1 r1v::dot(const r1v &x) const {
  return cblas_sdot(std::min(x._dimension, _dimension), x._data.data(), 1, _data.data(), 1);
}

template <> r2 r2v::dot(const r2v &x) const {
  return cblas_ddot(std::min(x._dimension, _dimension), x._data.data(), 1, _data.data(), 1);
}

template <> c1 c1v::dot(const c1v &x) const {
  c1 result;
  cblas_cdotc_sub(std::min(_dimension, x._dimension), _data.data(), 1, x._data.data(), 1, &result);
  return result;
}

template <> c2 c2v::dot(const c2v &x) const {
  c2 result;
  cblas_zdotc_sub(std::min(_dimension, x._dimension), _data.data(), 1, x._data.data(), 1, &result);
  return result;
}

} // namespace lila

PACKAGE_USE("COMMON-LISP");
NAMESPACE_PACKAGE_ASSOCIATION(lila, lila_pkg, "LILA");

#include <lila/vector.fwd.h>
#include <lila/vector.h>

#endif
