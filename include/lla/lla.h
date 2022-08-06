#ifndef LLA_H
#define LLA_H

#include <cblas.h>

namespace lla {

template <class T> class vector {
  using V = vector<T>;

  static std::allocator<T> alloc;

  std::size_t _dimension;
  bool _row;
  T *data;

public:
  vector(std::size_t dimension, bool row = false) : _dimension(dimension), _row(row) { data = alloc.allocate(_dimension); }
  vector(const V &other);

  ~vector() { alloc.deallocate(data, _dimension); }

  std::size_t dimension() { return _dimension; }
  bool row() { return _row; }

  T operator[](std::size_t i) const { return data[i]; }
  T &operator[](std::size_t i) { return data[i]; }

  V &operator*=(T x);
  V &operator+=(const V &x);
  V &operator-=(const V &x);

  V &transpose() {
    _row = !_row;
    return *this;
  }
};

template <class T> std::allocator<T> vector<T>::alloc = std::allocator<T>();

template <class T, class Allocator = std::allocator<T>> class matrix {
  std::size_t row_dimension, column_dimension;
  bool row;
  T *data;

public:
};

template <> vector<float>::vector(const vector<float> &other) {
  _row = other._row;
  data = alloc.allocate(_dimension = other._dimension);
  cblas_scopy(_dimension, other.data, 1, data, 1);
}

template <> vector<double>::vector(const vector<double> &other) {
  _row = other._row;
  data = alloc.allocate(_dimension = other._dimension);
  cblas_dcopy(_dimension, other.data, 1, data, 1);
}

template <> vector<float> &vector<float>::operator*=(float x) {
  cblas_sscal(_dimension, x, data, 1);
  return *this;
}

template <> vector<double> &vector<double>::operator*=(double x) {
  cblas_dscal(_dimension, x, data, 1);
  return *this;
}

template <> vector<float> &vector<float>::operator+=(const vector<float> &x) {
  cblas_saxpy(_dimension, 1.0, x.data, 1, data, 1);
  return *this;
}

template <> vector<double> &vector<double>::operator+=(const vector<double> &x) {
  cblas_daxpy(_dimension, 1.0, x.data, 1, data, 1);
  return *this;
}

template <> vector<float> &vector<float>::operator-=(const vector<float> &x) {
  cblas_saxpy(_dimension, -1.0, x.data, 1, data, 1);
  return *this;
}

template <> vector<double> &vector<double>::operator-=(const vector<double> &x) {
  cblas_daxpy(_dimension, -1.0, x.data, 1, data, 1);
  return *this;
}

} // namespace lla

#endif
