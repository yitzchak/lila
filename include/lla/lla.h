#ifndef LLA_H
#define LLA_H

#ifdef _TARGET_OS_DARWIN
#include <Accelerate/Accelerate.h>
typedef CBLAS_ORDER CBLAS_LAYOUT;
#else
#include <cblas.h>
#endif

namespace lla {

template <class T> class matrix;

template <class T> class vector {
  using V = vector<T>;
  using M = matrix<T>;

  std::size_t _dimension;
  std::vector<T> data;

  friend class matrix<T>;

public:
  vector(std::size_t dimension, const T &value = 0) : _dimension(dimension), data(dimension, value) {}
  vector(const V &other) : _dimension(other._dimension), data(other.data) {}

  std::size_t dimension() const { return _dimension; }

  T operator[](std::size_t i) const { return data[i]; }
  T &operator[](std::size_t i) { return data[i]; }

  V &operator*=(T x);
  V &operator+=(const V &x);
  V &operator-=(const V &x);

  T inner(const V &other) const;
  T dot(const V &other) const { return dot(other); }
  M outer(const V &other) const;
};

template <class T> class matrix {
  using V = vector<T>;
  using M = matrix<T>;

  std::size_t _row_dimension, _column_dimension;
  bool _column_major = false;
  std::vector<T> data;

  friend class vector<T>;

  std::size_t index(std::size_t row, std::size_t column) const {
    return _column_major ? column * _row_dimension + row : row * _column_dimension + column;
  }

  CBLAS_LAYOUT layout() const { return _column_major ? CblasColMajor : CblasRowMajor; }

public:
  matrix(std::size_t row_dimension, std::size_t column_dimension, const T &value = 0)
      : data(row_dimension * column_dimension, value), _row_dimension(row_dimension), _column_dimension(column_dimension) {}
  matrix(const M &other)
      : _column_major(other._column_major), _row_dimension(other._row_dimension), _column_dimension(other._column_dimension),
        data(other.data) {}

  std::size_t column_dimension() const { return _column_dimension; }
  std::size_t row_dimension() const { return _row_dimension; }
  
  T operator()(std::size_t row, std::size_t column) const { return data[index(row, column)]; }
  T &operator()(std::size_t row, std::size_t column) { return data[index(row, column)]; }

  void transpose() {
    _column_major = !_column_major;
    std::swap(_row_dimension, _column_dimension);
  }
};

template <> vector<float> &vector<float>::operator*=(float x) {
  cblas_sscal(dimension(), x, data.data(), 1);
  return *this;
}

template <> vector<double> &vector<double>::operator*=(double x) {
  cblas_dscal(dimension(), x, data.data(), 1);
  return *this;
}

template <> vector<float> &vector<float>::operator+=(const vector<float> &x) {
  cblas_saxpy(dimension(), 1.0f, x.data.data(), 1, data.data(), 1);
  return *this;
}

template <> vector<double> &vector<double>::operator+=(const vector<double> &x) {
  cblas_daxpy(dimension(), 1.0, x.data.data(), 1, data.data(), 1);
  return *this;
}

template <> vector<float> &vector<float>::operator-=(const vector<float> &x) {
  cblas_saxpy(dimension(), -1.0f, x.data.data(), 1, data.data(), 1);
  return *this;
}

template <> vector<double> &vector<double>::operator-=(const vector<double> &x) {
  cblas_daxpy(dimension(), -1.0, x.data.data(), 1, data.data(), 1);
  return *this;
}

template <> float vector<float>::inner(const vector<float> &other) const {
  return cblas_sdot(dimension(), other.data.data(), 1, data.data(), 1);
}

template <> double vector<double>::inner(const vector<double> &other) const {
  return cblas_ddot(dimension(), other.data.data(), 1, data.data(), 1);
}

template <> matrix<float> vector<float>::outer(const vector<float> &other) const {
  matrix<float> result(dimension(), dimension());

  cblas_sger(result.layout(), dimension(), dimension(), 1.0f, data.data(), 1, other.data.data(), 1, result.data.data(), dimension());

  return result;
}

} // namespace lla

#endif
