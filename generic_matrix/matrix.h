#ifndef _MINI_MATRIX_H
#define _MINI_MATRIX_H
#include <memory>
/*
  Mini matrix class
*/
template <typename T>
class Matrix {
 public:
  Matrix();
  Matrix(unsigned long i, unsigned long j);
  //~Matrix(void) { if(_data) {std::cerr << "D called\n"; delete [] _data;} }
  T& at(unsigned long i, unsigned long j);
  const T& at(unsigned long i, unsigned long j) const;
 private:
  std::unique_ptr<T[]> _data;
  unsigned long _width;
  unsigned long _height;
};

template <typename T>
Matrix<T>::Matrix()
{
  _width = 0;
  _height = 0;
  _data = nullptr;
}

template <typename T>
Matrix<T>::Matrix(unsigned long i, unsigned long j)
{
  _width = j;
  _height = i;
  _data = std::unique_ptr<T[]>(new T[i * j]);
}

template <typename T>
T& Matrix<T>::at(unsigned long i, unsigned long j)
{
  unsigned long col_offset = j * _height;
  return _data[col_offset + i];
}

template <typename T>
const T& Matrix<T>::at(unsigned long i, unsigned long j) const
{
  unsigned long col_offset = j * _height;
  const T& ret = _data[col_offset + i];
  return ret;
}

#endif
