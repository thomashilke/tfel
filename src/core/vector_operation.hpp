#ifndef _VECTOR_OPERATION_H_
#define _VECTOR_OPERATION_H_

#include <cassert>

#include <spikes/array.hpp>


template<typename T>
T dotp(const array<T>& a, const array<T>& b) {
  assert(a.get_rank() == 1);
  assert(b.get_rank() == 1);
  assert(a.get_size(0) == b.get_size(0));
  
  T result = {};

  for (std::size_t k(0); k < a.get_size(0); ++k)
    result += a.at(k) * b.at(k);
  
  return result;
}


template<typename T>
array<T> crossp(const array<T>& a, const array<T>& b) {
  assert(a.get_rank() == 1);
  assert(b.get_rank() == 1);
  assert(a.get_size(0) == 3);
  assert(a.get_size(0) == 3);
  
  array<T> result{3};

  result.at(0) = a.at(1) * b.at(2) - a.at(2) * b.at(1);
  result.at(1) = a.at(2) * b.at(0) - a.at(0) * b.at(2);
  result.at(2) = a.at(0) * b.at(1) - a.at(1) * b.at(0);

  return result;
}

#endif /* _VECTOR_OPERATION_H_ */
