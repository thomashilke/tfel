#ifndef _SIZE_H_
#define _SIZE_H_

#include "type_list.hpp"
#include "holder.hpp"


template<unsigned int i, typename T>
struct size_impl: size_impl<i + 1, typename tail<T>::type> {};

template<unsigned int i>
struct size_impl<i, empty_type_list> { typedef holder<unsigned int, i> type; };

template<typename T>
struct size { typedef typename size_impl<0, T>::type type; };

#endif /* _SIZE_H_ */
