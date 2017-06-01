#ifndef _HEAD_H_
#define _HEAD_H_

#include "type_list.hpp"


template<typename L>
struct head;

template<typename T, typename ... Ts>
struct head<type_list<T, Ts...> > { typedef T type; };

template<typename T>
struct head<type_list<T> > { typedef T type; };

#endif /* _HEAD_H_ */
