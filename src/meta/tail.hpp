#ifndef _TAIL_H_
#define _TAIL_H_

#include "type_list.hpp"


template<typename T>
struct tail;

template<typename T, typename ... Ts>
struct tail<type_list<T, Ts...> > {
  typedef type_list<Ts...> type;
};

template<typename T>
struct tail<type_list<T> > {
  typedef empty_type_list type;
};

#endif /* _TAIL_H_ */
