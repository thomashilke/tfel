#ifndef _AT_H_
#define _AT_H_

#include "type_list.hpp"
#include "head.hpp"
#include "tail.hpp"


template<unsigned int i, typename T>
struct at {
  typedef typename at<i - 1, typename tail<T>::type>::type type;
};

template<typename T>
struct at<0, T> {
  typedef typename head<T>::type type;
};

#endif /* _AT_H_ */
