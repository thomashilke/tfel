#ifndef _HOLDER_H_
#define _HOLDER_H_

template<typename T, T v>
struct holder {
  static const T value = v;
};

#endif /* _HOLDER_H_ */
