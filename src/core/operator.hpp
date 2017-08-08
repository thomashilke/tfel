#ifndef _OPERATOR_H_
#define _OPERATOR_H_

template<typename value_t>
struct add {
  using value_type = value_t;
  static value_type apply(value_type x, value_type y) { return x + y; }
};

template<typename value_t>
struct substract {
  using value_type = value_t;
  static value_type apply(value_type x, value_type y) { return x - y; }
};

template<typename value_t>
struct multiply {
  using value_type = value_t;
  static value_type apply(value_type x, value_type y) { return x * y; }
};

template<typename value_t>
struct divide {
  using value_type = value_t;
  static value_type apply(value_type x, value_type y) { return x / y; }
};

template<typename argument_type>
struct greater_than {
  using value_type = bool;
  static bool apply(argument_type x, argument_type y) { return x > y; }
};

template<typename argument_type>
struct greater_equal_than {
  using value_type = bool;
  static bool apply(argument_type x, argument_type y) { return x >= y; }
};

template<typename argument_type>
struct less_than {
  using value_type = bool;
  static bool apply(argument_type x, argument_type y) { return x < y; }
};

template<typename argument_type>
struct less_equal_than {
  using value_type = bool;
  static bool apply(argument_type x, argument_type y) { return x <= y; }
};

template<typename argument_type>
struct equal {
  using value_type = bool;
  static bool apply(argument_type x, argument_type y) { return x == y; }
};

template<typename argument_type>
struct not_equal {
  using value_type = bool;
  static bool apply(argument_type x, argument_type y) { return x != y; }
};


struct logical_and {
  using value_type = bool;
  static bool apply(bool x, bool y) { return x and y; }
};

struct logical_or {
  using value_type = bool;
  static bool apply(bool x, bool y) { return x or y; }
};

#endif /* _OPERATOR_H_ */
