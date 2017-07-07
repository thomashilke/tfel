#ifndef _INTEGRAND_EXPRESSION_H_
#define _INTEGRAND_EXPRESSION_H_

#include <cstddef>


template<std::size_t n, typename ... Ts> struct get_type;

template<std::size_t n, typename T, typename ... Ts>
struct get_type<n, T, Ts...> { typedef typename get_type<n - 1, Ts...>::type type; };

template<typename T, typename ... Ts>
struct get_type<0, T, Ts...> { typedef T type; };


template<std::size_t n, typename ... Ts> struct get_arg;

template<std::size_t n, typename T, typename ... Ts>
struct get_arg<n, T, Ts...> {
  static
  typename get_type<n, T, Ts...>::type
  get(T t, Ts ... ts) { return get_arg<n - 1, Ts...>::get(ts...); }
};

template<typename T, typename ... Ts>
struct get_arg<0, T, Ts...> {
  static T get(T t, Ts ... ts) { return t; }
};

template<typename T>
struct get_arg<0, T> {
  static T get(T t) { return t; }
};


template<std::size_t n, typename ... Ts>
typename get_type<n, Ts...>::type
argument(Ts... ts) {
  return get_arg<n, Ts...>::get(ts...);
}


template<std::size_t arg, std::size_t order, std::size_t derivative>
struct form {
  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return argument<arg>(ts...)[derivative];
  }
};

struct free_function {
  free_function(double (*f)(const double*)): f(f) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return f(x);
  }
  
private:
  typedef double (*function_type)(const double*);
  function_type f;
};

struct constant {
  constant(double c): value(c) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return value;
  }
  
private:
  const double value;
};


struct add { static double apply(double x, double y) { return x + y; } };
struct substract { static double apply(double x, double y) { return x - y; } };
struct multiply { static double apply(double x, double y) { return x * y; } };


template<typename expr_t>
struct expression {
  expr_t expr;
  expression(const expr_t& e): expr(e) {}

  template<typename ... Ts>
  double operator()(unsigned int k, const double* x, const double* x_hat,
		    Ts ... ts) const {
    return expr(k, x, x_hat, ts...);
  }
};

template<typename left, typename right, typename op>
struct binary_expression {
  binary_expression(const expression<left>& l, const expression<right>& r): l(l), r(r) {}

  template<typename ... Ts>
  double operator()(unsigned int k, const double* x, const double* x_hat,
		    Ts ... ts) const {
    return op::apply(l(k, x, x_hat, ts...),
		     r(k, x, x_hat, ts...));
  }
  
  expression<left> l;
  expression<right> r;
};


typedef expression<form<1, 0, 0> > trial_function_t;
typedef expression<form<0, 0, 0> > test_function_t;
typedef expression<free_function> free_function_t;
typedef expression<constant> constant_t;


template<std::size_t d, typename expression>
struct differentiate;

template<std::size_t d, std::size_t arg>
struct differentiate<d, form<arg, 0, 0> > {
  typedef form<arg, 1, d> type;
  
  static
  type initialize(const expression<form<arg, 0, 0> >&) { return type(); }
};

template<std::size_t d, typename left, typename right>
struct differentiate<d, binary_expression<left, right, multiply> > {
  typedef binary_expression<typename differentiate<d, left>::type, right, multiply> first_type;
  typedef binary_expression<left, typename differentiate<d, right>::type, multiply> second_type;
  typedef binary_expression<first_type, second_type, add> type;

  static
  type initialize(const binary_expression<left, right, multiply>& expr) {
    const first_type t1(differentiate<d, left>::initialize(expr.l), expr.r);
    const second_type t2(expr.l, differentiate<d, right>::initialize(expr.r));
    return type(t1, t2);
  }
};

template<std::size_t d, typename left, typename right, typename op>
struct differentiate<d, binary_expression<left, right, op> > {
  typedef binary_expression<typename differentiate<d, left>::type,
			    typename differentiate<d, right>::type, op>
  type;
  
  static
  type initialize(const binary_expression<left, right, op>& expr) {
    return type(differentiate<d, left>::initialize(expr.l),
		differentiate<d, right>::initialize(expr.r));
  }
};


template<std::size_t direction, typename expr>
expression<typename differentiate<direction, expr>::type>
d(const expression<expr>& e) {
  return differentiate<direction, expr>::initialize(e.expr);
};


template<typename left, typename right>
expression<binary_expression<left, right, multiply> >
operator*(const expression<left>& l, const expression<right>& r) {
  typedef expression<binary_expression<left, right, multiply> > ret;
  return ret(binary_expression<left, right, multiply>(l, r));
}

template<typename expr>
expression<binary_expression<expression<constant>, expr, multiply> >
operator*(double c, const expression<expr>& e) {
  typedef expression<binary_expression<expression<constant>, expr, multiply> > ret;
  return ret(binary_expression<expression<constant>,
	     expr,
	     multiply>(expression<constant>(constant(c)), e));
}

template<typename expr>
expression<binary_expression<expression<free_function>, expr, multiply> >
operator*(double(*f)(const double*), const expression<expr>& e) {
  typedef expression<binary_expression<expression<free_function>, expr, multiply> > ret;
  return ret(binary_expression<expression<free_function>,
	     expr,
	     multiply>(expression<free_function>(free_function(f)), e));
}

template<typename left, typename right>
expression<binary_expression<left, right, add> >
operator+(const expression<left>& l, const expression<right>& r) {
  typedef expression<binary_expression<left, right, add> > ret;
  return ret(binary_expression<left, right, add>(l, r));
}

template<typename left, typename right>
expression<binary_expression<left, right, substract> >
operator-(const expression<left>& l, const expression<right>& r) {
  typedef expression<binary_expression<left, right, substract> > ret;
  return ret(binary_expression<left, right, substract>(l, r));
}


#endif /* _INTEGRAND_EXPRESSION_H_ */
