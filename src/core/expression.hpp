#ifndef _INTEGRAND_EXPRESSION_H_
#define _INTEGRAND_EXPRESSION_H_

#include <cstddef>
#include <algorithm>
#include <functional>

#include "meta.hpp"
#include "fes.hpp"


template<std::size_t arg, std::size_t rnk, std::size_t derivative>
struct form {
  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    static_assert(sizeof...(Ts) >= arg, "");
    return argument<arg>(ts...)[derivative];
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const  {}

  static constexpr std::size_t rank = rnk;
  static constexpr std::size_t differential_order = derivative == 0 ? 0 : 1;
  static constexpr bool require_space_coordinates = false;
};

struct free_function {
  free_function(double (*f)(const double*)): f(f) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return cached_value;
  }

  double operator()(unsigned int k,
		    const double* x, const double* x_hat) const {
    return cached_value;
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {
    cached_value = f(x);
  }

  static constexpr std::size_t rank = 0;
  static constexpr std::size_t differential_order = 0;
  static constexpr bool require_space_coordinates = true;
  
private:
  typedef double (*function_type)(const double*);
  function_type f;
  mutable double cached_value = 0.0;
};


struct std_function {
  std_function(const std::function<double(const double*)>& f): f(f) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return cached_value;
  }

  double operator()(unsigned int k,
		    const double* x, const double* x_hat) const {
    return cached_value;
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {
    cached_value = f(x);
  }

  static constexpr std::size_t rank = 0;
  static constexpr std::size_t differential_order = 0;
  static constexpr bool require_space_coordinates = true;

private:
  std::function<double(const double*)> f;
  mutable double cached_value = 0.0;
};

template<typename fe>
struct finite_element_function {
public:
  finite_element_function(const typename finite_element_space<fe>::element& v): v(v) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return cached_value;
  }

  double operator()(unsigned int k,
		    const double* x, const double* x_hat) const {
    return cached_value;
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {
    cached_value = v.evaluate(k, x_hat);
  }

  static constexpr std::size_t rank = 0;
  static constexpr std::size_t differential_order = 0;
  static constexpr bool require_space_coordinates = false;
			     
private:
  const typename finite_element_space<fe>::element& v;
  mutable double cached_value = 0.0;
};

struct constant {
  constant(double c): value(c) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    Ts ... ts) const {
    return value;
  }

  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat) const {
    return value;
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {}

  static constexpr std::size_t rank = 0;
  static constexpr std::size_t differential_order = 0;
  static constexpr bool require_space_coordinates = false;
  
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

  double operator()(unsigned int k, const double* x, const double* x_hat) const {
    return expr(k, x, x_hat);
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {
    expr.prepare(k, x, x_hat);
  }

  static constexpr std::size_t rank = expr_t::rank;
  static constexpr std::size_t differential_order = expr_t::differential_order;
  static constexpr bool require_space_coordinates = expr_t::require_space_coordinates;
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

  double operator()(unsigned int k, const double* x, const double* x_hat) const {
    return op::apply(l(k, x, x_hat),
		     r(k, x, x_hat));
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {
    l.prepare(k, x, x_hat);
    r.prepare(k, x, x_hat);
  }
  
  static constexpr std::size_t rank = left::rank < right::rank ? right::rank : left::rank;
  static constexpr std::size_t differential_order = left::differential_order < right::differential_order ? right::differential_order : left::differential_order;

  static constexpr bool require_space_coordinates = left::require_space_coordinates or right::require_space_coordinates;
  
  expression<left> l;
  expression<right> r;
};

template<typename inner_expr>
struct composition {
  composition(double (*f)(double), const expression<inner_expr>& e): f(f), e(e) {}

  template<typename ... Ts>
  double operator()(unsigned int k, const double* x, const double* x_hat,
		    Ts ... ts) const {
    return f(e(k, x, x_hat, ts...));
  }

  double operator()(unsigned int k, const double* x, const double* x_hat) const {
    return f(e(k, x, x_hat));
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {
    e.prepare(k, x, x_hat);
  }

  static constexpr std::size_t rank = inner_expr::rank;
  static constexpr std::size_t differential_order = inner_expr::differential_order;
  static constexpr bool require_space_coordinates = inner_expr::require_space_coordinates;
  
  double (*f)(double);
  expression<inner_expr> e;
};

template<typename fe>
using fe_function_t = expression<finite_element_function<fe> >;

typedef expression<free_function> free_function_t;
typedef expression<constant> constant_t;
typedef expression<std_function> std_function_t;

template<typename fe>
fe_function_t<fe> make_expr(const typename finite_element_space<fe>::element& u) {
  return fe_function_t<fe>(finite_element_function<fe>(u) );
}
free_function_t make_expr(double(*f)(const double*)) {
  return free_function_t(free_function(f));
}
constant_t make_expr(double c) {
  return constant_t(constant(c));
}

std_function_t make_expr(const std::function<double(const double*)>& f) {
  return std_function(f);
}

template<std::size_t d, typename expression>
struct differentiate;

template<std::size_t d, std::size_t arg, std::size_t rnk>
struct differentiate<d, form<arg, rnk, 0> > {
  typedef form<arg, rnk, d> type;
  
  static
  type initialize(const expression<form<arg, rnk, 0> >&) { return type(); }
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


template<typename inner_expr>
expression<composition<inner_expr> > compose(double (*f)(double), const expression<inner_expr>& e) {
  return composition<inner_expr>(f, e.expr);
}

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


template<std::size_t n, std::size_t n_max>
struct expression_call_wrapper {
  template<typename form_t, typename ... As>
  static double call(form_t form, const double** arg_array, As ... as) {
    return expression_call_wrapper<n + 1, n_max>::template call<form_t,
								As...,
								const double*>(form,
									       arg_array + 1,
									       as...,
									       *arg_array);
  }
};

template<std::size_t n_max>
struct expression_call_wrapper<n_max, n_max> {
  template<typename form_t, typename ... As>
  static double call(form_t form, const double** arg_array, As ... as) {
    return form(as...);
  }
};

#endif /* _INTEGRAND_EXPRESSION_H_ */
