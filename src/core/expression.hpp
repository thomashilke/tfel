#ifndef _INTEGRAND_EXPRESSION_H_
#define _INTEGRAND_EXPRESSION_H_

#include <cstddef>
#include <algorithm>
#include <functional>

#include "meta.hpp"
#include "fes.hpp"
#include "operator.hpp"
#include "fe_value_manager.hpp"


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

template<typename fe, std::size_t d = 0>
struct finite_element_function {
  template<typename fe_p, std::size_t d_p>
  friend struct finite_element_function;
  
public:
  finite_element_function(const typename finite_element_space<fe>::element& v)
    : v(v), fe_values(1), xq_hat{1, fe::cell_type::n_dimension} {}

  template<std::size_t dd>
  finite_element_function(const finite_element_function<fe, dd>& f)
    : v(f.v), fe_values(1), xq_hat{1, fe::cell_type::n_dimension} {}

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
    std::copy(x_hat, x_hat + xq_hat.get_size(1), xq_hat.get_data());
    fe_values.set_points(xq_hat);
    if (d > 0)
      fe_values.prepare(v.get_finite_element_space().get_mesh().get_jmt(k));

    const std::size_t n_dof(fe::n_dof_per_element);
    const array<double>& phi(fe_values.template get_values<0>());
    const finite_element_space<fe>& fes(v.get_finite_element_space());
    const array<double>& coefficients(v.get_coefficients());

    cached_value = 0.0;
    for (unsigned int i(0); i < n_dof; ++i)
      cached_value += coefficients.at(fes.get_dof(k, i)) * phi.at(0, i, d);
  }

  static constexpr std::size_t rank = 0;
  static constexpr std::size_t differential_order = d;
  static constexpr bool require_space_coordinates = false;
			     
private:
  const typename finite_element_space<fe>::element& v;
  mutable fe_value_manager<type_list<fe> > fe_values;
  mutable array<double> xq_hat;
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

template<typename value_t, typename mesh_t>
struct mesh_data_component {
  mesh_data_component(const mesh_data<value_t, mesh_t>& data, std::size_t c):
    data(data), c(c) {}

  template<typename ... Ts>
  double operator()(unsigned int k,
                    const double* x,
                    const double* x_hat,
                    Ts ... ts) const {
    return data.evaluate(k, x_hat, c);
  }
  
  template<typename ... Ts>
  double operator()(unsigned int k,
		    const double* x, const double* x_hat) const {
    return data.evaluate(k, x_hat, c);
  }

  void prepare(unsigned int k, const double* x, const double* x_hat) const {}

  static constexpr std::size_t rank = 0;
  static constexpr std::size_t differential_order = 0;
  static constexpr bool require_space_coordinates = false;

private:
  const mesh_data<value_t, mesh_t>& data;
  std::size_t c;
};


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

template<typename value_t, typename mesh_t>
using mesh_data_component_t = expression<mesh_data_component<value_t, mesh_t> >;

template<typename fe>
fe_function_t<fe>
make_expr(const typename finite_element_space<fe>::element& u) {
  return fe_function_t<fe>(finite_element_function<fe>(u) );
}

free_function_t
make_expr(double(*f)(const double*)) {
  return free_function_t(free_function(f));
}

constant_t
make_expr(double c) {
  return constant_t(constant(c));
}

std_function_t
make_expr(const std::function<double(const double*)>& f) {
  return std_function(f);
}

template<typename value_t, typename mesh_t>
mesh_data_component_t<value_t, mesh_t>
make_expr(const mesh_data<value_t, mesh_t>& data, std::size_t c) {
  return mesh_data_component_t<value_t, mesh_t> (mesh_data_component<value_t, mesh_t>(data, c));
}

template<std::size_t d, typename expression>
struct differentiate;

template<std::size_t d, std::size_t arg, std::size_t rnk>
struct differentiate<d, form<arg, rnk, 0> > {
  typedef form<arg, rnk, d> type;
  
  static
  type initialize(const expression<form<arg, rnk, 0> >&) { return type(); }
};

template<std::size_t d, typename fe>
struct differentiate<d, finite_element_function<fe, 0> > {
  typedef finite_element_function<fe, d> type;

  static
  type initialize(const finite_element_function<fe, 0> & e) {
    return type(e);
  }
};

template<std::size_t d, typename left, typename right>
struct differentiate<d, binary_expression<left, right, multiply<double>> > {
  typedef binary_expression<typename differentiate<d, left>::type, right, multiply<double>> first_type;
  typedef binary_expression<left, typename differentiate<d, right>::type, multiply<double>> second_type;
  typedef binary_expression<first_type, second_type, add<double>> type;

  static
    type initialize(const binary_expression<left, right, multiply<double>>& expr) {
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
expression<binary_expression<left, right, multiply<double>> >
operator*(const expression<left>& l, const expression<right>& r) {
  typedef expression<binary_expression<left, right, multiply<double>> > ret;
  return ret(binary_expression<left, right, multiply<double>>(l, r));
}

template<typename expr>
expression<binary_expression<expression<constant>, expr, multiply<double>> >
operator*(double c, const expression<expr>& e) {
  typedef expression<binary_expression<expression<constant>, expr, multiply<double>> > ret;
  return ret(binary_expression<expression<constant>,
	     expr,
	     multiply<double>>(expression<constant>(constant(c)), e));
}

template<typename expr>
expression<binary_expression<expression<free_function>, expr, multiply<double>> >
operator*(double(*f)(const double*), const expression<expr>& e) {
  typedef expression<binary_expression<expression<free_function>, expr, multiply<double>> > ret;
  return ret(binary_expression<expression<free_function>,
	     expr,
	     multiply<double>>(expression<free_function>(free_function(f)), e));
}

template<typename left, typename right>
expression<binary_expression<left, right, add<double>> >
operator+(const expression<left>& l, const expression<right>& r) {
  typedef expression<binary_expression<left, right, add<double>> > ret;
  return ret(binary_expression<left, right, add<double>>(l, r));
}

template<typename left, typename right>
expression<binary_expression<left, right, substract<double>> >
operator-(const expression<left>& l, const expression<right>& r) {
  typedef expression<binary_expression<left, right, substract<double>> > ret;
  return ret(binary_expression<left, right, substract<double>>(l, r));
}


template<std::size_t n, std::size_t n_max>
struct expression_call_wrapper {
  template<typename form_t, typename ... As>
  static double call(const form_t& form, const double** arg_array, As&& ... as) {
    // leave the type deduction of the call method to the compiler...
    return expression_call_wrapper<n + 1, n_max>::template call</*form_t,
								As...,
								const double**/>(form,
									       arg_array + 1,
									       std::forward<As>(as)...,
									       *arg_array);
  }
};

template<std::size_t n_max>
struct expression_call_wrapper<n_max, n_max> {
  template<typename form_t, typename ... As>
  static double call(const form_t& form, const double** arg_array, As&& ... as) {
    return form(std::forward<As>(as)...);
  }
};

#endif /* _INTEGRAND_EXPRESSION_H_ */
