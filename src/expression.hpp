#ifndef _INTEGRAND_EXPRESSION_H_
#define _INTEGRAND_EXPRESSION_H_

#include <cstddef>


template<std::size_t order, std::size_t derivative>
struct test_function {
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    const double* phi, const double* psi) const {
    return phi[derivative];
  }
};

template<std::size_t order, std::size_t derivative>
struct trial_function {
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    const double* phi, const double* psi) const {
    return psi[derivative];
  }
};

struct free_function {
  free_function(double (*f)(const double*)): f(f) {}
  
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
	      const double* phi, const double* psi) const {
    return f(x);
  }
  
private:
  typedef double (*function_type)(const double*);
  function_type f;
};

struct constant {
  constant(double c): value(c) {}
  
  double operator()(unsigned int k,
		    const double* x, const double* x_hat,
		    const double* phi, const double* psi) const {
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
  double operator()(unsigned int k, const double* x, const double* x_hat,
		    const double* phi, const double* psi) const {
    return expr(k, x, x_hat, phi, psi);
  }
};

template<typename left, typename right, typename op>
struct binary_expression {
  binary_expression(const expression<left>& l, const expression<right>& r): l(l), r(r) {}

  double operator()(unsigned int k, const double* x, const double* x_hat,
		    const double* phi, const double* psi) const {
    return op::apply(l(k, x, x_hat, phi, psi),
		     r(k, x, x_hat, phi, psi));
  }
  
  expression<left> l;
  expression<right> r;
};


typedef expression<test_function<0, 0> > test_function_t;
typedef expression<trial_function<0, 0> > trial_function_t;
typedef expression<free_function> free_function_t;
typedef expression<constant> constant_t;


template<std::size_t d, typename expression>
struct differentiate;

template<std::size_t d>
struct differentiate<d, test_function<0, 0> > {
  typedef test_function<1, d> type;
  
  static
  type initialize(const expression<test_function<0, 0> >&) { return type(); }
};

template<std::size_t d>
struct differentiate<d, trial_function<0, 0> > {
  typedef trial_function<1, d> type;
  
  static
  type initialize(const expression<trial_function<0, 0> >&) { return type(); }
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
