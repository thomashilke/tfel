#ifndef _PROJECTOR_H_
#define _PROJECTOR_H_

#include <functional>

#include "fe.hpp"
#include "quadrature.hpp"
#include "form.hpp"
#include "timer.hpp"

namespace projector {

  template<typename fe_type, typename quadrature_type, typename expr_t>
  typename finite_element_space<fe_type>::element
  l2(const expression<expr_t>& expr, const finite_element_space<fe_type>& fes) {

    typedef finite_element_space<fe_type> fes_type;
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      a += integrate<quadrature_type>(u * v, fes.get_mesh());
    }

    linear_form<fes_type> f(fes); {
      auto v(f.get_test_function());
      f += integrate<quadrature_type>(expr * v, fes.get_mesh());
    }

    return a.solve(f);
  }
  
  template<typename fe_type, typename quadrature_type>
  typename finite_element_space<fe_type>::element
  l2(const std::function<double(const double*)>& fun, const finite_element_space<fe_type>& fes) {
    return l2<fe_type, quadrature_type>(make_expr(fun), fes);
  }

  template<typename fe_type>
  typename finite_element_space<fe_type>::element
  lagrange(const std::function<double(const double*)>& fun, const finite_element_space<fe_type>& fes) {
    static_assert(fe_type::is_lagrangian,
		  "lagrange projector is only defined for lagrangian finite element space.");
    
    array<double> coefficients{fes.get_dof_number()};

    for (std::size_t n(0); n < fes.get_dof_number(); ++n) {
      const auto x(fes.get_dof_space_coordinate(n));
      coefficients.at(n) = fun(&x.at(0, 0));
    }

    return typename finite_element_space<fe_type>::element(fes, coefficients);
  }

}

#endif /* _PROJECTOR_H_ */
