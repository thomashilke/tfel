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
    //timer t;
    typedef finite_element_space<fe_type> fes_type;
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      a += integrate<quadrature_type>(u * v, fes.get_mesh());
    }
    //std::cout << std::setw(40) << std::right << "projector::l2::bilinear form: " << t.tic() << " [ms]\n";
    linear_form<fes_type> f(fes); {
      auto v(f.get_test_function());
      f += integrate<quadrature_type>(expr * v, fes.get_mesh());
    }
    //std::cout << std::setw(40) << std::right << "projector::l2::linear form: " << t.tic() << " [ms]\n";

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
    array<double> coefficients{fes.get_dof_number()};

    for (std::size_t n(0); n < fes.get_dof_number(); ++n) {
      const auto x(fes.get_dof_space_coordinate(n));
      coefficients.at(n) = fun(&x.at(0, 0));
    }

    return typename finite_element_space<fe_type>::element(fes, coefficients);
  }

  /* 
  // This cannot work since the return value hold a reference to the fes.
  template<typename fe_type, typename quadrature_type>
  typename finite_element_space<fe_type>::element
  l2(double (*fun)(const double*), const mesh<typename fe_type::cell_type>& m) {
    typedef finite_element_space<fe_type> fes_type;
    fes_type fes(m);
    return l2<fe_type, quadrature_type>(fun, fes);
  }
  */
}

#endif /* _PROJECTOR_H_ */
