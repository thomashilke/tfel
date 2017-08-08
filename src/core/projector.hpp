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
    static_assert(expr_t::rank == 0, "");
    
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

  template<typename fe_type, typename expr_t>
  typename finite_element_space<fe_type>::element
  lagrange(const expression<expr_t>& expr, const finite_element_space<fe_type>& fes) {
    using cell_type = typename fe_type::cell_type;
    
    static_assert(expr_t::rank == 0, "");
    static_assert(fe_type::is_lagrangian,
		  "lagrange projector is only defined for lagrangian finite element space.");

    const mesh<cell_type>& m(fes.get_mesh());
    
    array<double> coefficients{fes.get_dof_number()};
    for (std::size_t i(0); i < fes.get_dof_number(); ++i) {
      const std::size_t i_hat(fes.get_dof_local_id(i));
      const std::size_t k(fes.get_dof_element(i));

      array<double> x_hat{1, cell_type::n_dimension};
      std::copy(
	&fe_type::x[i_hat][0],
	&fe_type::x[i_hat][0] + cell_type::n_dimension,
	&x_hat.at(0, 0));

      array<double> x(cell_type::map_points_to_space_coordinates(
        m.get_vertices(),
	m.get_elements(),
	k,
	x_hat));

      expr.prepare(k, &x.at(0,0), &x_hat.at(0,0));
      coefficients.at(i) = expr(k, &x.at(0,0), &x_hat.at(0,0));
    }
    
    return typename finite_element_space<fe_type>::element(fes, coefficients);
  }

}

#endif /* _PROJECTOR_H_ */
