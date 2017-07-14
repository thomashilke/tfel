#ifndef _PROJECTOR_H_
#define _PROJECTOR_H_

#include "form.hpp"
#include "timer.hpp"

namespace projector {

  template<typename fe_type, typename quadrature_type>
  typename finite_element_space<fe_type>::element
  l2(double (*fun)(const double*), const finite_element_space<fe_type>& fes) {
    timer t;
    typedef finite_element_space<fe_type> fes_type;
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      a += integrate<quadrature_type>(u * v, fes.get_mesh());
    }
    std::cout << std::setw(40) << std::right << "projector::l2::bilinear form: " << t.tic() << " [ms]\n";
    linear_form<fes_type> f(fes); {
      auto v(f.get_test_function());
      f += integrate<quadrature_type>(fun * v, fes.get_mesh());
    }
    std::cout << std::setw(40) << std::right << "projector::l2::linear form: " << t.tic() << " [ms]\n";

    return a.solve(f);
  }

  template<typename fe_type>
  typename finite_element_space<fe_type>::element
  lagrange(double (*fun)(const double*), const finite_element_space<fe_type>& fes) {
    typename finite_element_space<fe_type>::element fun_h;

    for (std::size_t n(0); n < fes.get_dof_number(); ++n) {
      // TODO
    }
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