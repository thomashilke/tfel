#ifndef _LINEAR_FORM_H_
#define _LINEAR_FORM_H_

template<typename test_fes_type>
class linear_form {
public:
  linear_form(const test_fes_type& te_fes)
    : test_fes(te_fes), f{te_fes.get_dof_number()} {
    f.fill(0.0);
  }

  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 1, "linear_form expects rank-1 expression.");
    
    typedef typename test_fes_type::fe_type test_fe_type;
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    const auto& m(integration_proxy.m);
    const std::size_t dim(m.get_embedding_space_dimension());

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);


    // storage for the point-wise basis function evaluation
    const std::size_t n_dof(test_fe_type::n_dof_per_element);
    array<double> psi{n_q, n_dof, dim + 1};

    
    // loop over the elements
    for (unsigned int k(0); k < m.get_element_number(); ++k) {
      const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
      const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
									m.get_elements(),
									k, xq_hat));

      // prepare the basis function on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i)
	  psi.at(q, i, 0) = test_fe_type::phi(i, &xq_hat.at(q, 0));
      }
      
      // prepare the basis function derivatives on the quadrature points
      const array<double> jmt(m.get_jmt(k));
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i) {
	  for (std::size_t n(0); n < dim; ++n) {
	    psi.at(q, i, 1 + n) = 0.0;
	    for (std::size_t k(0); k < dim; ++k)
	      psi.at(q, i, 1 + n) += jmt.at(n, k) * test_fe_type::dphi(k, i, &xq_hat.at(q, 0));
	  }
	}
      }

      // evaluate the weak form
      const double volume(m.get_cell_volume(k));
      for (unsigned int j(0); j < n_dof; ++j) {
	double rhs_el(0.0);
	for (unsigned int q(0); q < n_q; ++q) {
	  rhs_el += volume * omega.at(q) * integration_proxy.f(k, &xq.at(q, 0),
							       &xq_hat.at(q, 0),
							       &psi.at(q, j, 0));
	}
	f.at(test_fes.get_dof(integration_proxy.get_global_element_id(k), j)) += rhs_el;
      }
    }
  }

  expression<form<0,1,0> > get_test_function() const { return form<0,1,0>(); }

  const array<double>& get_coefficients() const { return f; }

  void clear() {
    std::fill(f.get_data(),
	      f.get_data() + f.get_size(0),
	      0.0);
  }
  
  void show(std::ostream& stream) {
    stream << "rhs = [" << f.at(0);
    for (std::size_t j(1); j < f.get_size(0); ++j)
      stream << "; " << f.at(j);
    stream << "];" << std::endl;
  }
private:
  const test_fes_type& test_fes;
  
  array<double> f;
};

#endif /* _LINEAR_FORM_H_ */
