#ifndef _COMPOSITE_LINEAR_FORM_H_
#define _COMPOSITE_LINEAR_FORM_H_

template<typename te_cfe_type>
class linear_form<composite_finite_element_space<te_cfe_type> > {
public:
  using test_cfes_type = composite_finite_element_space<te_cfe_type>;
  using test_cfe_type = typename test_cfes_type::cfe_type;
  using test_fe_list = typename test_cfe_type::fe_list;
  using fe_list = test_fe_list;
  using unique_fe_list = unique_t<fe_list>;
  using fe_cell_type = typename get_element_at_t<0, unique_fe_list>::cell_type;

  static const std::size_t n_test_component = test_cfe_type::n_component;

  linear_form(const test_cfes_type& te_cfes)
    : test_cfes(te_cfes),
      f{te_cfes.get_total_dof_number()} {
    f.fill(0.0);

    std::size_t test_global_dof_number[n_test_component];
    fill_array_with_return_values<std::size_t,
				  get_dof_number_impl<test_cfes_type>,
				  0,
				  n_test_component>::template fill<const test_cfes_type&>(&test_global_dof_number[0],
											  test_cfes);

    test_global_dof_offset = std::vector<std::size_t>(n_test_component, 0ul);
    std::partial_sum(test_global_dof_number,
		     test_global_dof_number + n_test_component - 1,
		     test_global_dof_offset.begin() + 1);
  }


  template<typename B_INFO>
  struct evaluate_block {
    using T = get_element_at_t<0, B_INFO>;
    typedef typename T::quadrature_type quadrature_type;

    static const std::size_t m = get_element_at_t<1, B_INFO>::value;

    using test_fe_type = get_element_at_t<m, typename test_cfe_type::fe_list>;
    using linear_form_type = linear_form<composite_finite_element_space<te_cfe_type> >;

    
    static void call(linear_form_type& linear_form,
		     const get_element_at_t<0, B_INFO>& integration_proxy,
		     const std::size_t k,
		     const array<double>& omega,
		     const array<double>& xq,
		     const array<double>& xq_hat,
		     const fe_value_manager<unique_fe_list>& fe_values,
		     const fe_value_manager<unique_fe_list>& fe_zvalues) {
      const double* psi[n_test_component];
      
      const std::size_t n_q(quadrature_type::n_point);
      const std::size_t n_test_dof(test_fe_type::n_dof_per_element);

      // evaluate the weak form
      const double volume(integration_proxy.m.get_cell_volume(k));
      for (unsigned int i(0); i < n_test_dof; ++i) {
	double rhs_el(0.0);
	for (unsigned int q(0); q < n_q; ++q) {
	  select_function_valuation<test_fe_list, m, unique_fe_list>(psi, q, i,
								     fe_values, fe_zvalues);
	  integration_proxy.f.prepare(k, &xq.at(q, 0), &xq_hat.at(q, 0));
	  rhs_el += volume * omega.at(q)
	    * expression_call_wrapper<0, n_test_component>::call(integration_proxy.f, psi,
								 k, &xq.at(q, 0), &xq_hat.at(q, 0));
	}
	linear_form.f.at(linear_form.test_cfes.template get_dof<m>(integration_proxy.get_global_element_id(k), i)
			 + linear_form.test_global_dof_offset[m]) += rhs_el;
      }
    }
  };

  
  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 1, "linear_form expects rank-2 expression.");

    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;
    using form_type = typename T::form_type;

    // m is the mesh over which we integrate
    const auto& m(integration_proxy.m);

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the quadrature points;
    array<double> xq_hat{n_q, fe_cell_type::n_dimension};
    array<double> xq{n_q, fe_cell_type::n_dimension};
    
    // storage for the point-wise basis function evaluation
    fe_value_manager<unique_fe_list> fe_values(n_q), fe_zvalues(n_q);
    fe_zvalues.clear();

    if (T::point_set_number == 1) {
      xq_hat = integration_proxy.get_quadrature_points(0);
      fe_values.set_points(xq_hat);
    }
    
    for (unsigned int k(0); k < m.get_element_number(); ++k) {

      // prepare the quadrature points
      if (T::point_set_number > 1) {
	xq_hat = integration_proxy.get_quadrature_points(k);
	fe_values.set_points(xq_hat);
      }

      if (form_type::require_space_coordinates)
	cell_type::map_points_to_space_coordinates(xq, m.get_vertices(),
						   m.get_elements(),
						   k, xq_hat);
      
      // prepare the basis function values
      if (form_type::differential_order > 0) {
	const array<double> jmt(m.get_jmt(k));
	fe_values.prepare(jmt);
      }


      /*
       *  Compile time loop over all the blocks
       */
      using test_blocks_il = wrap_t<type_list, make_integral_list_t<std::size_t, n_test_component> >;
      using block_info = append_to_each_element_t<T, test_blocks_il>;

      call_for_each<evaluate_block, block_info>::call(*this, integration_proxy,
						      k,
						      omega,
						      xq, xq_hat,
						      fe_values, fe_zvalues);
    }
  }

  template<std::size_t n>
  expression<form<n, 1, 0> > get_test_function() const {
    return form<n, 1, 0>();
  }

  void clear() {
    f.fill(0.0);
  }
  
  void show(std::ostream& stream) 
  {
    stream << "rhs = [" << f.at(0);
    for (std::size_t j(1); j < f.get_size(0); ++j)
      stream << "; " << f.at(j);
    stream << "];" << std::endl;
  }

  const array<double>& get_coefficients() const { return f; }
    
private:
  const test_cfes_type& test_cfes;
  std::vector<std::size_t> test_global_dof_offset;

  array<double> f;
};

#endif /* _COMPOSITE_LINEAR_FORM_H_ */
