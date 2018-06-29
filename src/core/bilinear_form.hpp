#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

enum class algebraic_block {test_block, trial_block};

template<typename test_fes_type, typename trial_fes_type>
class bilinear_form {
public:
  bilinear_form(const test_fes_type& te_fes,
		const trial_fes_type& tr_fes,
                std::size_t algebraic_equation_number = 0,
                std::size_t algebraic_dof_number = 0)
    : test_fes(te_fes), trial_fes(tr_fes),
      a(te_fes.get_dof_number() + algebraic_equation_number,
        tr_fes.get_dof_number() + algebraic_dof_number),
      a_eq_number(algebraic_equation_number),
      a_dof_number(algebraic_dof_number) {
    clear();
  }

  ~bilinear_form() {}

  template<typename T>
  void operator+=(T integration_proxy) {
    static_assert(T::form_type::rank == 2, "linear_form expects rank-2 expression.");
      
    typedef typename test_fes_type::fe_type test_fe_type;
    typedef typename trial_fes_type::fe_type trial_fe_type;
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;
    using form_type = typename T::form_type;

    const auto& m(integration_proxy.m);

    
    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the quadrature points
    array<double> xq_hat{n_q, test_fe_type::cell_type::n_dimension};
    array<double> xq{n_q, test_fe_type::cell_type::n_dimension};
    
    // storage for the point-wise basis function evaluation
    using fe_list = type_list<test_fe_type, trial_fe_type>;
    using unique_fe_list = unique_t<fe_list>;

    const std::size_t test_fe_index(get_index_of_element<test_fe_type, unique_fe_list>::value);
    const std::size_t trial_fe_index(get_index_of_element<trial_fe_type, unique_fe_list>::value);
    fe_value_manager<unique_fe_list> fe_values(n_q);
    
    if (T::point_set_number == 1) {
      xq_hat = integration_proxy.get_quadrature_points(0);
      fe_values.set_points(xq_hat);
      xq.fill(0.0);
    }

    const std::size_t n_test_dof(test_fe_type::n_dof_per_element);
    const std::size_t n_trial_dof(trial_fe_type::n_dof_per_element);
    array<double> a_el{n_test_dof, n_trial_dof};

    // loop over the elements
    for (unsigned int k(0); k < m.get_cell_number(); ++k) {
      a_el.fill(0.0);

      // prepare the quadrature points if necessary
      if (T::point_set_number > 1) {
	xq_hat = integration_proxy.get_quadrature_points(k);
	fe_values.set_points(xq_hat);
      }

      if (form_type::require_space_coordinates)
	cell_type::map_points_to_space_coordinates(xq,
						   m.get_vertices(),
						   m.get_cells(),
						   k, xq_hat);

      if (form_type::differential_order == 1ul) {
	// prepare the basis function values
	const array<double>& jmt(m.get_jmt(k));
	fe_values.prepare(jmt); 
      }
      
      const array<double>& psi(fe_values.template get_values<test_fe_index>());
      const array<double>& phi(fe_values.template get_values<trial_fe_index>());

      // evaluate the weak form
      const double volume(m.get_cell_volume(k));
      for (unsigned int q(0); q < n_q; ++q) {
	integration_proxy.f.prepare(k, &xq.at(q, 0ul), &xq_hat.at(q, 0ul));
	
	for (unsigned int i(0); i < n_test_dof; ++i) {
	  for (unsigned int j(0); j < n_trial_dof; ++j) {
	    a_el.at(i, j) += omega.at(q) * integration_proxy.f(k,
							       &xq.at(q, 0), &xq_hat.at(q, 0),
							       &psi.at(q, i, 0),
							       &phi.at(q, j, 0));
	  }
	}
      }
      
      for (unsigned int i(0); i < n_test_dof; ++i)
	for (unsigned int j(0); j < n_trial_dof; ++j) {
	  accumulate(test_fes.get_dof(integration_proxy.get_global_cell_id(k), i),
		     trial_fes.get_dof(integration_proxy.get_global_cell_id(k), j),
		     a_el.at(i, j) * volume);
        }
    }
  }

  expression<form<0,1,0> > get_test_function() const { return form<0,1,0>(); }
  expression<form<1,2,0> > get_trial_function() const { return form<1,2,0>(); }

  template<algebraic_block block>
  class algebraic_block_handle;
  
  algebraic_block_handle<algebraic_block::test_block> algebraic_test_block(std::size_t a_dof_id);
  algebraic_block_handle<algebraic_block::trial_block> algebraic_trial_block(std::size_t a_eq_id);
  double& algebraic_block(std::size_t a_eq, std::size_t a_dof);
  
  typename trial_fes_type::element solve(const linear_form<test_fes_type>& form,
                                         solver::basic_solver& s,
                                         dictionary* result = nullptr) const {
    // Add the value of the dirichlet dof in the right hand side
    array<double> f{trial_fes.get_dof_number() + a_dof_number};
    std::copy(&form.get_coefficients().at(0),
	      &form.get_coefficients().at(0) + test_fes.get_dof_number(),
	      &f.at(0));
    std::copy(form.get_constraint_values().begin(),
	      form.get_constraint_values().end(),
	      &f.at(0) + test_fes.get_dof_number());
    
    for (const auto& i: test_fes.get_dirichlet_dof_values()) {
      const auto x(trial_fes.get_dof_space_coordinate(i.first));
      f.at(i.first) = i.second;
    }
    
    s.set_operator(a);
    array<double> x{trial_fes.get_dof_number() + a_dof_number};
    dictionary r;
    s.solve(f, x, r);

    if (result)
      *result = r;
    
    if (a_dof_number == 0) {
      return typename trial_fes_type::element(trial_fes, x);
    } else {
      array<double> coefficients{trial_fes.get_dof_number()};
      std::copy(&x.at(0),
		&x.at(0) + trial_fes.get_dof_number(),
		&coefficients.at(0));
      
      return typename trial_fes_type::element(trial_fes, coefficients);
    }
  }

  void clear() {
    a.clear();

    // Add the identity equations for each dirichlet dof
    for (const auto& i: test_fes.get_dirichlet_dof_values())
      a.add(i.first, i.first, 1.0);
  }
  
  void show(std::ostream& stream) const {
    //a.show(stream);
  }

  void export_data(std::ostream& stream) const {
    //a.export_data(stream);
  }
  
private:
  const test_fes_type& test_fes;
  const trial_fes_type& trial_fes;
  
  sparse_matrix a;
  std::size_t a_eq_number;
  std::size_t a_dof_number;

  void accumulate(std::size_t i, std::size_t j, double value) {
    if (test_fes.get_dirichlet_dof_values().count(i) == 0)
      a.add(i, j, value);
  }
};

template<typename test_fes_type, typename trial_fes_type>
template<algebraic_block block>
class bilinear_form<test_fes_type, trial_fes_type>::algebraic_block_handle {
  using bilinear_form_type = bilinear_form<test_fes_type, trial_fes_type>;
  using test_fe_type = typename test_fes_type::fe_type;
  using trial_fe_type = typename trial_fes_type::fe_type;

  using fe_list = type_list<test_fe_type, trial_fe_type>;
  using unique_fe_list = unique_t<fe_list>;

  static const std::size_t test_fe_index = get_index_of_element<test_fe_type, unique_fe_list>::value;
  static const std::size_t trial_fe_index = get_index_of_element<trial_fe_type, unique_fe_list>::value;

public:
  algebraic_block_handle(bilinear_form_type& f, std::size_t id)
    : b_form(f), id(id) {}

  std::size_t get_id() const { return id; }
  
  template<typename T>
  void operator+=(const T& integration_proxy) {
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    const auto& m(integration_proxy.m);

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the point-wise basis function evaluation
    fe_value_manager<unique_fe_list> fe_values(n_q), fe_zvalues(n_q);
    fe_zvalues.clear();

    if (block == algebraic_block::trial_block) {
      const std::size_t n_trial_dof(trial_fe_type::n_dof_per_element);
      array<double> a_el{n_trial_dof};
    
      // loop over the elements
      for (unsigned int k(0); k < m.get_cell_number(); ++k) {
        a_el.fill(0.0);
      
        // prepare the quadrature points
        const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
        const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
                                                                          m.get_cells(),
                                                                          k, xq_hat));

        // prepare the basis function values
        const array<double> jmt(m.get_jmt(k));
        fe_values.set_points(xq_hat);
        fe_values.prepare(jmt);
      
        const array<double>& psi(fe_zvalues.template get_values<test_fe_index>());
        const array<double>& phi(fe_values.template get_values<trial_fe_index>());

        // evaluate the weak form
        const double volume(m.get_cell_volume(k));
        for (unsigned int q(0); q < n_q; ++q) {
          integration_proxy.f.prepare(k, &xq.at(q, 0), &xq_hat.at(q, 0));
	
          for (unsigned int j(0); j < n_trial_dof; ++j) {
            a_el.at(j) += omega.at(q) * integration_proxy.f(k,
                                                            &xq.at(q, 0), &xq_hat.at(q, 0),
                                                            &psi.at(q, 0, 0),
                                                            &phi.at(q, j, 0));
          }
        }
        for (unsigned int j(0); j < n_trial_dof; ++j)
          b_form.accumulate(b_form.test_fes.get_dof_number() + id,
                            b_form.trial_fes.get_dof(integration_proxy.get_global_cell_id(k), j),
                            a_el.at(j) * volume);
      }

    } else if (block == algebraic_block::test_block) {
      const std::size_t n_test_dof(test_fe_type::n_dof_per_element);
      array<double> a_el{n_test_dof};

      // loop over the elements
      for (unsigned int k(0); k < m.get_cell_number(); ++k) {
        a_el.fill(0.0);
      
        // prepare the quadrature points
        const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
        const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
                                                                          m.get_cells(),
                                                                          k, xq_hat));
        // prepare the basis function values
        const array<double> jmt(m.get_jmt(k));
        fe_values.set_points(xq_hat);
        fe_values.prepare(jmt);
      
        const array<double>& psi(fe_values.template get_values<test_fe_index>());
        const array<double>& phi(fe_zvalues.template get_values<trial_fe_index>());

        // evaluate the weak form
        const double volume(m.get_cell_volume(k));
        for (unsigned int q(0); q < n_q; ++q) {
          integration_proxy.f.prepare(k, &xq.at(q, 0), &xq_hat.at(q, 0));
          for (unsigned int i(0); i < n_test_dof; ++i) {
            a_el.at(i) += omega.at(q) * integration_proxy.f(k,
                                                            &xq.at(q, 0), &xq_hat.at(q, 0),
                                                            &psi.at(q, i, 0),
                                                            &phi.at(q, 0, 0));
          }
        }
        for (unsigned int i(0); i < n_test_dof; ++i)
          b_form.accumulate(b_form.test_fes.get_dof(integration_proxy.get_global_cell_id(k), i),
                            b_form.trial_fes.get_dof_number() + id,
                            a_el.at(i) * volume);
      }
    }
  }

private:
  bilinear_form_type& b_form;
  std::size_t id;
};

template<typename test_fes_type, typename trial_fes_type>
typename bilinear_form<test_fes_type, trial_fes_type>::template algebraic_block_handle<algebraic_block::test_block>
bilinear_form<test_fes_type, trial_fes_type>::algebraic_test_block(std::size_t id) {
  return algebraic_block_handle<algebraic_block::test_block>(*this, id);
}

template<typename test_fes_type, typename trial_fes_type>
typename bilinear_form<test_fes_type, trial_fes_type>::template algebraic_block_handle<algebraic_block::trial_block>
bilinear_form<test_fes_type, trial_fes_type>::algebraic_trial_block(std::size_t id) {
  return algebraic_block_handle<algebraic_block::trial_block>(*this, id);
}

template<typename test_fes_type, typename trial_fes_type>
double& bilinear_form<test_fes_type, trial_fes_type>::algebraic_block(std::size_t id_1, std::size_t id_2) {
  return a.
    get(test_fes.get_dof_number() + id_1,
	     trial_fes.get_dof_number() + id_2);
}


#endif /* _BILINEAR_FORM_H_ */
