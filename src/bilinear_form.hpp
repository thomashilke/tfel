#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

template<typename test_fes_type, typename trial_fes_type>
class bilinear_form {
public:
  bilinear_form(const test_fes_type& te_fes,
		const trial_fes_type& tr_fes)
    : test_fes(te_fes), trial_fes(tr_fes),
      a(te_fes.get_dof_number(), tr_fes.get_dof_number()) {}

  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 2, "linear_form expects rank-2 expression.");
      
    typedef typename test_fes_type::fe_type test_fe_type;
    typedef typename trial_fes_type::fe_type trial_fe_type;
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    const auto& m(integration_proxy.m);

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the point-wise basis function evaluation
    using fe_list = type_list<test_fe_type, trial_fe_type>;
    using unique_fe_list = unique_t<fe_list>;

    const std::size_t test_fe_index(get_index_of_element<test_fe_type, unique_fe_list>::value);
    const std::size_t trial_fe_index(get_index_of_element<trial_fe_type, unique_fe_list>::value);
    fe_value_manager<unique_fe_list> fe_values((n_q));


    // loop over the elements
    for (unsigned int k(0); k < m.get_element_number(); ++k) {
      // prepare the quadrature points
      const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
      const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
									m.get_elements(),
									k, xq_hat));
      // prepare the basis function values
      const array<double> jmt(m.get_jmt(k));
      fe_values.prepare(jmt, xq_hat);
      
      const array<double>& psi(fe_values.template get_values<test_fe_index>());
      const array<double>& phi(fe_values.template get_values<trial_fe_index>());

      // evaluate the weak form
      const std::size_t n_test_dof(test_fe_type::n_dof_per_element);
      const std::size_t n_trial_dof(trial_fe_type::n_dof_per_element);
      const double volume(m.get_cell_volume(k));
      for (unsigned int i(0); i < n_test_dof; ++i) {
	for (unsigned int j(0); j < n_trial_dof; ++j) {
	  double a_el(0.0);
	  for (unsigned int q(0); q < n_q; ++q) {
	    a_el += volume * omega.at(q) * integration_proxy.f(k,
							       &xq.at(q, 0), &xq_hat.at(q, 0),
							       &psi.at(q, i, 0),
							       &phi.at(q, j, 0));
	  }
	  accumulate(test_fes.get_dof(integration_proxy.get_global_element_id(k), j),
		     trial_fes.get_dof(integration_proxy.get_global_element_id(k), i), a_el);
	}
      }
    }
  }

  expression<form<0,1,0> > get_test_function() const { return form<0,1,0>(); }
  expression<form<1,2,0> > get_trial_function() const { return form<1,2,0>(); }

  typename trial_fes_type::element solve(const linear_form<test_fes_type>& form) {
    // Add the identity equations for each dirichlet dof
    for (const auto& i: test_fes.get_dirichlet_dof())
      a.accumulate(i, i, 1.0);

    // Add the value of the dirichlet dof in the right hand side
    array<double> f(form.get_coefficients());
    for (const auto& i: test_fes.get_dirichlet_dof()) {
      const auto x(trial_fes.get_dof_space_coordinate(i));
      f.at(i) = trial_fes.boundary_value(&x.at(0, 0));
    }
    
    linear_solver s;
    auto petsc_gmres_ilu(s.get_solver(solver::petsc,
				      method::gmres,
				      preconditioner::ilu,
				      test_fes.get_dof_number()));
    if (true) {
      // Convert to CRS format
      std::vector<int>
	row(a.get_equation_number() + 1),
	col(a.get_elements().size());
      std::vector<double>
	val(a.get_elements().size());
    
      std::size_t row_id(0), val_id(0);
      row[0] = row_id;
      for (const auto& v: a.get_elements()) {
	while (row_id < v.first.first) {
	  ++row_id;
	  row[row_id] = val_id;
	}
	col[val_id] = v.first.second;
	val[val_id] = v.second;
	++val_id;
      }
      row.back() = val_id;

      // count the number of non-zero per row
      std::vector<int> nnz(a.get_equation_number());
      for (std::size_t n(0); n < a.get_equation_number(); ++n)
	nnz[n] = row[n + 1] - row[n];
      
      petsc_gmres_ilu->preallocate(&nnz[0]);

      if(true) {
	// Assemble line by line
	for (std::size_t row_id(0); row_id < row.size() - 1; ++row_id) {
	  if (row[row_id + 1] > row[row_id])
	    petsc_gmres_ilu->add_row(row_id,
				     nnz[row_id],
				     &col[row[row_id]],
				     &val[row[row_id]]);
	}
      } else {
	// Assemble element by element
	for (const auto& v: a.get_elements())
	  petsc_gmres_ilu->add_value(v.first.first, v.first.second, v.second);
      }
    }

    petsc_gmres_ilu->assemble();

    typename trial_fes_type::element result(trial_fes, petsc_gmres_ilu->solve(f));
    delete petsc_gmres_ilu; petsc_gmres_ilu = nullptr;
    return result;
  }

  void clear() {
    a.clear();
  }
  
  void show(std::ostream& stream) {
    a.show(stream);
  }
  
private:
  const test_fes_type& test_fes;
  const trial_fes_type& trial_fes;
  
  sparse_linear_system a;

  void accumulate(std::size_t i, std::size_t j, double value) {
    if (test_fes.get_dirichlet_dof().count(i) == 0)
      a.accumulate(i, j, value);
  }
};

#endif /* _BILINEAR_FORM_H_ */
