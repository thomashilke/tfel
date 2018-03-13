#ifndef _BILINEAR_FORM_H_
#define _BILINEAR_FORM_H_

template<typename test_fes_type, typename trial_fes_type>
class bilinear_form {
public:
  class constraint_handle;
  
  bilinear_form(const test_fes_type& te_fes,
		const trial_fes_type& tr_fes)
    : test_fes(te_fes), trial_fes(tr_fes),
      a(te_fes.get_dof_number(), tr_fes.get_dof_number()),
      petsc(
	linear_solver().get_solver(
	  solver::petsc,
	  method::gmres,
	  preconditioner::ilu)), constraint_number(0) {
    clear();
  }

  ~bilinear_form() {
    delete petsc; petsc = nullptr;
  }

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
	const array<double> jmt(m.get_jmt(k));
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

    dirty = true;
  }

  expression<form<0,1,0> > get_test_function() const { return form<0,1,0>(); }
  expression<form<1,2,0> > get_trial_function() const { return form<1,2,0>(); }

  constraint_handle new_constraint();
  void assemble_constraint(const constraint_handle& l_1, const constraint_handle l_2, double value);
  
  void prepare_solver() {
    petsc->set_size(test_fes.get_dof_number() + constraint_number);
    
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
      
    petsc->preallocate(&nnz[0]);

    if(true) {
      // Assemble line by line
      for (std::size_t row_id(0); row_id < row.size() - 1; ++row_id) {
	if (row[row_id + 1] > row[row_id])
	  petsc->add_row(row_id,
			 nnz[row_id],
			 &col[row[row_id]],
			 &val[row[row_id]]);
      }
    } else {
      // Assemble element by element
      for (const auto& v: a.get_elements())
	petsc->add_value(v.first.first, v.first.second, v.second);
    }
    petsc->assemble();
    //petsc->show();
    
    
    dirty = false;
  }
  
  typename trial_fes_type::element solve(const linear_form<test_fes_type>& form) {
    // Add the value of the dirichlet dof in the right hand side
    array<double> f{trial_fes.get_dof_number() + constraint_number};
    std::copy(&form.get_coefficients().at(0),
	      &form.get_coefficients().at(0) + trial_fes.get_dof_number(),
	      &f.at(0));
    std::copy(form.get_constraint_values().begin(),
	      form.get_constraint_values().end(),
	      &f.at(0) + trial_fes.get_dof_number());
    
    for (const auto& i: test_fes.get_dirichlet_dof()) {
      const auto x(trial_fes.get_dof_space_coordinate(i));
      f.at(i) = trial_fes.boundary_value(&x.at(0, 0));
    }
    /*std::cout << "f = [";
    for (std::size_t k(0); k < f.get_size(0); ++k) {
      std::cout << f.at(k) << std::endl;
    }
    std::cout << "];\n";*/
    
    if (dirty)
      prepare_solver();
    
    if (constraint_number == 0) {
      return typename trial_fes_type::element(trial_fes, petsc->solve(f));
    } else {
      array<double> solution(petsc->solve(f));
      array<double> coefficients{trial_fes.get_dof_number()};
      std::copy(&solution.at(0),
		&solution.at(0) + trial_fes.get_dof_number(),
		&coefficients.at(0));
      
      return typename trial_fes_type::element(trial_fes, coefficients);
    }
  }

  void clear() {
    a.clear();

    // Add the identity equations for each dirichlet dof
    for (const auto& i: test_fes.get_dirichlet_dof())
      a.accumulate(i, i, 1.0);

    constraint_number = 0;
  }
  
  void show(std::ostream& stream) const {
    a.show(stream);
  }

  void export_data(std::ostream& stream) const {
    a.export_data(stream);
  }
  
private:
  const test_fes_type& test_fes;
  const trial_fes_type& trial_fes;
  
  sparse_linear_system a;
  linear_solver_impl::solver_base* petsc;
  bool dirty = true;
  std::size_t constraint_number;

  void accumulate(std::size_t i, std::size_t j, double value) {
    if (test_fes.get_dirichlet_dof().count(i) == 0)
      a.accumulate(i, j, value);
  }
};

template<typename test_fes_type, typename trial_fes_type>
class bilinear_form<test_fes_type, trial_fes_type>::constraint_handle {
  using bilinear_form_type = bilinear_form<test_fes_type, trial_fes_type>;
  using test_fe_type = typename test_fes_type::fe_type;
  using trial_fe_type = typename trial_fes_type::fe_type;

  using fe_list = type_list<test_fe_type, trial_fe_type>;
  using unique_fe_list = unique_t<fe_list>;

  static const std::size_t test_fe_index = get_index_of_element<test_fe_type, unique_fe_list>::value;
  static const std::size_t trial_fe_index = get_index_of_element<trial_fe_type, unique_fe_list>::value;

public:
  constraint_handle(bilinear_form_type& f, std::size_t constraint_id)
    : b_form(f), constraint_id(constraint_id) {}

  std::size_t get_id() const { return constraint_id; }
  
  template<typename T>
  void operator-=(const T& integration_proxy) {
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
	b_form.accumulate(b_form.test_fes.get_dof_number() + constraint_id,
			  b_form.trial_fes.get_dof(integration_proxy.get_global_cell_id(k), j),
			  a_el.at(j) * volume);
    }

    b_form.dirty = true;
  }

  template<typename T>
  void operator|=(const T& integration_proxy) {
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
			  b_form.trial_fes.get_dof_number() + constraint_id,
			  a_el.at(i) * volume);
    }

    b_form.dirty = true;
  }

private:
  bilinear_form_type& b_form;
  std::size_t constraint_id;
};


template<typename test_fes_type, typename trial_fes_type>
typename bilinear_form<test_fes_type, trial_fes_type>::constraint_handle
bilinear_form<test_fes_type, trial_fes_type>::new_constraint() {
  constraint_number += 1;
  a.set_sizes(test_fes.get_dof_number() + constraint_number,
	      trial_fes.get_dof_number() + constraint_number);
  return constraint_handle(*this, constraint_number - 1);
}


template<typename test_fes_type, typename trial_fes_type>
void bilinear_form<test_fes_type, trial_fes_type>::assemble_constraint(const constraint_handle& l_1, const constraint_handle l_2, double value) {
  accumulate(test_fes.get_dof_number() + l_1.get_id(),
	     trial_fes.get_dof_number() + l_2.get_id(), value);
}


#endif /* _BILINEAR_FORM_H_ */
