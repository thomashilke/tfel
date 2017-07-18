
#include <type_traits>
#include <cstddef>
#include <iostream>
#include <typeinfo>


#include "../src/mesh.hpp"
#include "../src/fe.hpp"
#include "../src/fes.hpp"
#include "../src/meta.hpp"
#include "../src/fe_value_manager.hpp"
#include "../src/form.hpp"
#include "../src/quadrature.hpp"
#include "../src/export.hpp"
#include "../src/fes.hpp"


/*
 *  A composite finite element is simply a wrapper around a typelist
 *  with an interface to access each of the finite element types in the list.
 */
template<typename ... fe_pack>
struct composite_finite_element {
  using fe_list = type_list<fe_pack...>;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;

  static const std::size_t n_component = list_size<fe_list>::value;
};


/*
 *  Declaration of a composite finite element space.
 *  A composite finite element space is a wrapper around an std::tuple
 *  where each component is a simple finite element space.
 */
template<typename cfe_type>
class composite_finite_element_space;


template<typename cfe_type, std::size_t n, std::size_t n_max>
struct dof_number_sum_impl {
  static std::size_t call(const composite_finite_element_space<cfe_type>& cfes) {
    return cfes.template get_dof_number<n>() + dof_number_sum_impl<cfe_type, n + 1, n_max>::call(cfes);
  }
};

template<typename cfe_type, std::size_t n_max>
struct dof_number_sum_impl<cfe_type, n_max, n_max> {
  static std::size_t call(const composite_finite_element_space<cfe_type>& cfes) {
    return 0;
  }
};


template<typename R, typename F, std::size_t n, std::size_t n_max>
struct fill_array_with_return_values {
  template<typename ... As>
  static void fill(R* ptr, As... as) {
    *ptr = F::template call<n>(as...);
    fill_array_with_return_values<R, F, n + 1, n_max>::template fill<As...>(ptr + 1, as...);
  }
};

template<typename R, typename F, std::size_t n_max>
struct fill_array_with_return_values<R, F, n_max, n_max> {
  template<typename ... As>
  static void fill(R* ptr, As... as) {}
};


template<typename ... fe_pack>
class composite_finite_element_space<composite_finite_element<fe_pack...> > {
public:
  struct element;
  using cfe_type = composite_finite_element<fe_pack...>;
  using cell_list = unique_t<transform<get_cell_type, typename cfe_type::fe_list> >;
  using cell_type = get_element_at_t<0, cell_list>;
  using fe_list = type_list<fe_pack...>;

  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-value"

  composite_finite_element_space(const mesh<cell_type>& m)
    : fe_instances((sizeof(fe_pack), m)...) {}
  
  #pragma clang diagnostic pop

  template<std::size_t n>
  void set_dirichlet_boundary_condition(const submesh<cell_type>& dm) {
    std::get<n>(fe_instances).set_dirichlet_boundary_condition(dm);
  }
  
  std::size_t get_total_dof_number() const {
    return dof_number_sum_impl<cfe_type, 0, cfe_type::n_component>::call(*this);
  }
  
  template<std::size_t n>
  std::size_t get_dof_number() const {
    return std::get<n>(fe_instances).get_dof_number();
  }

  template<std::size_t n>
  unsigned int get_dof(std::size_t k, std::size_t i) const {
    return std::get<n>(fe_instances).get_dof(k, i);
  }

  template<std::size_t n>
  const std::unordered_set<unsigned int>& get_dirichlet_dof() const {
    return std::get<n>(fe_instances).get_dirichlet_dof();
  }
	
  template<std::size_t n>
  const std::vector<std::set<cell::subdomain_type> > get_subdomain_list() const {
    return std::get<n>(fe_instances).get_subdomain_list();
  }

  template<std::size_t n>
  void show(std::ostream& stream) {
    std::get<n>(fe_instances).show(stream);
  }
    
  const mesh<cell_type>& get_mesh() const {
    return std::get<0>(fe_instances).get_mesh();
  }

  template<std::size_t n>
  double boundary_value(const double* x) const {
    return std::get<n>(fe_instances).boundary_value(x);
  }

  template<std::size_t n>
  array<double> get_dof_space_coordinate(unsigned int i) const {
    return std::get<n>(fe_instances).get_dof_space_coordinate(i);
  }

  template<std::size_t n>
  const finite_element_space<get_element_at_t<n, fe_list> >& get_finite_element_space() const {
    return std::get<n>(fe_instances);
  }

private:
  std::tuple<finite_element_space<fe_pack>...> fe_instances;
};


template<typename cfes_type>
struct get_dof_number_impl {
  template<std::size_t n>
  static std::size_t call(const cfes_type& cfes) {
    return cfes.template get_dof_number<n>();
  }
};


template<typename ... fe_pack>
struct composite_finite_element_space<composite_finite_element<fe_pack...> >::element {
  using cfe_type = composite_finite_element<fe_pack...>;
  using cfes_type = composite_finite_element_space<cfe_type>;
  using cell_type = typename cfes_type::cell_type;
  using fe_list = type_list<fe_pack...>;
  
  static const std::size_t n_component = sizeof...(fe_pack);

  element(const cfes_type& cfes)
    : cfes(cfes), coefficients{cfes.get_total_dof_number()} {
    setup_offsets();
  }

  element(const cfes_type& cfes,
	  const array<double>& a)
    : cfes(cfes), coefficients(a) {
      setup_offsets();
  }

  element(const cfes_type& cfes,
	  array<double>&& a)
    : cfes(cfes), coefficients(a) {
      setup_offsets();
  }
  
  element(const element& e)
    : cfes(e.cfes), coefficients(e.coefficients) {
      setup_offsets();
  }

  ~element() {}

  const cfes_type& get_finite_element_space() const { return cfes; }
  const mesh<cell_type>& get_mesh() const { return cfes.get_mesh(); }
  const array<double>& get_coefficients() const { return coefficients; }

  template<std::size_t n>
  typename finite_element_space<get_element_at_t<n, fe_list> >::element get_component() const {
    array<double> cf{dof_numbers[n]};
    std::copy(&coefficients.at(dof_offsets[n]),
	      &coefficients.at(dof_offsets[n]) + dof_numbers[n],
	      &cf.at(0));
    return typename finite_element_space<get_element_at_t<n, fe_list> >::element(cfes.template get_finite_element_space<n>(),
										 cf);
  }

private:
  const cfes_type& cfes;
  array<double> coefficients;
  std::vector<std::size_t> dof_numbers, dof_offsets;
  
private:
  void setup_offsets() {
    dof_numbers = std::vector<std::size_t>(n_component, 0);
    fill_array_with_return_values<std::size_t,
				  get_dof_number_impl<cfes_type>,
				  0,
				  n_component>::template fill<const cfes_type&>(&dof_numbers[0], cfes);
    
    dof_offsets = std::vector<std::size_t>(n_component, 0ul);
    std::partial_sum(dof_numbers.begin(),
		     dof_numbers.begin() + n_component - 1,
		     dof_offsets.begin() + 1);
  }
};

template<typename fe_list, std::size_t n, std::size_t m, typename unique_fe_list>
struct select_function_valuation_impl {
  static void call(const double** ptr, std::size_t q, std::size_t i,
		   const fe_value_manager<unique_fe_list>& values,
		   const fe_value_manager<unique_fe_list>& zvalues) {
    const std::size_t fe_index = get_index_of_element<head_t<fe_list>, unique_fe_list>::value;
    
    *ptr = &(zvalues.template get_values<fe_index>().at(q,0,0));
    select_function_valuation_impl<tail_t<fe_list>, n + 1, m, unique_fe_list>::call(ptr + 1, q, i, values, zvalues);
  }
};

template<typename fe_list, std::size_t m, typename unique_fe_list>
struct select_function_valuation_impl<fe_list, m, m, unique_fe_list> {
  static void call(const double** ptr, std::size_t q, std::size_t i,
		   const fe_value_manager<unique_fe_list>& values,
		   const fe_value_manager<unique_fe_list>& zvalues) {
    const std::size_t fe_index = get_index_of_element<head_t<fe_list>, unique_fe_list>::value;
    
    *ptr = &(values.template get_values<fe_index>().at(q,i,0));
    select_function_valuation_impl<tail_t<fe_list>, m + 1, m, unique_fe_list>::call(ptr + 1, q, i, values, zvalues);
  }
};

template<std::size_t n, std::size_t m, typename unique_fe_list>
struct select_function_valuation_impl<type_list<>, n, m, unique_fe_list> {
  static void call(const double** ptr, std::size_t q, std::size_t i,
		   const fe_value_manager<unique_fe_list>& values,
		   const fe_value_manager<unique_fe_list>& zvalues) {}
};



template<typename fe_list, std::size_t m, typename unique_fe_list>
void select_function_valuation(const double** ptr, std::size_t q, std::size_t i,
			       const fe_value_manager<unique_fe_list>& values,
			       const fe_value_manager<unique_fe_list>& zvalues) {
  select_function_valuation_impl<fe_list, 0, m, unique_fe_list>::call(ptr, q, i, values, zvalues);
}



template<typename te_cfe_type>
class linear_form<composite_finite_element_space<te_cfe_type> > {
public:
  using test_cfes_type = composite_finite_element_space<te_cfe_type>;
  using test_cfe_type = typename test_cfes_type::cfe_type;
  using test_fe_list = typename test_cfe_type::fe_list;
  using fe_list = test_fe_list;
  using unique_fe_list = unique_t<fe_list>;

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

    // m is the mesh over which we integrate
    const auto& m(integration_proxy.m);

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the point-wise basis function evaluation
    fe_value_manager<unique_fe_list> fe_values(n_q), fe_zvalues(n_q);
    fe_zvalues.clear();
    
    for (unsigned int k(0); k < m.get_element_number(); ++k) {
      // prepare the quadrature points
      const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
      const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
									m.get_elements(),
									k, xq_hat));
      // prepare the basis function values
      const array<double> jmt(m.get_jmt(k));
      fe_values.prepare(jmt, xq_hat);


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


template<typename te_cfe_type, typename tr_cfe_type>
class bilinear_form<composite_finite_element_space<te_cfe_type>,
		    composite_finite_element_space<tr_cfe_type> > {
public:
  using bilinear_form_type = bilinear_form<composite_finite_element_space<te_cfe_type>,
					   composite_finite_element_space<tr_cfe_type> >;
  
  using test_cfes_type = composite_finite_element_space<te_cfe_type>;
  using trial_cfes_type = composite_finite_element_space<tr_cfe_type>;

  using test_cfe_type = typename test_cfes_type::cfe_type;
  using trial_cfe_type = typename trial_cfes_type::cfe_type;

  using test_fe_list = typename test_cfe_type::fe_list;
  using trial_fe_list = typename trial_cfe_type::fe_list;
  using fe_list = cat_list_t<typename test_cfe_type::fe_list,
			     typename test_cfe_type::fe_list>;
  using unique_fe_list = unique_t<fe_list>;

  // both cfe should have the same number of components
  static const std::size_t n_test_component = test_cfe_type::n_component;
  static const std::size_t n_trial_component = trial_cfe_type::n_component;


  template<typename IC>
  struct handle_dirichlet_dof_equations {
    static const std::size_t m = IC::value;
    static void call(const bilinear_form_type& bilinear_form, sparse_linear_system& a) {
      for (const auto& i: bilinear_form.test_cfes.template get_dirichlet_dof<m>()) {
	a.accumulate(bilinear_form.test_global_dof_offset[m] + i,
		     bilinear_form.test_global_dof_offset[m] + i,
		     1.0);
      }      
    }
  };
  
  bilinear_form(const test_cfes_type& te_cfes,
		const trial_cfes_type& tr_cfes)
    : test_cfes(te_cfes), trial_cfes(tr_cfes),
      a(te_cfes.get_total_dof_number(),
	tr_cfes.get_total_dof_number()) {
    std::size_t test_global_dof_number[n_test_component];
    fill_array_with_return_values<std::size_t,
				  get_dof_number_impl<test_cfes_type>,
				  0,
				  n_test_component>::template fill<const test_cfes_type&>(&test_global_dof_number[0],
											  test_cfes);
    std::size_t trial_global_dof_number[n_trial_component];
    fill_array_with_return_values<std::size_t,
				  get_dof_number_impl<trial_cfes_type>,
				  0,
				  n_trial_component>::template fill<const trial_cfes_type&>(&trial_global_dof_number[0],
											    trial_cfes);

    test_global_dof_offset = std::vector<std::size_t>(n_test_component, 0ul);
    std::partial_sum(test_global_dof_number,
		     test_global_dof_number + n_test_component - 1,
		     test_global_dof_offset.begin() + 1);
    
    trial_global_dof_offset = std::vector<std::size_t>(n_trial_component, 0ul);
    std::partial_sum(trial_global_dof_number,
		     trial_global_dof_number + n_trial_component - 1,
		     trial_global_dof_offset.begin() + 1);

    // we need to specify the equation for the dirichlet dof
    call_for_each<handle_dirichlet_dof_equations, make_integral_list_t<std::size_t, n_test_component> >::call(*this, a);
  }


  template<typename B_INFO>
  struct evaluate_block {
    using T = get_element_at_t<0, B_INFO>;
    typedef typename T::quadrature_type quadrature_type;

    static const std::size_t m = get_element_at_t<1, B_INFO>::value;
    static const std::size_t n = get_element_at_t<2, B_INFO>::value;

    using test_fe_type = get_element_at_t<m, typename test_cfe_type::fe_list>;
    using trial_fe_type = get_element_at_t<n, typename trial_cfe_type::fe_list>;

    static void call(bilinear_form_type& bilinear_form,
		     const get_element_at_t<0, B_INFO>& integration_proxy,
		     const std::size_t k,
		     const array<double>& omega,
		     const array<double>& xq,
		     const array<double>& xq_hat,
		     const fe_value_manager<unique_fe_list>& fe_values,
		     const fe_value_manager<unique_fe_list>& fe_zvalues) {
      const double* psi_phi[n_test_component + n_trial_component];

      const std::size_t n_q(quadrature_type::n_point);
      const std::size_t n_test_dof(test_fe_type::n_dof_per_element);
      const std::size_t n_trial_dof(trial_fe_type::n_dof_per_element);

      // evaluate the weak form
      const double volume(integration_proxy.m.get_cell_volume(k));
      for (unsigned int i(0); i < n_test_dof; ++i) {
	for (unsigned int j(0); j < n_trial_dof; ++j) {
	  double a_el(0.0);
	  for (unsigned int q(0); q < n_q; ++q) {
	    select_function_valuation<test_fe_list, m, unique_fe_list>(psi_phi, q, i,
								       fe_values, fe_zvalues);
	    select_function_valuation<trial_fe_list, n, unique_fe_list>(psi_phi + n_test_component, q, j,
									fe_values, fe_zvalues);

	    a_el += volume * omega.at(q)
	      * (expression_call_wrapper<0, n_test_component + n_trial_component>
		 ::call(integration_proxy.f, psi_phi,
			k, &xq.at(q, 0), &xq_hat.at(q, 0)));
	  }
	  bilinear_form.accumulate_in_block<m, n>(bilinear_form.test_cfes.
						  template get_dof<m>(integration_proxy.
								      get_global_element_id(k), i),
						  bilinear_form.trial_cfes.
						  template get_dof<n>(integration_proxy.
								      get_global_element_id(k), j), a_el);
	}
      }
    }
  };

  
  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 2, "bilinear_form expects rank-2 expression.");

    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    // m is the mesh over which we integrate
    const auto& m(integration_proxy.m);

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the point-wise basis function evaluation
    fe_value_manager<unique_fe_list> fe_values(n_q), fe_zvalues(n_q);
    fe_zvalues.clear();
    
    for (unsigned int k(0); k < m.get_element_number(); ++k) {
      // prepare the quadrature points
      const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
      const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
									m.get_elements(),
									k, xq_hat));
      // prepare the basis function values
      const array<double> jmt(m.get_jmt(k));
      fe_values.prepare(jmt, xq_hat);
      

      /*
       *  Compile time loop over all the blocks
       */
      using test_blocks_il = make_integral_list_t<std::size_t, n_test_component>;
      using trial_blocks_il = make_integral_list_t<std::size_t, n_trial_component>;
      using block_list = tensor_product_of_lists_t<test_blocks_il, trial_blocks_il>;
      using block_info = append_to_each_element_t<T, block_list>;

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

  template<std::size_t n>
  expression<form<n + n_test_component, 2, 0> > get_trial_function() const {
    return form<n + n_test_component, 2, 0>();
  }

  
  template<typename IC>
  struct handle_dirichlet_dof_values {
    static const std::size_t m = IC::value;
    static void call(const bilinear_form_type& bilinear_form, array<double>& f) {
      for (const auto& i: bilinear_form.test_cfes.template get_dirichlet_dof<m>()) {
	const auto x(bilinear_form.trial_cfes.template get_dof_space_coordinate<m>(i));
	f.at(bilinear_form.test_global_dof_offset[m] + i) = bilinear_form.trial_cfes.template boundary_value<m>(&x.at(0, 0));
      }      
    }
  };
    
  typename trial_cfes_type::element solve(const linear_form<test_cfes_type>& form) const {
    array<double> f(form.get_coefficients());
    call_for_each<handle_dirichlet_dof_values, make_integral_list_t<std::size_t, n_test_component> >::call(*this, f);
    
    linear_solver s;
    auto petsc_gmres_ilu(s.get_solver(solver::petsc,
				      method::gmres,
				      preconditioner::ilu,
				      test_cfes.get_total_dof_number()));
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
      for (std::size_t n(0); n < a.get_equation_number(); ++n) {
	nnz[n] = row[n + 1] - row[n];
      }
      
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
    //petsc_gmres_ilu->show();

    typename trial_cfes_type::element result(trial_cfes, petsc_gmres_ilu->solve(f));
    delete petsc_gmres_ilu; petsc_gmres_ilu = nullptr;
    return result;
  }

  
  void clear();

  void show(std::ostream& stream);
  
private:
  const test_cfes_type& test_cfes;
  const trial_cfes_type& trial_cfes;

  std::vector<std::size_t> test_global_dof_offset;
  std::vector<std::size_t> trial_global_dof_offset;

  sparse_linear_system a;

  template<std::size_t m, std::size_t n>
  void accumulate_in_block(std::size_t i, std::size_t j, double value) {
    if (test_cfes.template get_dirichlet_dof<m>().count(i) == 0)
      a.accumulate(i + test_global_dof_offset[m],
		   j + trial_global_dof_offset[n],
		   value);
  }
};


void test_meta() {
  {
    using tl = type_list<double, char, int>;
    using tlt = transform<identity, tl>;

    using tl_element_1 = get_element_at_t<0, tl>;
    using tl_element_2 = get_element_at_t<1, tl>;
    using tl_element_3 = get_element_at_t<2, tl>;
    //using tl_element_3 = get_t<3, tl>; // will not compile

    print_type<tl_element_1>("element 0: ");
    print_type<tl_element_2>("element 1: ");
    print_type<tl_element_3>("element 2: ");

    print_type<tlt>("transformed by the identity: ");

    using i_tl = make_index_list_t<tl>;
    print_type<i_tl>("index list generated from tl: ");
    
    std::cout << is_member<int, tl>::value << std::endl;
    std::cout << is_member<double, tl>::value << std::endl;
    std::cout << is_member<char, tl>::value << std::endl;
    std::cout << is_member<short, tl>::value << std::endl;

    std::cout << "index of double is: " << get_index_of_element<double, tl>::value << std::endl;
    std::cout << "index of char is: " << get_index_of_element<char, tl>::value << std::endl;
    std::cout << "index of int is: " << get_index_of_element<int, tl>::value << std::endl;
  }

  {
    using tl = type_list<double, char, double>;

    print_type<unique_t<tl> >();
  }
}


void test_meta_integral_sequence() {
  using is = make_integral_sequence_t<std::size_t, 2>;
  using il = sequence_to_list_t<is>;
  using ill = wrap_t<type_list, il>;
  using ic = integral_constant<int, 7>;
  using ic_ill = append_to_each_element_t<ic, ill>;
  using tp = tensor_product_of_lists_t<il, il>;
  
  print_type<il>("list of integral_constant: ");

  print_type<wrap_t<type_list, il> >("wrap il elements in type_list's: ");
  print_type<ic_ill>("ic appended: ");

  print_type<tp>("tensor_product of the list with itself: ");
  print_type<append_to_each_element_t<int, tp> >("tensor product with int appended: ");
}


double f_0(const double* x) { return 1.0; }
double f_1(const double* x) { return -1.0; }

void test_cfe() {
  using cell_type = cell::triangle;
  using fe_0_type = finite_element::triangle_lagrange_p1;
  using fe_1_type = finite_element::triangle_lagrange_p1;
  using fe_type = composite_finite_element<fe_0_type, fe_1_type>;
  using fes_type = composite_finite_element_space<fe_type>;
  
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 10, 10));
  submesh<cell_type> dm(m.get_boundary_submesh());
  fes_type fes(m);

  fes.set_dirichlet_boundary_condition<0>(dm);
  fes.set_dirichlet_boundary_condition<1>(dm);
  
  //fes.show<0>(std::cout);
  //fes.show<1>(std::cout);


  bilinear_form<fes_type, fes_type> a(fes, fes); {
    const auto u_0(a.get_trial_function<0>());
    const auto u_1(a.get_trial_function<1>());
    const auto v_0(a.get_test_function<0>());
    const auto v_1(a.get_test_function<1>());

    a += integrate<quad::triangle::qf5pT>(d<1>(u_0) * d<1>(v_0) + d<2>(u_0) * d<2>(v_0)
					  + d<1>(u_1) * d<1>(v_1) + d<2>(u_1) * d<2>(v_1), m);
  }

  linear_form<fes_type> f(fes); {
    const auto v_0(f.get_test_function<0>());
    const auto v_1(f.get_test_function<1>());

    f += integrate<quad::triangle::qf5pT>(f_0 * v_0 + f_1 * v_1, m);
  }
  //f.show(std::cout);

  const fes_type::element x(a.solve(f));

  exporter::ensight6<fe_0_type>("laplacien_0", x.get_component<0>(), "solution");
  exporter::ensight6<fe_1_type>("laplacien_1", x.get_component<1>(), "solution");
}


void test_fe_value_manager() {
  using fe_list = type_list<finite_element::triangle_lagrange_p0,
			    finite_element::triangle_lagrange_p1>;

  const std::size_t n_q(3);
  array<double> xq{3, 2};
  xq.at(0, 0) = 0.0;
  xq.at(0, 0) = 0.0;
    
  xq.at(0, 0) = 1.0;
  xq.at(0, 0) = 0.0;
	
  xq.at(0, 0) = 0.0;
  xq.at(0, 0) = 1.0;

  array<double> jmt{2, 2};
  jmt.at(0, 0) = 1.0;
  jmt.at(0, 1) = 0.0;
  jmt.at(1, 0) = 0.0;
  jmt.at(1, 1) = 1.0;
  
  fe_value_manager<fe_list> fe_values(n_q);
  fe_values.prepare(jmt, xq);

  std::cout << fe_values.get_values<0ul>().at(0,0,0) << std::endl;
  const array<double>& phi(fe_values.get_values<0>());
  const array<double>& psi(fe_values.get_values<1>());
  
  std::cout << phi.at(0,0,0) << std::endl;
  std::cout << psi.at(0,0,0) << std::endl;
}


void test_cfes() {
  
}


int main(int argc, char *argv[]) {
  try {
  //test_meta();
  test_cfe();
  //test_fe_value_manager();
  //test_cfes();
  //test_meta_integral_sequence();
  } catch (std::string& e) {
    std::cout << e << std::endl;
  }

  return 0;
}
