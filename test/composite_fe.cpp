
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

  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-value"

  composite_finite_element_space(const mesh<cell_type>& m)
    : fe_instances((sizeof(fe_pack), m)...) {}
  
  #pragma clang diagnostic pop


  std::size_t get_total_dof_number() const {
    return dof_number_sum_impl<cfe_type, 0, cfe_type::n_component - 1>::call(*this);
  }
  
  template<std::size_t n>
  std::size_t get_dof_number() const {
    return std::get<n>(fe_instances).get_dof_number();
  }

  template<std::size_t n>
  unsigned int get_dof(std::size_t k, std::size_t i) const {
    return std::get<n>(fe_instances).get_dof();
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

private:
  std::tuple<finite_element_space<fe_pack>...> fe_instances;
};


template<typename ... fe_pack>
struct composite_finite_element_space<composite_finite_element<fe_pack...> >::element {
  using cfe_type = composite_finite_element<fe_pack...>;
  using cfes_type = composite_finite_element_space<cfe_type>;
  using cell_type = typename cfes_type::cell_type;
  
  element(const cfes_type& cfes)
    : cfes(cfes), coefficients{cfes.get_total_dof_number()} {}

  element(const cfes_type& cfes,
	  const array<double>& a)
    : cfes(cfes), coefficients(a) {}

  element(const cfes_type& cfes,
	  array<double>&& a)
    : cfes(cfes), coefficients(a) {}
  
  element(const element& e)
    : cfes(e.cfes), coefficients(e.coefficients) {}

  ~element() {}

  const cfes_type& get_finite_element_space() const { return cfes; }
  const mesh<cell_type>& get_mesh() const { return cfes.get_mesh(); }
  const array<double>& get_coefficients() const { return coefficients; }

  template<std::size_t n>
  const array<double>& get_component() const;

private:
  const cfes_type& cfes;
  array<double> coefficients;
};


template<typename te_cfe_type, typename tr_cfe_type>
class bilinear_form<composite_finite_element_space<te_cfe_type>,
		    composite_finite_element_space<tr_cfe_type> > {
public:
  using test_cfes_type = composite_finite_element_space<te_cfe_type>;
  using trial_cfes_type = composite_finite_element_space<tr_cfe_type>;

  using test_cfe_type = typename test_cfes_type::cfe_type;
  using trial_cfe_type = typename trial_cfes_type::cfe_type;

  // both cfe should have the same number of components
  static const std::size_t n_component = test_cfe_type::n_component;


  template<typename cfes_type>
  struct mp_get_dof_number {
    template<std::size_t n>
    static std::size_t call(const cfes_type& cfes) {
      return cfes.template get_dof_number<n>();
    }
  };
  
  bilinear_form(const test_cfes_type& te_cfes,
		const trial_cfes_type& tr_cfes)
    : test_cfes(te_cfes), trial_cfes(tr_cfes),
      a(te_cfes.get_total_dof_number(),
	tr_cfes.get_total_dof_number()) {

    std::size_t test_global_dof_number[n_component];
    fill_array_with_return_values<std::size_t,
				  mp_get_dof_number<test_cfes_type>,
				  0,
				  n_component>::template fill<const test_cfes_type&>(&test_global_dof_number[0],
										     test_cfes);
    std::size_t trial_global_dof_number[n_component];
    fill_array_with_return_values<std::size_t,
				  mp_get_dof_number<trial_cfes_type>,
				  0,
				  n_component>::template fill<const trial_cfes_type&>(&trial_global_dof_number[0],
										      trial_cfes);

    test_global_dof_offset = std::vector<std::size_t>(n_component, 0ul);
    std::partial_sum(test_global_dof_number,
		     test_global_dof_number + n_component - 1,
		     test_global_dof_offset.begin() + 1);
    
    trial_global_dof_offset = std::vector<std::size_t>(n_component, 0ul);
    std::partial_sum(trial_global_dof_number,
		     trial_global_dof_number + n_component - 1,
		     trial_global_dof_offset.begin() + 1);
  }
  
  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 2, "linear_form expects rank-2 expression.");

    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    // m is the mesh over which we integrate
    const auto& m(integration_proxy.m);

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the point-wise basis function evaluation
    using fe_list = cat_list_t<typename test_cfe_type::fe_list,
			       typename test_cfe_type::fe_list>;
    using unique_fe_list = unique_t<fe_list>;
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
      
      /*
      //for m = 1:n_components
      // for n = 1:n_components
      double psi_phi[2 * n_component];
      fill_argument::fill<m, test_fe_list>(&psi[0], fe_values, fe_zvalues);
      fill_argument::fill<m, trial_fe_list>(&psi[n_component], fe_values, fe_zvalues);
      evaluate_weak_block<m, n>(integration_factor, psi_phi);
      //  endfor
      //endfor
      */
      
    }
  }

  template<typename T, std::size_t m, std::size_t n>
  void evaluate_block(const T& integration_proxy,
		      std::size_t k,
		      const array<double>& omega,
		      const array<double>& xq,
		      const array<double>& xq_hat,
		      const double** psi_phi) {

    typedef typename T::quadrature_type quadrature_type;
    
    using test_fe_type = get_element_at_t<m, typename test_cfe_type::fe_list>;
    using trial_fe_type = get_element_at_t<n, typename trial_cfe_type::fe_list>;
    
    const std::size_t n_q(quadrature_type::n_point);
    const std::size_t n_test_dof(test_fe_type::n_dof_per_element);
    const std::size_t n_trial_dof(trial_fe_type::n_dof_per_element);

    // evaluate the weak form
    const double volume(integration_proxy.m.get_cell_volume(k));
    for (unsigned int i(0); i < n_test_dof; ++i) {
	for (unsigned int j(0); j < n_trial_dof; ++j) {
	  double a_el(0.0);
	  for (unsigned int q(0); q < n_q; ++q) {
	    a_el += volume * omega.at(q)
	      * expression_call_wrapper<0, 2 * n_component>::call(integration_proxy.f, k,
								  &xq.at(q, 0), &xq_hat.at(q, 0),
								  psi_phi);
	  }
	  accumulate_in_block<m, n>(test_cfes.template get_dof<m>(integration_proxy.get_global_element_id(k), j),
				    trial_cfes.template get_dof<n>(integration_proxy.get_global_element_id(k), i), a_el);
	}
      }

  }

  template<std::size_t n>
  expression<form<n, 0, 0> > get_test_function() const {
    return form<n, 0, 0>();
  }

  template<std::size_t n>
  expression<form<n + n_component, 0, 0> > get_trial_function() const {
    return form<n + n_component, 0, 0>();
  }

  typename trial_cfes_type::element solve(const linear_form<test_cfes_type>& form) const;

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


void test_cfe() {
  using cell_type = cell::triangle;
  using fe_0_type = finite_element::triangle_lagrange_p0;
  using fe_1_type = finite_element::triangle_lagrange_p1;
  using fe_type = composite_finite_element<fe_0_type, fe_1_type>;
  using fes_type = composite_finite_element_space<fe_type>;
  
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 3, 3));
  fes_type fes(m);

  //fes.show<0>(std::cout);
  //fes.show<1>(std::cout);

  bilinear_form<fes_type, fes_type> a(fes, fes);
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

  test_meta();
  test_cfe();
  test_fe_value_manager();
  test_cfes();
  
  return 0;
}
