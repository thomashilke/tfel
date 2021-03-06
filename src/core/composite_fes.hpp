#ifndef _COMPOSITE_FES_H_
#define _COMPOSITE_FES_H_

#include "fes.hpp"
#include "composite_fe.hpp"


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

template<typename cfes_type>
struct get_dof_number_impl {
  template<std::size_t n>
  static std::size_t call(const cfes_type& cfes) {
    return cfes.template get_dof_number<n>();
  }
};


template<typename ... fe_pack>
class composite_finite_element_space<composite_finite_element<fe_pack...> > {
public:
  struct element;
  using fe_list = type_list<fe_pack...>;
  using cfe_type = composite_finite_element<fe_pack...>;
  using cell_list = unique_t<transform<get_cell_type, fe_list> >;
  using cell_type = get_element_at_t<0, cell_list>;
  
  static_assert(list_size<cell_list>::value == 1, "");
  
  template<std::size_t n>
  using fes_type = finite_element_space<get_element_at_t<n, fe_list> >;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;

  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-value"

  composite_finite_element_space(const fe_mesh<cell_type>& m)
    : fe_instances((sizeof(fe_pack), m)...) {}
  
  #pragma clang diagnostic pop

  template<std::size_t n, typename c_cell_type>
  void add_dirichlet_boundary(const submesh<cell_type, c_cell_type>& dm,
                              double value = 0.0) {
    std::get<n>(fe_instances).add_dirichlet_boundary(dm, value);
  }

  template<std::size_t n, typename c_cell_type>
  void add_dirichlet_boundary(const submesh<cell_type, c_cell_type>& dm,
                              const std::function<double(const double*)>& f_bc) {
    std::get<n>(fe_instances).add_dirichlet_boundary(dm, f_bc);
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
  const std::map<unsigned int, double>& get_dirichlet_dof_values() const {
    return std::get<n>(fe_instances).get_dirichlet_dof_values();
  }
	
  template<std::size_t n>
  const std::vector<std::set<cell::subdomain_type> > get_subdomain_list() const {
    return std::get<n>(fe_instances).get_subdomain_list();
  }

  template<std::size_t n>
  void show(std::ostream& stream) {
    std::get<n>(fe_instances).show(stream);
  }
    
  const fe_mesh<cell_type>& get_mesh() const {
    return std::get<0>(fe_instances).get_mesh();
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


template<std::size_t component, std::size_t n_component>
struct setup_coefficients {
  template<typename fes_element_type, typename ... Es>
  static void
  call(array<double>& coefficients,
       const std::vector<std::size_t>& dof_numbers,
       const std::vector<std::size_t>& dof_offsets,
       fes_element_type& element,
       Es&& ... es) {
    std::copy(&element.get_coefficients().at(0),
	      &element.get_coefficients().at(0) + dof_numbers[component],
	      &coefficients.at(dof_offsets[component]));
    setup_coefficients<component + 1, n_component>::call(coefficients, dof_numbers, dof_offsets,
					    std::forward<Es>(es)...);
  }
};

template<std::size_t n_component>
struct setup_coefficients<n_component, n_component> {
  static void
  call(array<double>& coefficients,
       const std::vector<std::size_t>& dof_numbers,
       const std::vector<std::size_t>& dof_offsets) {}
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


  element(const cfes_type& cfes,
	  const typename finite_element_space<fe_pack>::element& ... elements)
    : cfes(cfes), coefficients{cfes.get_total_dof_number()} {
    setup_offsets();
    setup_coefficients<0, n_component>::call(coefficients, dof_numbers, dof_offsets, elements...);
  }

  ~element() {}

  const cfes_type& get_finite_element_space() const { return cfes; }
  const fe_mesh<cell_type>& get_mesh() const { return cfes.get_mesh(); }
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

  typename composite_finite_element_space<composite_finite_element<fe_pack...> >::element&
  operator=(const typename composite_finite_element_space<composite_finite_element<fe_pack...> >::element& e) {
    if (&cfes != &(e.cfes))
      throw std::string("Assigment of elements between different finite element spaces is not supported.");
    coefficients = e.coefficients;
    dof_numbers = e.dof_numbers;
    dof_offsets = e.dof_offsets;

    return *this;
  }

  template<std::size_t n, std::size_t n_max, typename cfes_type>
  struct restrict_impl {
    using submesh_type = submesh<cell_type, cell_type>;
    
    template<typename ... Cs>
    static typename cfes_type::element
    call(const submesh_type& sm, const cfes_type& sm_cfes, const typename cfes_type::element& element,
	 Cs&& ... cs) {
      return restrict_impl<n + 1, n_max, cfes_type>::call(
        sm, sm_cfes, element, std::forward<Cs>(cs)...,
	element.template get_component<n>().restrict(sm_cfes.template get_finite_element_space<n>(), sm));
    }
  };

  template<std::size_t n, typename cfes_type>
  struct restrict_impl<n, n, cfes_type> {
    using submesh_type = submesh<cell_type, cell_type>;
    
    template<typename ... Cs>
    static typename cfes_type::element
    call(const submesh_type& sm, const cfes_type& sm_cfes, const typename cfes_type::element& element,
	 Cs&& ... cs) {
      return typename cfes_type::element(sm_cfes, std::forward<Cs>(cs)...);
    }
  };
  
  typename composite_finite_element_space<composite_finite_element<fe_pack...> >::element
  restrict(const composite_finite_element_space<composite_finite_element<fe_pack...> >& sm_cfes,
	   const submesh<cell_type, cell_type>& sm) const {
    return restrict_impl<0, n_component, cfes_type>::call(sm, sm_cfes, *this);
  }

  
  template<std::size_t n, std::size_t n_max, typename cfes_type>
  struct extend_impl {
    using submesh_type = submesh<cell_type, cell_type>;
    
    template<typename ... Cs>
    static typename cfes_type::element
    call(const submesh_type& sm, const cfes_type& m_cfes, const typename cfes_type::element& element,
	 Cs&& ... cs) {
      return extend_impl<n + 1, n_max, cfes_type>::call(
        sm, m_cfes, element, std::forward<Cs>(cs)...,
	element.template get_component<n>().extend(m_cfes.template get_finite_element_space<n>(), sm));
    }
  };

  template<std::size_t n, typename cfes_type>
  struct extend_impl<n, n, cfes_type> {
    using submesh_type = submesh<cell_type, cell_type>;
    
    template<typename ... Cs>
    static typename cfes_type::element
    call(const submesh_type& sm, const cfes_type& m_cfes, const typename cfes_type::element& element,
	 Cs&& ... cs) {
      return typename cfes_type::element(m_cfes, std::forward<Cs>(cs)...);
    }
  };

  typename composite_finite_element_space<composite_finite_element<fe_pack...> >::element
  extend(const composite_finite_element_space<composite_finite_element<fe_pack...> >& m_cfes,
	 const submesh<cell_type, cell_type>& sm) const {
    return extend_impl<0, n_component, cfes_type>::call(sm, m_cfes, *this);
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



template<std::size_t n, std::size_t n_max, typename cfes_type>
struct to_composite_p1_impl {
  using cell_type = typename cfes_type::cell_type;
  
  template<typename ... Es>
  static typename cfes_type::element
  unwrap(const cfes_type& cfes, const mesh_data<double, fe_mesh<cell_type> >& data, Es&& ... es) {
    return to_composite_p1_impl<n + 1, n_max, cfes_type>::
             template unwrap<Es..., typename cfes_type::template fes_type<n>::element>(
               cfes,
               data,
               std::forward<Es>(es)...,
               to_p1_finite_element_function<cell_type>(
                 data[n],
                 cfes.template get_finite_element_space<n>()));
  }
};

template<std::size_t n, typename cfes_type>
struct  to_composite_p1_impl<n, n, cfes_type> {
  using cell_type = typename cfes_type::cell_type;
  
  template<typename ... Es>
  static typename cfes_type::element
  unwrap(const cfes_type& cfes, const mesh_data<double, fe_mesh<cell_type> >& data, Es&& ... es) {
    return typename cfes_type::element(cfes, std::forward<Es>(es)...);
  }
};


  
template<std::size_t n, typename cell_type>
typename composite_finite_element_space<
  make_composite_finite_element_t<n, typename cell_type::fe::lagrange_p1>
>::element
to_composite_p1_finite_element_function(const mesh_data<double, fe_mesh<cell_type> >& data,
                                        const composite_finite_element_space<
                                          make_composite_finite_element_t<n, typename cell_type::fe::lagrange_p1>
                                        >& cfes) {


  using cfes_type = composite_finite_element_space<
      make_composite_finite_element_t<n, typename cell_type::fe::lagrange_p1>
    >;

  return to_composite_p1_impl<0, n, cfes_type>::unwrap(cfes, data);

}


template<std::size_t n, std::size_t n_max, typename cfes_type>
struct to_composite_p0_impl {
  using cell_type = typename cfes_type::cell_type;
  
  template<typename ... Es>
  static typename cfes_type::element
  unwrap(const cfes_type& cfes, const mesh_data<double, fe_mesh<cell_type> >& data, Es&& ... es) {
    return to_composite_p0_impl<n + 1, n_max, cfes_type>::
             template unwrap<Es..., typename cfes_type::template fes_type<n>::element>(
               cfes,
               data,
               std::forward<Es>(es)...,
               to_p0_finite_element_function<cell_type>(
                 data[n],
                 cfes.template get_finite_element_space<n>()));
  }
};

template<std::size_t n, typename cfes_type>
struct  to_composite_p0_impl<n, n, cfes_type> {
  using cell_type = typename cfes_type::cell_type;
  
  template<typename ... Es>
  static typename cfes_type::element
  unwrap(const cfes_type& cfes, const mesh_data<double, fe_mesh<cell_type> >& data, Es&& ... es) {
    return typename cfes_type::element(cfes, std::forward<Es>(es)...);
  }
};


template<std::size_t n, typename cell_type>
typename composite_finite_element_space<
    make_composite_finite_element_t<n, typename cell_type::fe::lagrange_p0>
  >::element
to_composite_p0_finite_element_function(const mesh_data<double, fe_mesh<cell_type> >& data,
                                        const composite_finite_element_space<
                                          make_composite_finite_element_t<n, typename cell_type::fe::lagrange_p0>
                                        >& cfes) {
  using cfes_type = composite_finite_element_space<
      make_composite_finite_element_t<n, typename cell_type::fe::lagrange_p0>
    >;

  return to_composite_p1_impl<0, n, cfes_type>::unwrap(cfes, data);
}



template<std::size_t n, std::size_t n_max, typename result_type, typename element_type>
struct to_mesh_cell_impl {
  using cell_type = typename element_type::cell_type;
  using fe_type = get_element_at_t<n, typename element_type::fe_list>;
  
  static void
  unwrap(result_type& result, const element_type& v) {
    mesh_data<double, fe_mesh<cell_type> > comp(to_mesh_cell_data<fe_type>(v.template get_component<n>()));
    for (std::size_t k(0); k < v.get_mesh().get_element_number(); ++k)
      result.value(k, n) = comp.value(k, 0);

    to_mesh_cell_impl<n + 1, n_max, result_type, element_type>::unwrap(result, v);
  }
};

template<std::size_t n, typename result_type, typename element_type>
struct to_mesh_cell_impl<n, n, result_type, element_type> {
  static void unwrap(result_type& result, const element_type& v) {}
};

template<typename fe_type>
mesh_data<double, fe_mesh<typename fe_type::cell_type> >
to_mesh_cell_data(const typename composite_finite_element_space<fe_type>::element& v) {
  throw std::string("not implemented yet");
  
  using cell_type = typename fe_type::cell_type;
  using cfe_type = fe_type;
  using cfes_type = composite_finite_element_space<fe_type>;
  using result_type = mesh_data<double, fe_mesh<cell_type> >;
  using element_type = typename cfes_type::element;

  mesh_data<double, fe_mesh<cell_type> > result(v.get_mesh(),
                                             mesh_data_kind::cell,
                                             cfe_type::n_component);
  
  to_mesh_cell_impl<0, cfe_type::n_component, result_type, element_type>::unwrap(result, v);

  return result;
}


template<std::size_t n, std::size_t n_max, typename result_type, typename element_type>
struct to_mesh_vertex_impl {
  using cell_type = typename element_type::cell_type;
  using fe_type = get_element_at_t<n, typename element_type::fe_list>;
  
  static void
  unwrap(result_type& result, const element_type& v) {
    mesh_data<double, fe_mesh<cell_type> > comp(to_mesh_vertex_data<fe_type>(v.template get_component<n>()));
    for (std::size_t k(0); k < v.get_mesh().get_vertex_number(); ++k)
      result.value(k, n) = comp.value(k, 0);

    to_mesh_vertex_impl<n + 1, n_max, result_type, element_type>::unwrap(result, v);
  }
};

template<std::size_t n, typename result_type, typename element_type>
struct to_mesh_vertex_impl<n, n, result_type, element_type> {
  static void unwrap(result_type& result, const element_type& v) {}
};

template<typename fe_type>
mesh_data<double, fe_mesh<typename fe_type::cell_type> >
to_mesh_vertex_data(const typename composite_finite_element_space<fe_type>::element& v) {
  static_assert(is_continuous<fe_type>::value,
		"conversion to vertex data is only defined "
		"for continuous finite element spaces.");

  using cell_type = typename fe_type::cell_type;
  using cfe_type = fe_type;
  using cfes_type = composite_finite_element_space<fe_type>;
  using result_type = mesh_data<double, fe_mesh<cell_type> >;
  using element_type = typename cfes_type::element;

  mesh_data<double, fe_mesh<cell_type> > result(v.get_mesh(),
                                             mesh_data_kind::vertex,
                                             cfe_type::n_component);
  
  to_mesh_vertex_impl<0, cfe_type::n_component, result_type, element_type>::unwrap(result, v);

  return result;
}


#endif /* _COMPOSITE_FES_H_ */
