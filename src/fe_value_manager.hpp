#ifndef FE_VALUE_MANAGER_H
#define FE_VALUE_MANAGER_H

#include <tuple>

#include "meta.hpp"


template<template<typename> class F, typename TL, typename ... As>
struct call_for_each;

template<template<typename> class F, typename T, typename ... Ts, typename ... As>
struct call_for_each<F, type_list<T, Ts...>, As...> {
  static void call(As ... as) {
    F<T>::call(as...);
    call_for_each<F, type_list<Ts...>, As...>::call(as...);
  }
};

template<template<typename> class F, typename ... As>
struct call_for_each<F, type_list<>, As...> {
  static void call(As ... as) {}
};
  

/*
 *  Metafunction: return the cell_type of a fe_type
 */
struct get_cell_type {
  template<typename fe_type>
  struct apply: is_type<typename fe_type::cell_type> {};
};


template<typename FE_LIST>
struct fe_value_manager;

template<typename ... fe_pack>
struct fe_value_manager<type_list<fe_pack...> > {
  using fe_list = type_list<fe_pack...>;
  using cell_list = unique_t<transform<get_cell_type, fe_list> >;
  using cell_type = get_element_at_t<0, cell_list>;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;
  
  fe_value_manager(std::size_t n_quadrature_point):
    values({n_quadrature_point,
	  fe_pack::n_dof_per_element,
	  fe_pack::cell_type::n_dimension + 1}...) {
    static_assert(list_size<cell_list>::value == 1,
		  "the cells types must be homogeneous");
  }

  void prepare(const array<double>& jmt,
	       const array<double>& xq) {
    using index_sequence = make_index_list_t<fe_list>;
    call_for_each<prepare_impl, index_sequence,
		  values_type&,
		  const array<double>&,
		  const array<double>&>::call(values, jmt, xq);
  }

  void clear() {
    using index_sequence = make_index_list_t<fe_list>;
    call_for_each<clear_impl, index_sequence, values_type&>::call(values);
  }
  
  template<std::size_t n>
  const array<double>& get_values() const {
    return std::get<n>(values);
  }
  
private:
  template<typename A, typename B> struct return_2nd: is_type<B> {};
  template<typename A, typename B> using return_2nd_t = typename return_2nd<A, B>::type;

  using values_type = std::tuple<return_2nd_t<fe_pack, array<double> >... >;
  
  values_type values;

  template<typename integral_value>
  struct prepare_impl;
  
  template<std::size_t n>
  struct prepare_impl<integral_constant<std::size_t, n> > {
    static void call(values_type& values,
		     const array<double>& jmt,
		     const array<double>& xq) {
      const std::size_t n_q(xq.get_size(0));
      
      array<double>& phi(std::get<n>(values));
    
      const std::size_t n_dim(cell_type::n_dimension);
      const std::size_t n_dof(fe_type<n>::n_dof_per_element);
    
      // prepare the basis function on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i)
	  phi.at(q, i, 0) = fe_type<n>::phi(i, &xq.at(q, 0));
      }
      
      // prepare the basis function derivatives on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i) {
	  for (std::size_t s(0); s < n_dim; ++s) {
	    phi.at(q, i, 1 + s) = 0.0;
	    for (std::size_t t(0); t < n_dim; ++t)
	      phi.at(q, i, 1 + s) += jmt.at(s, t) * fe_type<n>::dphi(t, i, &xq.at(q, 0));
	  }
	}
      }
    }
  };


  template<typename integral_value>
  struct clear_impl;

  template<std::size_t n>
  struct clear_impl<integral_constant<std::size_t, n> > {
    static void call(values_type& values) {
      array<double>& phi(std::get<n>(values));
      phi.fill(0.0);
    }
  };
};



#endif /* FE_VALUE_MANAGER_H */
