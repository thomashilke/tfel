#ifndef FE_VALUE_MANAGER_H
#define FE_VALUE_MANAGER_H

#include <tuple>

#include <spikes/array.hpp>

#include "meta.hpp"


template<template<typename> class F, typename TL>
struct call_for_each;

template<template<typename> class F, typename T, typename ... Ts>
struct call_for_each<F, type_list<T, Ts...> > {
  template<typename ... As>
  static void call(As&& ... as) {
    F<T>::call(as...);
    call_for_each<F, type_list<Ts...> >::call(std::forward<As>(as)...);
  }
};

template<template<typename> class F>
struct call_for_each<F, type_list<> > {
  template<typename ... As>
  static void call(As&& ... as) {}
};
  

/*
 *  Metafunction: return the cell_type of a fe_type
 */
struct get_cell_type {
  template<typename fe_type>
  struct apply: is_type<typename fe_type::cell_type> {};
};


template<typename fe_list>
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
	    fe_pack::cell_type::n_dimension + 1}...),
    values_hat({n_quadrature_point,
	        fe_pack::n_dof_per_element,
	  fe_pack::cell_type::n_dimension + 1}...),
    xq_hat{n_quadrature_point, cell_type::n_dimension} {
    
    static_assert(list_size<cell_list>::value == 1,
		  "the cells types must be homogeneous");
  }

  void set_points(const array<double>& xq_hat) {
    this->xq_hat = xq_hat;

    using index_sequence = make_integral_list_t<std::size_t, sizeof...(fe_pack)>;
    call_for_each<prepare_hat_impl, index_sequence>::call(values, values_hat, xq_hat);
  }
    
  
  void prepare(const array<double>& jmt) {
    using index_sequence = make_integral_list_t<std::size_t, sizeof...(fe_pack)>;
    call_for_each<prepare_impl, index_sequence
		  >::call(values, values_hat, jmt, xq_hat);
  }

  void clear() {
    using index_sequence = make_integral_list_t<std::size_t, sizeof...(fe_pack)>;
    call_for_each<clear_impl, index_sequence>::call(values);
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
  values_type values_hat;

  array<double> xq_hat;

  template<typename integral_value>
  struct prepare_impl;
  
  template<std::size_t n>
  struct prepare_impl<integral_constant<std::size_t, n> > {
    static void call(values_type& values,
		     const values_type& values_hat,
		     const array<double>& jmt,
		     const array<double>& xq_hat) {
      
      array<double>& phi(std::get<n>(values));
      const array<double>& phi_hat(std::get<n>(values_hat));
    
      const std::size_t n_q(xq_hat.get_size(0));
      const std::size_t n_dim(cell_type::n_dimension);
      const std::size_t n_dof(fe_type<n>::n_dof_per_element);
    
      // prepare the basis function derivatives on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i) {
	  phi.at(q, i, 0) = phi_hat.at(q, i, 0);
	  for (std::size_t s(0); s < n_dim; ++s) {
	    phi.at(q, i, 1 + s) = 0.0;
	    for (std::size_t t(0); t < n_dim; ++t)
	      phi.at(q, i, 1 + s) += jmt.at(s, t) * phi_hat.at(q, i, 1 + t);
	  }
	}
      }
    }
  };


  template<typename integral_value>
  struct prepare_hat_impl;
  
  template<std::size_t n>
  struct prepare_hat_impl<integral_constant<std::size_t, n> > {
    static void call(values_type& values,
		     values_type& values_hat,
		     const array<double>& xq_hat) {

      array<double>& phi(std::get<n>(values));
      array<double>& phi_hat(std::get<n>(values_hat));
    
      const std::size_t n_q(xq_hat.get_size(0));
      const std::size_t n_dim(cell_type::n_dimension);
      const std::size_t n_dof(fe_type<n>::n_dof_per_element);
    
      // prepare the basis function on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i) {
	  phi_hat.at(q, i, 0) = fe_type<n>::phi(i, &xq_hat.at(q, 0));
	  phi.at(q, i, 0) = phi_hat.at(q, i, 0);
	  for (std::size_t s(0); s < n_dim; ++s) {
	    phi_hat.at(q, i, 1 + s) = fe_type<n>::dphi(s, i, &xq_hat.at(q, 0));
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
