#ifndef _COMPOSITE_FE_H_
#define _COMPOSITE_FE_H_

#include "meta.hpp"
#include "fe_value_manager.hpp"


/*
 *  A composite finite element is simply a wrapper around a typelist
 *  with an interface to access each of the finite element types in the list.
 */
template<typename ... fe_pack>
struct composite_finite_element {
  using fe_list = type_list<fe_pack...>;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;

  using cell_list = unique_t<transform<get_cell_type, fe_list> >;
  using cell_type = get_element_at_t<0, cell_list>;
  
  static const std::size_t n_component = list_size<fe_list>::value;
};

template<std::size_t n, typename U, typename ... Ts>
struct mk_composite_fe_impl: is_type<typename mk_composite_fe_impl<n - 1, U, Ts..., U>::type > {};

template<typename U, typename ... Ts>
struct mk_composite_fe_impl<0, U, Ts...>: is_type<composite_finite_element<Ts...> > {};

template<std::size_t n, typename fe_type>
struct make_composite_finite_element
  : is_type<
      typename mk_composite_fe_impl<
	         n, fe_type>::type > {};

template<std::size_t n, typename fe_type>
using make_composite_finite_element_t = typename make_composite_finite_element<n, fe_type>::type;

#endif /* _COMPOSITE_FE_H_ */
