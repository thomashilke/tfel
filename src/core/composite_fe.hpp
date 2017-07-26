#ifndef _COMPOSITE_FE_H_
#define _COMPOSITE_FE_H_

#include "meta.hpp"


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

#endif /* _COMPOSITE_FE_H_ */
