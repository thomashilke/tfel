#ifndef _ELEMENT_DIAMETER_H_
#define _ELEMENT_DIAMETER_H_

#include "mesh.hpp"
#include "fes.hpp"
#include "fe.hpp"


template<typename cell_type>
struct per_element_fe_type {
  using cell_list = type_list<cell::edge, cell::triangle, cell::tetrahedron>;
  using fe_list = type_list<finite_element::edge_lagrange_p0,
			    finite_element::triangle_lagrange_p0,
			    finite_element::tetrahedron_lagrange_p0>;

  using type = get_element_at_t<get_index_of_element<cell_type, cell_list>::value, fe_list>;
};

template<typename cell_type>
typename finite_element_space<typename per_element_fe_type<cell_type>::type>::element
build_element_diameter_function(const mesh<cell_type>& m,
				const finite_element_space<typename per_element_fe_type<cell_type>::type>& fes) {
  array<double> coef{fes.get_dof_number()};
  for (std::size_t k(0); k < m.get_element_number(); ++k)
    coef.at(fes.get_dof(k, 0)) = m.get_element_diameter(k);
    
  return typename finite_element_space<typename per_element_fe_type<cell_type>::type>::element(fes, coef);
}


#endif /* _ELEMENT_DIAMETER_H_ */

