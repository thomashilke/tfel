#ifndef _ELEMENT_DIAMETER_H_
#define _ELEMENT_DIAMETER_H_

#include "mesh.hpp"
#include "fes.hpp"
#include "fe.hpp"


template<typename cell_type>
typename finite_element_space<typename cell_type::fe::lagrange_p0>::element
build_element_diameter_function(const fe_mesh<cell_type>& m,
				const finite_element_space<typename cell_type::fe::lagrange_p0>& fes) {
  array<double> coef{fes.get_dof_number()};
  for (std::size_t k(0); k < m.get_element_number(); ++k)
    coef.at(fes.get_dof(k, 0)) = m.get_element_diameter(k);
    
  return typename finite_element_space<typename cell_type::fe::lagrange_p0>::element(fes, coef);
}


#endif /* _ELEMENT_DIAMETER_H_ */

