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
  for (std::size_t k(0); k < m.get_cell_number(); ++k)
    coef.at(fes.get_dof(k, 0)) = m.get_cell_diameter(k);
    
  return typename finite_element_space<typename cell_type::fe::lagrange_p0>::element(fes, coef);
}


template<typename mesh_type>
mesh_data<double, mesh_type> cell_diameter(const mesh_type& m) {
  array<double> h{m.get_cell_number(), 1};
  for (std::size_t k(0); k < m.get_cell_number(); ++k)
    h.at(k, 1) = m.get_cell_diameter(k);

  return mesh_data<double, mesh_type>(m, mesh_data_kind::cell, std::move(h));
}

#endif /* _ELEMENT_DIAMETER_H_ */

