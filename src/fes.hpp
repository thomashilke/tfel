#ifndef _FES_H_
#define _FES_H_

#include "mesh.hpp"

template<typename fe>
class finite_element_space {
public:
  typedef fe fe_type;
  typedef typename fe_type::cell_type cell_type;

  finite_element_space(const mesh<cell_type>& m)
    : dof_map{m.get_element_number(),
              fe_type::n_node_per_element} {}

private:
  array<unsigned int> dof_map;
};

#endif /* _FES_H_ */
