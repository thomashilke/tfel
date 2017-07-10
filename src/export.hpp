#ifndef _EXPORT_H_
#define _EXPORT_H_

#include <fstream>

namespace exporter {
  void ascii(const std::string& filename,
	     const typename finite_element_space<finite_element::edge_lagrange_p1>::element& v) {
    std::ofstream file(filename.c_str(), std::ios::out);

    const mesh<cell::edge>& m(v.get_mesh());
    const std::set<cell::subdomain_type>& nodes(v.get_finite_element_space().get_subdomain_list().front());
    std::size_t dof_id(0);
    for (const auto& node: nodes) {
      unsigned int vertex_id(*node.begin());
      for (unsigned int i(0); i < m.get_embedding_space_dimension(); ++i)
	file << m.get_vertices().at(vertex_id, i) << " " << v.get_components().at(dof_id);
      file << '\n';
      ++dof_id;
    }
  }
}

#endif /* _EXPORT_H_ */
