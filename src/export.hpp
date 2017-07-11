#ifndef _EXPORT_H_
#define _EXPORT_H_

#include <fstream>
#include <type_traits>

namespace exporter {

  template<typename fe>
  void ascii(const std::string& filename,
	     const typename finite_element_space<fe>::element& v) {
    typedef fe fe_type;
    typedef typename fe_type::cell_type cell_type;
    static_assert(std::is_same<cell_type, cell::edge>::value, "ASCII export is only defined for variables defined on edge meshes.");
    static_assert(finite_element::is_continuous<fe_type>::value, "ASCII export is only defined for continuous finite element space elements.");
    static_assert(fe_type::n_dof_per_subdomain(0) > 0, "ASCII export is only defined for finite element spaces with dof associated to the mesh's nodes.");
    
    std::ofstream file(filename.c_str(), std::ios::out);
    file.precision(12);

    const mesh<cell_type>& m(v.get_mesh());
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
