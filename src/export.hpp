#ifndef _EXPORT_H_
#define _EXPORT_H_

#include <fstream>
#include <type_traits>

#include "fe.hpp"
#include "fes.hpp"
#include "quadrature.hpp"

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


  namespace ensight6_detail {
    template<typename fe>
    void write_static_case_file(std::ostream& stream,
				const std::string& geometry_filename,
				const std::string& variable_name,
				const std::string& variable_filename);

    template<>
    void write_static_case_file<finite_element::edge_lagrange_p1>(std::ostream& stream,
								  const std::string& geometry_filename,
								  const std::string& variable_name,
								  const std::string& variable_filename) {
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      stream << "scalar per node: " << variable_name
	     << " " << variable_filename << '\n';
    }

    template<>
    void write_static_case_file<finite_element::edge_lagrange_p0>(std::ostream& stream,
 								  const std::string& geometry_filename,
								  const std::string& variable_name,
								  const std::string& variable_filename) {
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      stream << "scalar per element: " << variable_name
	     << " " << variable_filename << '\n';
    }

    template<>
    void write_static_case_file<finite_element::triangle_lagrange_p1>(std::ostream& stream,
								      const std::string& geometry_filename,
								      const std::string& variable_name,
								      const std::string& variable_filename) {
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      stream << "scalar per node: " << variable_name
	     << " " << variable_filename << '\n';

    }

    template<>
    void write_static_case_file<finite_element::triangle_lagrange_p0>(std::ostream& stream,
								      const std::string& geometry_filename,
								      const std::string& variable_name,
								      const std::string& variable_filename) {
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      stream << "scalar per element: " << variable_name
	     << " " << variable_filename << '\n';

    }

    template<typename cell_type>
    struct ensight_element_type_name;
    
    template<>
    struct ensight_element_type_name<cell::edge> { static constexpr const char* value = "bar2"; };

    template<>
    struct ensight_element_type_name<cell::triangle> { static constexpr const char* value = "tria3"; };

    
    template<typename cell_type>
    void write_static_geometry_file(std::ostream& stream, const mesh<cell_type>& m, const std::string& mesh_name, const std::string& part_name) {
      stream << mesh_name << "\n\n";
      stream << "node id off\n";
      stream << "element id off\n";
      stream << "coordinates\n";
      stream << std::setw(8) << std::right << m.get_vertices().get_size(0) << '\n';

      for (std::size_t n(0); n < m.get_vertices().get_size(0); ++n) {
	for (std::size_t k(0); k < m.get_vertices().get_size(1); ++k) {
	  stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific << m.get_vertices().at(n, k);
	}
	for (std::size_t k(m.get_vertices().get_size(1)); k < 3; ++k) {
	  stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific << 0.0;
	}
	stream << '\n';
      }

      stream << "part" << std::setw(8) << std::right << 1 << '\n';
      stream << part_name << '\n';
      stream << ensight_element_type_name<cell_type>::value << '\n';
      stream << std::setw(8) << std::right << m.get_elements().get_size(0) << '\n';

      for (std::size_t n(0); n < m.get_elements().get_size(0); ++n) {
	for (std::size_t k(0); k < cell_type::n_vertex_per_element ; ++k)
	  stream << std::setw(8) << std::right << m.get_elements().at(n, k) + 1;
	stream << '\n';
      }
    }

    template<typename fe>
    void write_static_variable_file(std::ostream& stream,
				    const typename finite_element_space<fe>::element& v,
				    const std::string& var_name);

    template<>
    void write_static_variable_file<finite_element::edge_lagrange_p0>(std::ostream& stream,
								      const typename finite_element_space<finite_element::edge_lagrange_p0>::element& v,
								      const std::string& var_name) {
      stream << var_name;
      for (std::size_t k(0); k < v.get_components().get_size(0); ++k) {
	if (k % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_components().at(v.get_finite_element_space().get_dof(k, 0));
      }
    }

    template<>
    void write_static_variable_file<finite_element::triangle_lagrange_p0>(std::ostream& stream,
									  const typename finite_element_space<finite_element::triangle_lagrange_p0>::element& v,
									  const std::string& var_name) {
      stream << var_name;
      for (std::size_t k(0); k < v.get_components().get_size(0); ++k) {
	if (k % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_components().at(v.get_finite_element_space().get_dof(k, 0));
      }
    }


    
    template<>
    void write_static_variable_file<finite_element::edge_lagrange_p1>(std::ostream& stream,
								      const typename finite_element_space<finite_element::edge_lagrange_p1>::element& v,
								      const std::string& var_name) {
      stream << var_name;
      for (std::size_t n(0); n < v.get_components().get_size(0); ++n) {
	if (n % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_components().at(n);
      }
    }

    template<>
    void write_static_variable_file<finite_element::triangle_lagrange_p1>(std::ostream& stream,
									  const typename finite_element_space<finite_element::triangle_lagrange_p1>::element& v,
									  const std::string& var_name) {
      stream << var_name;
      for (std::size_t n(0); n < v.get_components().get_size(0); ++n) {
	if (n % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_components().at(n);
      }
    }
  }
  
  template<typename fe>
  void ensight6(const std::string& filename,
		const typename finite_element_space<fe>::element& v,
		const std::string& var_name) {
    const std::string
      case_filename(filename + ".case"),
      geometry_filename(filename + ".geom"),
      variable_filename(filename + ".var." + var_name);
    
    std::ofstream
      case_file(case_filename.c_str(), std::ios::out),
      geometry_file(geometry_filename.c_str(), std::ios::out),
      variable_file(variable_filename.c_str(), std::ios::out);

    ensight6_detail::write_static_case_file<fe>(case_file,
						geometry_filename,
						var_name,
						variable_filename);
    ensight6_detail::write_static_geometry_file<typename fe::cell_type>(geometry_file,
									v.get_mesh(),
									"mesh", "part1");
    ensight6_detail::write_static_variable_file<fe>(variable_file, v, var_name);
  }
}



#endif /* _EXPORT_H_ */
