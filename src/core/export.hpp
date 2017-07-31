#ifndef _EXPORT_H_
#define _EXPORT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <type_traits>

#include "fe.hpp"
#include "fes.hpp"
#include "quadrature.hpp"
#include "meta.hpp"

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
	file << m.get_vertices().at(vertex_id, i) << " " << v.get_coefficients().at(dof_id);
      file << '\n';
      ++dof_id;
    }
  }


  namespace ensight6_detail {
    enum class variable_type {nodal, elemental};
    
    template<typename fes_element_type>
    struct variable_type_from_element {
    private:
      using fe_list = type_list<finite_element::edge_lagrange_p0,
				finite_element::edge_lagrange_p1,
				finite_element::edge_lagrange_p1_bubble,
				finite_element::triangle_lagrange_p0,
				finite_element::triangle_lagrange_p1,
				finite_element::triangle_lagrange_p1_bubble,
				finite_element::tetrahedron_lagrange_p0,
				finite_element::tetrahedron_lagrange_p1>;
      static constexpr variable_type variable_types[] = {variable_type::elemental, variable_type::nodal, variable_type::nodal,
							 variable_type::elemental, variable_type::nodal, variable_type::nodal,
							 variable_type::elemental, variable_type::nodal};
      using fe_type = typename fes_element_type::fe_type;

    public:
      static constexpr variable_type value = variable_types[get_index_of_element<fe_type, fe_list>::value];
    };
    
    void write_static_case_file(std::ostream& stream,
				const std::string& geometry_filename,
				const std::vector<std::string>& var_filenames,
				const std::vector<std::string>& var_names,
				const std::vector<variable_type>& var_types) {
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      for (std::size_t i(0); i < var_names.size(); ++i) {
	if (var_types[i] == variable_type::nodal) {
	  stream << "scalar per node: " << var_names[i]
		 << " " << var_filenames[i] << '\n';
	} else if (var_types[i] == variable_type::elemental) {
	  stream << "scalar per element: " << var_names[i]
		 << " " << var_filenames[i] << '\n';
	}
      }
    }
      
    template<typename cell_type>
    struct ensight_element_type_name {
    private:
      using cell_list = type_list<cell::point, cell::edge, cell::triangle, cell::tetrahedron>;
      static constexpr const char* values[] = {"point", "bar2", "tria3", "tetra4"};
    public:
      static constexpr const char* value = values[get_index_of_element<cell_type, cell_list>::value];
    };
    
    template<>
    struct ensight_element_type_name<cell::edge> { static constexpr const char* value = "bar2"; };

    template<>
    struct ensight_element_type_name<cell::triangle> { static constexpr const char* value = "tria3"; };


    template<typename fe_type> struct variable_section_item;
    template<> struct variable_section_item<finite_element::triangle_lagrange_p0> { static constexpr const char* value = "scalar per element: "; };
    template<> struct variable_section_item<finite_element::triangle_lagrange_p1> { static constexpr const char* value = "scalar per node: "; };
    template<> struct variable_section_item<finite_element::tetrahedron_lagrange_p0> { static constexpr const char* value = "scalar per element: "; };
    template<> struct variable_section_item<finite_element::tetrahedron_lagrange_p1> { static constexpr const char* value = "scalar per node: "; };


    template<typename cell_type>
    void write_static_geometry_file(std::ostream& stream, const std::string& mesh_name, const std::string& part_name,
				    const mesh<cell_type>& m) {
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
    
    template<typename fes_element_type, typename ... As>
    void write_static_geometry_file(std::ostream& stream, const std::string& mesh_name, const std::string& part_name,
				    const fes_element_type& v, As&& ...) {
      using cell_type = typename fes_element_type::fe_type::cell_type;
      const mesh<cell_type>& m(v.get_mesh());
      write_static_geometry_file(stream, mesh_name, part_name, m);
    }


    template<typename fes_element_type>
    void write_static_variable_file(std::ostream& stream,
				    const fes_element_type& v,
				    const std::string& var_name);

    template<>
    void write_static_variable_file<finite_element_space<finite_element::edge_lagrange_p0>::element>
    (std::ostream& stream,
     const finite_element_space<finite_element::edge_lagrange_p0>::element& v,
     const std::string& var_name) {
      stream << var_name << '\n';
      stream << "part" << std::setw(8) << std::right << 1 << '\n';
      stream << ensight_element_type_name<cell::edge>::value;
      for (std::size_t k(0); k < v.get_coefficients().get_size(0); ++k) {
	if (k % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(v.get_finite_element_space().get_dof(k, 0));
      }
    }

    template<>
    void write_static_variable_file<finite_element_space<finite_element::triangle_lagrange_p0>::element>
    (std::ostream& stream,
     const typename finite_element_space<finite_element::triangle_lagrange_p0>::element& v,
     const std::string& var_name) {
      stream << var_name << '\n';
      stream << "part" << std::setw(8) << std::right << 1 << '\n';
      stream << ensight_element_type_name<cell::triangle>::value;
      for (std::size_t k(0); k < v.get_coefficients().get_size(0); ++k) {
	if (k % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(v.get_finite_element_space().get_dof(k, 0));
      }
    }


    
    template<>
    void write_static_variable_file<finite_element_space<finite_element::edge_lagrange_p1>::element>
    (std::ostream& stream,
     const typename finite_element_space<finite_element::edge_lagrange_p1>::element& v,
     const std::string& var_name) {
      stream << var_name;
      for (std::size_t n(0); n < v.get_coefficients().get_size(0); ++n) {
	if (n % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(n);
      }
    }

    template<>
    void write_static_variable_file<finite_element_space<finite_element::triangle_lagrange_p1>::element>
    (std::ostream& stream,
     const typename finite_element_space<finite_element::triangle_lagrange_p1>::element& v,
     const std::string& var_name) {
      stream << var_name;
      for (std::size_t n(0); n < v.get_coefficients().get_size(0); ++n) {
	if (n % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(n);
      }
    }


    template<>
    void write_static_variable_file<finite_element_space<finite_element::triangle_lagrange_p1_bubble>::element>
    (std::ostream& stream,
     const typename finite_element_space<finite_element::triangle_lagrange_p1_bubble>::element& v,
     const std::string& var_name) {
      stream << var_name;
      for (std::size_t n(0); n < v.get_mesh().get_vertex_number(); ++n) {
	if (n % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(n);
      }
    }


    template<>
    void write_static_variable_file<finite_element_space<finite_element::tetrahedron_lagrange_p0>::element>
    (std::ostream& stream,
     const typename finite_element_space<finite_element::tetrahedron_lagrange_p0>::element& v,
     const std::string& var_name) {
      stream << var_name << '\n';
      stream << "part" << std::setw(8) << std::right << 1 << '\n';
      stream << ensight_element_type_name<cell::tetrahedron>::value;
      for (std::size_t k(0); k < v.get_coefficients().get_size(0); ++k) {
	if (k % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(v.get_finite_element_space().get_dof(k, 0));
      }
    }


    template<>
    void write_static_variable_file<finite_element_space<finite_element::tetrahedron_lagrange_p1>::element>
    (std::ostream& stream,
     const typename finite_element_space<finite_element::tetrahedron_lagrange_p1>::element& v,
     const std::string& var_name) {
      stream << var_name;
      for (std::size_t n(0); n < v.get_coefficients().get_size(0); ++n) {
	if (n % 6 == 0)
	  stream << '\n';
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << v.get_coefficients().at(n);
      }
    }


    template<typename fes_element_type, typename ... As>
    void write_static_variable_file_dispatch(std::vector<std::string>::iterator filename_it,
					     std::vector<std::string>::iterator name_it,
					     std::vector<variable_type>::iterator vt_it,
					     const std::string& filename,
					     const fes_element_type& v, const std::string& var_name,
					     As&& ... as) {
      const std::string variable_filename(filename + ".var." + var_name);
      std::ofstream variable_file(variable_filename.c_str(), std::ios::out);
      
      write_static_variable_file(variable_file, v, var_name);
      *filename_it = variable_filename;
      *name_it = var_name;
      *vt_it = variable_type_from_element<fes_element_type>::value;
      
      write_static_variable_file_dispatch(filename_it + 1, name_it + 1, vt_it + 1, filename, std::forward<As>(as)...);
    }

    void write_static_variable_file_dispatch(std::vector<std::string>::iterator filename_it,
					     std::vector<std::string>::iterator name_it,
					     std::vector<variable_type>::iterator vt_it,
					     const std::string& filename) {}
  }


  
  template<typename ... As>
  void ensight6(const std::string& filename, 
		As&& ... as) {
    const std::string
      case_filename(filename + ".case"),
      geometry_filename(filename + ".geom");

    std::vector<std::string> var_filenames(sizeof...(As)/2);
    std::vector<std::string> var_names(sizeof...(As)/2);
    std::vector<ensight6_detail::variable_type> var_types(sizeof...(As)/2);
    
    std::ofstream
      case_file(case_filename.c_str(), std::ios::out),
      geometry_file(geometry_filename.c_str(), std::ios::out);

    ensight6_detail::write_static_variable_file_dispatch(var_filenames.begin(), var_names.begin(), var_types.begin(),
							 filename, std::forward<As>(as)...);

    ensight6_detail::write_static_case_file(case_file,
					    geometry_filename,
					    var_filenames,
					    var_names,
					    var_types);
    
    ensight6_detail::write_static_geometry_file(geometry_file, "mesh", "part1", std::forward<As>(as)...);
  }

  template<typename cell_type>
  void ensight6_geometry(const std::string& filename,
			 const mesh<cell_type>& m) {
    const std::string
      case_filename(filename + ".case"),
      geometry_filename(filename + ".geom");

    std::vector<std::string> var_filenames;
    std::vector<std::string> var_names;
    std::vector<ensight6_detail::variable_type> var_types;
    
    std::ofstream
      case_file(case_filename.c_str(), std::ios::out),
      geometry_file(geometry_filename.c_str(), std::ios::out);

    ensight6_detail::write_static_case_file(case_file,
					    geometry_filename,
					    var_filenames,
					    var_names,
					    var_types);
    
    ensight6_detail::write_static_geometry_file(geometry_file, "mesh", "part1", m);
  }

  template<typename fe>
  class ensight6_transient {
  public:
    typedef fe fe_type;
    typedef typename fe::cell_type cell_type;
    typedef finite_element_space<fe> fes_type;
    
    ensight6_transient(const std::string& filename,
		       const mesh<cell_type>& m,
		       const std::string& variable_name)
      : filename(filename),
	variable_name(variable_name) {
      std::string geometry_filename(filename + ".geom");
      std::ofstream geometry_file(geometry_filename.c_str(), std::ios::out);
      ensight6_detail::write_static_geometry_file<cell_type>(geometry_file,
							     "mesh", "part1", m);
    }

    ~ensight6_transient() {
      std::string case_filename(filename + ".case");
      std::ofstream case_file(case_filename.c_str(), std::ios::out);
      write_transient_case_file(case_file);
    }
    
    void export_time_step(double time,
			  const typename fes_type::element& var) {
      std::ostringstream variable_filename;
      variable_filename << filename << ".var." << variable_name << "."
			<< std::setfill('0') << std::setw(6) << std::right
			<< times.size();
      
      std::ofstream variable_file(variable_filename.str().c_str(), std::ios::out);
      ensight6_detail::write_static_variable_file(variable_file, var, variable_name);

      times.push_back(time);
    }

  private:
    std::string filename;
    std::string variable_name;
    std::vector<double> times;

  private:
    void write_transient_case_file(std::ostream& stream) {
      std::string geometry_filename(filename + ".geom");
      std::ostringstream variable_filename_pattern;
      variable_filename_pattern << filename << ".var." << variable_name << "."
				<< std::string(6, '*');
      
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << 1 << " " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      stream << ensight6_detail::variable_section_item<fe_type>::value << " " << 1 << " "
	     << variable_name
	     << " " << variable_filename_pattern.str() << "\n\n";

      stream << "TIME\n";
      stream << "time set: 1\n";
      stream << "number of steps: " << times.size() << '\n';
      stream << "filename start number: " << 0 << '\n';
      stream << "filename increment: " << 1 << '\n';
      stream << "time values: ";
      
      for (auto t: times)
	stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
	       << t << '\n';
    }
  };
}



#endif /* _EXPORT_H_ */
