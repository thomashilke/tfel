#ifndef _EXPORT_H_
#define _EXPORT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <type_traits>

#include "quadrature.hpp"
#include "meta.hpp"
#include "mesh_data.hpp"


namespace exporter {

  template<typename mesh_type>
  void ascii(const std::string& filename,
	     const mesh_data<double, mesh_type>& data) {
    typedef typename mesh_type::cell_type cell_type;
    const mesh_type& m(data.get_mesh());
    
    std::ofstream file(filename.c_str(), std::ios::out);
    file.precision(12);


    switch (data.get_kind()) {
    case mesh_data_kind::vertex: {

    for (std::size_t k(0); k < m.get_vertex_number(); ++k) {
      for (std::size_t n(0); n < m.get_embedding_space_dimension(); ++n)
	file << m.get_vertices().at(k, n) << " ";
      
      
      for (unsigned int i(0); i < data.get_component_number(); ++i)
	file << data.get_values().at(k, i) << " ";
      file << '\n';
    }
      
    }
      break;
      
    case mesh_data_kind::cell: {

    for (std::size_t k(0); k < m.get_element_number(); ++k) {
      const auto bc(cell_type::barycenter(m.get_vertices(), m.get_elements(), k));
      for (std::size_t n(0); n < m.get_embedding_space_dimension(); ++n)
	file << bc.at(n) << " ";
      
      for (unsigned int i(0); i < data.get_component_number(); ++i)
	file << data.get_values().at(k, i) << " ";
      file << '\n';
    }
      
    }
      break;
    }
  }


  namespace ensight6_detail {


    inline
    const char* variable_section_item(mesh_data_kind t,
				      std::size_t n_component) {
      switch (n_component) {
      case 1:
	switch (t) {
	case mesh_data_kind::cell: return "scalar per element";
	case mesh_data_kind::vertex: return "scalar per node";
	}
      case 2:
      case 3:
	switch (t) {
	case mesh_data_kind::cell: return "vector per element";
	case mesh_data_kind::vertex: return "vector per node";
	}
      default:
	throw std::string("unsupported component number");
      }
    }
    
    void write_static_case_file(std::ostream& stream,
				const std::string& geometry_filename,
				const std::vector<std::string>& var_filenames,
				const std::vector<std::string>& var_names,
				const std::vector<std::pair<mesh_data_kind, std::size_t> >& var_types) {
      stream << "FORMAT\ntype: ensight\n\n";
      stream << "GEOMETRY\n";
      stream << "model: " << geometry_filename << "\n\n";
      stream << "VARIABLE\n";
      for (std::size_t i(0); i < var_names.size(); ++i) {
        stream << variable_section_item(var_types[i].first, var_types[i].second) << ": " << var_names[i]
               << " " << var_filenames[i] << '\n';

      }
    }
      
    template<typename cell_type>
    struct ensight_cell_name {
    private:
      using cell_list = type_list<cell::point, cell::edge, cell::triangle, cell::tetrahedron>;
      static constexpr const char* values[] = {"point", "bar2", "tria3", "tetra4"};
    public:
      static constexpr const char* value = values[get_index_of_element<cell_type, cell_list>::value];
    };


      

    template<typename mesh_type>
    void write_static_geometry_file(std::ostream& stream, const std::string& mesh_name, const std::string& part_name,
				    const mesh_type& m) {
      using cell_type = typename mesh_type::cell_type;
      
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
      stream << ensight_cell_name<cell_type>::value << '\n';
      stream << std::setw(8) << std::right << m.get_elements().get_size(0) << '\n';

      for (std::size_t n(0); n < m.get_elements().get_size(0); ++n) {
	for (std::size_t k(0); k < cell_type::n_vertex_per_element ; ++k)
	  stream << std::setw(8) << std::right << m.get_elements().at(n, k) + 1;
	stream << '\n';
      }
    }
    
    template<typename mesh_type, typename ... As>
    void write_static_geometry_file(std::ostream& stream, const std::string& mesh_name, const std::string& part_name,
				    const mesh_data<double, mesh_type>& data, As&& ...) {
      const mesh_type& m(data.get_mesh());
      write_static_geometry_file(stream, mesh_name, part_name, m);
    }


    template<typename mesh_type>
    void write_static_variable_file(std::ostream& stream,
				    const mesh_data<double, mesh_type>& data,
				    const std::string& var_name) {
      std::size_t fill(1);
      switch (data.get_component_number()) {
      case 1:
	fill = 1; break;
      case 2:
      case 3:
	fill = 3; break;
      default:
	throw std::string("write_static_variable_file: unsupported component number");
      }
      
      switch (data.get_kind()) {
      case mesh_data_kind::cell: {
	stream << var_name << '\n';
	stream << "part" << std::setw(8) << std::right << 1 << '\n';
	stream << ensight_cell_name<typename mesh_type::cell_type>::value;
	for (std::size_t k(0); k < data.get_values().get_size(0); ++k) {
          for (std::size_t n(0); n < data.get_component_number(); ++n) {
            std::size_t value_id(k * fill + n);
            if (value_id % 6 == 0)
              stream << '\n';
            stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
                   << data.value(k, n);
          }
	  for (std::size_t n(data.get_component_number()); n < fill; ++n) {
	    std::size_t value_id(k * fill + n);
	    if (value_id % 6 == 0)
              stream << '\n';
	    stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
                   << 0.0;
	  }
	}
      }
	break;
	
      case mesh_data_kind::vertex:  {
	stream << var_name;
	for (std::size_t k(0); k < data.get_values().get_size(0); ++k) {
          for (std::size_t n(0); n < data.get_component_number(); ++n) {
            std::size_t value_id(k * fill + n);
            if (value_id % 6 == 0)
              stream << '\n';
            stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
                   << data.value(k, n);
          }
	  for (std::size_t n(data.get_component_number()); n < fill; ++n) {
	    std::size_t value_id(k * fill + n);
	    if (value_id % 6 == 0)
              stream << '\n';
	    stream << std::setw(12) << std::right << std::setprecision(5) << std::scientific
                   << 0.0;
	  }
	}
      }
	break;
      }
    }

    void write_static_variable_file_dispatch(std::vector<std::string>::iterator filename_it,
					     std::vector<std::string>::iterator name_it,
					     std::vector<std::pair<mesh_data_kind, std::size_t> >::iterator vt_it,
					     const std::string& filename) {}

    template<typename mesh_type, typename ... As>
    void write_static_variable_file_dispatch(std::vector<std::string>::iterator filename_it,
					     std::vector<std::string>::iterator name_it,
					     std::vector<std::pair<mesh_data_kind, std::size_t> >::iterator vt_it,
					     const std::string& filename,
					     const mesh_data<double, mesh_type>& data,
					     const std::string& var_name,
					     As&& ... as) {
      const std::string variable_filename(filename + ".var." + var_name);
      std::ofstream variable_file(variable_filename.c_str(), std::ios::out);
      
      write_static_variable_file(variable_file, data, var_name);
      *filename_it = variable_filename;
      *name_it = var_name;
      *vt_it = std::make_pair(data.get_kind(), data.get_component_number());
      
      write_static_variable_file_dispatch(filename_it + 1, name_it + 1, vt_it + 1, filename, std::forward<As>(as)...);
    }

  }


  
  template<typename ... As>
  void ensight6(const std::string& filename, 
		As&& ... as) {
    const std::string
      case_filename(filename + ".case"),
      geometry_filename(filename + ".geom");

    std::vector<std::string> var_filenames(sizeof...(As)/2);
    std::vector<std::string> var_names(sizeof...(As)/2);
    std::vector<std::pair<mesh_data_kind, std::size_t> > var_types(sizeof...(As)/2);
    
    std::ofstream
      case_file(case_filename.c_str(), std::ios::out),
      geometry_file(geometry_filename.c_str(), std::ios::out);

    ensight6_detail::write_static_variable_file_dispatch(var_filenames.begin(),
							 var_names.begin(),
							 var_types.begin(),
							 filename, std::forward<As>(as)...);

    ensight6_detail::write_static_case_file(case_file,
					    geometry_filename,
					    var_filenames,
					    var_names,
					    var_types);
    
    ensight6_detail::write_static_geometry_file(geometry_file, "mesh", "part1", std::forward<As>(as)...);
  }

  template<typename mesh_type>
  void ensight6_geometry(const std::string& filename,
			 const mesh_type& m) {
    const std::string
      case_filename(filename + ".case"),
      geometry_filename(filename + ".geom");

    std::vector<std::string> var_filenames;
    std::vector<std::string> var_names;
    std::vector<std::pair<mesh_data_kind, std::size_t> > var_types;
    
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

  template<typename mesh_type>
  class ensight6_transient {
  public:
    using cell_type = typename mesh_type::cell_type;
    
    
    ensight6_transient(const std::string& filename,
		       const mesh_type& m,
		       const std::string& variable_name)
      : filename(filename),
	variable_name(variable_name),
	kind(mesh_data_kind::cell), n_component(0) {
      std::string geometry_filename(filename + ".geom");
      std::ofstream geometry_file(geometry_filename.c_str(), std::ios::out);
      ensight6_detail::write_static_geometry_file(geometry_file,
						  "mesh", "part1", m);
    }

    ~ensight6_transient() {
      std::string case_filename(filename + ".case");
      std::ofstream case_file(case_filename.c_str(), std::ios::out);
      write_transient_case_file(case_file);
    }
    
    void export_time_step(double time,
			  const mesh_data<double, mesh_type>& data) {
      std::ostringstream variable_filename;
      variable_filename << filename << ".var." << variable_name << "."
			<< std::setfill('0') << std::setw(6) << std::right
			<< times.size();
      
      std::ofstream variable_file(variable_filename.str().c_str(), std::ios::out);
      ensight6_detail::write_static_variable_file(variable_file, data, variable_name);

      times.push_back(time);
      kind = data.get_kind();
      n_component = data.get_component_number();
    }

  private:
    std::string filename;
    std::string variable_name;
    mesh_data_kind kind;
    std::size_t n_component;
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
      stream << ensight6_detail::variable_section_item(kind, n_component) << ": " << 1 << " "
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
