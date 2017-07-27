#ifndef IMPORTER_H
#define IMPORTER_H

#include <alucelldb.hpp>

#include "../core/cell.hpp"
#include "../core/mesh.hpp"
#include "../core/fe.hpp"
#include "../core/fes.hpp"

namespace importer {
  namespace alucell {
    template<typename cell_type>
    mesh<cell_type> mesh(const std::string& db_filename,
			 const std::string& mesh_name) {
      ::alucell::database_read_access db(db_filename);
      ::alucell::database_index index(&db);

      const unsigned int
	mesh_nodes_id(index.get_variable_id(mesh_name + "_nodes")),
	mesh_elems_id(index.get_variable_id(mesh_name + "_nodes"));

      if (db.get_variable_type(mesh_nodes_id) != ::alucell::data_type::real_array
	  or db.get_variable_type(mesh_nodes_id) != ::alucell::data_type::int_array)
	throw std::string("importer::alucell::mesh: incompatible variable type");

      ::alucell::variable::array<double> n(&db, mesh_nodes_id);
      ::alucell::variable::array<int> e(&db, mesh_elems_id);

      if (n.get_components() != 3)
	throw std::string("importer::alucell::mesh: the nodes array is expected to have 3 components");

      if (e.get_components() != cell_type::n_vertex_per_element)
	throw std::string("importer::alucell::mesh: the element array is expected to have the same "
			  "number of components as the number of vertices per element.");
	  
      array<double> vertices{n.get_size(), cell_type::n_dimension};
      array<unsigned int> elements{e.get_size(), e.get_components()};

      for (unsigned int i(0); i < vertices.get_size(0); ++i)
	for (unsigned int j(0); j < vertices.get_size(1); ++j)
	  vertices.at(i, j) = n.get_value(i, j);
      
      for (unsigned int i(0); i < elements.get_size(0); ++i)
	for (unsigned int j(0); j < elements.get_size(1); ++j)
	  elements.at(i, j) = e.get_value(i, j);
      
      return mesh<cell_type>(&vertices.at(0,0), vertices.get_size(0), vertices.get_size(1),
			     &elements.at(0,0), elements.get_size(0));
    }
    
    template<typename fe_type>
    typename finite_element_space<fe_type>::element
    variable(const std::string& db_filename,
	     const std::string& mesh_name,
	     const std::string& var_name,
	     const finite_element_space<fe_type>& fes) {
      ::alucell::database_read_access db(db_filename);
      ::alucell::database_index index(&db);

      const unsigned int
	variable_id(index.get_variable_id(mesh_name + "_" + var_name));

      if (db.get_variable_type(variable_id) != ::alucell::data_type::real_array)
	throw std::string("");

      ::alucell::variable::array<double> v(&db, variable_id);

      if (v.get_components() != 1)
	throw std::string("");
      
      array<double> coefficients{v.get_size()};
      v.get_data(&coefficients.at(0));
      
      return typename finite_element_space<fe_type>::element(fes, coefficients);
    }

  }
}

#endif /* IMPORTER_H */
