#ifndef IMPORTER_H
#define IMPORTER_H

#include <alucelldb/alucelldb.hpp>


#include "../core/cell.hpp"
#include "../core/mesh.hpp"
#include "../core/fe.hpp"
#include "../core/fes.hpp"
#include "../core/composite_fe.hpp"
#include "../core/composite_fes.hpp"

namespace importer {
  namespace alucell {
    template<typename cell_type>
    mesh<cell_type> mesh(const std::string& db_filename,
			 const std::string& mesh_name) {
      ::alucell::database_read_access db(db_filename);
      ::alucell::database_index index(&db);

      const unsigned int
	mesh_nodes_id(index.get_variable_id(mesh_name + "_nodes")),
	mesh_elems_id(index.get_variable_id(mesh_name + "_elems")),
	mesh_refs_id(index.get_variable_id(mesh_name + "_refs"));
      
      if (db.get_variable_type(mesh_nodes_id) != ::alucell::data_type::real_array
	  or db.get_variable_type(mesh_elems_id) != ::alucell::data_type::element_array
	  or db.get_variable_type(mesh_refs_id) != ::alucell::data_type::int_array)
	throw std::string("importer::alucell::mesh: incompatible variable type");

      ::alucell::variable::array<double> n(&db, mesh_nodes_id);
      ::alucell::variable::array<int> e(&db, mesh_elems_id);
      ::alucell::variable::array<int> r(&db, mesh_refs_id);

      if (n.get_components() != 3)
	throw std::string("importer::alucell::mesh: the nodes array is expected to have 3 components");

      if (e.get_components() != cell_type::n_vertex_per_element)
	throw std::string("importer::alucell::mesh: the element array is expected to have the same "
			  "number of components as the number of vertices per element.");

      if (r.get_components() != 1)
	throw std::string();

      if (r.get_size() != e.get_size())
	throw std::string();
	  
      array<double> vertices{n.get_size(), cell_type::n_dimension};
      array<unsigned int> elements{e.get_size(), e.get_components()};
      array<unsigned int> references{e.get_size()};

      for (unsigned int i(0); i < vertices.get_size(0); ++i)
	for (unsigned int j(0); j < vertices.get_size(1); ++j)
	  vertices.at(i, j) = n.get_value(i, j);
      
      for (unsigned int i(0); i < elements.get_size(0); ++i) {
	references.at(i) = r.get_value(i, 0);
	for (unsigned int j(0); j < elements.get_size(1); ++j) 
	  elements.at(i, j) = e.get_value(i, j) - 1;
      }
      
      return ::mesh<cell_type>(&vertices.at(0,0), vertices.get_size(0), vertices.get_size(1),
			       &elements.at(0,0), elements.get_size(0),
			       &references.at(0));
    }




    template<typename fe_type>
    struct variable_impl {
      using result_type = typename finite_element_space<fe_type>::element;

      static
      result_type call(const std::string& db_filename,
		       const std::string& mesh_name,
		       const std::string& var_name,
		       const finite_element_space<fe_type>& fes) {
	static_assert(std::is_same<fe_type, finite_element::tetrahedron_lagrange_p1>::value or
		      std::is_same<fe_type, finite_element::triangle_lagrange_p1>::value or
		      std::is_same<fe_type, finite_element::edge_lagrange_p1>::value, "");

	/*
	 *  P0 field are not yet supported, since they require a permutation of the coefficients.
	 */
	
	::alucell::database_read_access db(db_filename);
	::alucell::database_index index(&db);

	const unsigned int
	  variable_id(index.get_variable_id(mesh_name + "_" + var_name));

	if (db.get_variable_type(variable_id) != ::alucell::data_type::real_array)
	  throw std::string("importer::alucell::variable: only real variables are supported");

	::alucell::variable::array<double> v(&db, variable_id);

	if (v.get_components() != 1)
	  throw std::string("importer::alucell::variable: only scalar variables are supported");
      
	array<double> coefficients{v.get_size()};
	v.get_data(&coefficients.at(0));
      
	return typename finite_element_space<fe_type>::element(fes, coefficients);
      }
    };

    template<typename ... fe_pack>
    struct variable_impl<composite_finite_element<fe_pack...> > {
      using cfe_type = composite_finite_element<fe_pack...>;
      using cfes_type = composite_finite_element_space<cfe_type>;

      using fe_list = typename cfe_type::fe_list;
      using unique_fe_list = unique_t<fe_list>;
      using fe_type = get_element_at_t<0, unique_fe_list>;
      
      using result_type = typename cfes_type::element;
      
      static
      result_type call(const std::string& db_filename,
		       const std::string& mesh_name,
		       const std::string& var_name,
		       const composite_finite_element_space<cfe_type>& cfes) {
	static_assert(list_size<unique_fe_list>::value == 1 and
		      (std::is_same<fe_type, finite_element::tetrahedron_lagrange_p1>::value or
		       std::is_same<fe_type, finite_element::triangle_lagrange_p1>::value or
		       std::is_same<fe_type, finite_element::edge_lagrange_p1>::value), "");

	/*
	 *  P0 field are not yet supported, since they require a permutation of the coefficients.
	 */
	
	::alucell::database_read_access db(db_filename);
	::alucell::database_index index(&db);

	const unsigned int
	  variable_id(index.get_variable_id(mesh_name + "_" + var_name));

	if (db.get_variable_type(variable_id) != ::alucell::data_type::real_array)
	  throw std::string("importer::alucell::variable: only real variables are supported");

	::alucell::variable::array<double> v(&db, variable_id);

	if (v.get_components() != cfe_type::n_component)
	  throw std::string("importer::alucell::variable: wrong number of components");
	
      
	array<double>
	  data{v.get_size(), cfe_type::n_component},
	  coefficients{v.get_size() * cfe_type::n_component};
	
        v.get_data(&data.at(0, 0));

	for (std::size_t n(0); n < cfe_type::n_component; ++n)
	  for (std::size_t k(0); k < v.get_size(); ++k)
	    coefficients.at(n * v.get_size() + k) = data.at(k, n);

	return result_type(cfes, coefficients);
      }
    };


    template<typename fe_type, typename fes_type>
    typename variable_impl<fe_type>::result_type
    variable(const std::string& db_filename,
	     const std::string& mesh_name,
	     const std::string& var_name,
	     const fes_type& fes) {
      return variable_impl<fe_type>::call(db_filename, mesh_name, var_name, fes);
    }

  }
}

#endif /* IMPORTER_H */
