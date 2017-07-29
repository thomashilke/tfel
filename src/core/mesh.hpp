#ifndef _MESH_H_
#define _MESH_H_

#include <algorithm>
#include <map>
#include <vector>
#include <ostream>
#include <cassert>
#include <functional>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "vector_operation.hpp"


template<typename cell>
class mesh;


template<typename p_cell_type, typename cell_t = typename p_cell_type::boundary_cell_type>
class submesh {
public:
  typedef p_cell_type parent_cell_type;
  typedef cell_t cell_type;

  submesh(const mesh<parent_cell_type>& parent,
	  const array<unsigned int>& el,
	  const array<unsigned int>& el_id,
	  const array<unsigned int>& sd_id)
    : m(parent),
      elements(el),
      parent_element_id(el_id),
      parent_subdomain_id(sd_id) {}
  
  submesh(const mesh<parent_cell_type>& parent)
    : m(parent),
      elements{0, 0},
      parent_element_id{0},
      parent_subdomain_id{0} {}
  std::size_t get_embedding_space_dimension() const { return m.get_embedding_space_dimension(); }
  double get_cell_volume(std::size_t k) const { return cell_type::get_cell_volume(m.get_vertices(), elements, k); }
  std::size_t get_element_number() const { return elements.get_size(0); }
  array<double> get_jmt(std::size_t k) const { return m.get_jmt(parent_element_id.at(k)); }
  std::size_t get_subdomain_id(std::size_t k) const { return parent_subdomain_id.at(k); }
  std::size_t get_parent_element_id(std::size_t k) const { return parent_element_id.at(k); }

  const array<double>& get_vertices() const { return m.get_vertices(); }
  const array<unsigned int>& get_elements() const { return elements; }
  
  submesh<parent_cell_type> query_elements(const std::function<bool(const double*)>& f) const {
    std::vector<bool> selected_elements(get_element_number(), false);
    
    const array<double>& vertices(m.get_vertices());
    for (std::size_t k(0); k < get_element_number(); ++k) {
      for (std::size_t n(0); n < vertices.get_size(1); ++n)
	selected_elements.at(k) = selected_elements.at(k) || f(&vertices.at(elements.at(k, n), 0));
    }
    
    return submesh_from_selection(selected_elements);
  }

  submesh<parent_cell_type> inflow_boundary(const std::function<array<double>(const double*)>& b) const {
    std::vector<bool> selected_elements(get_element_number(), false);

    for (std::size_t k(0); k < get_element_number(); ++k)
      for (std::size_t n(0); n < get_vertices().get_size(1); ++n) {
	const auto normal(cell_type::normal(m.get_vertices(), elements, k));
	const auto b_val(b(&m.get_vertices().at(elements.at(k, n), 0)));
	const double dot_product(dotp(normal, b_val));
	const bool is_inflow_element(dot_product < 0.0);

	selected_elements.at(k) = selected_elements.at(k) || is_inflow_element;
      }


    return submesh_from_selection(selected_elements);
  }
  

  void show(std::ostream& stream) const {
    for (std::size_t k(0); k < get_element_number(); ++k) {
      stream << "element " << k << ": ";
      for (std::size_t i(0); i < elements.get_size(1); ++i) {
	stream << elements.at(k, i) << " ";
      }
      stream << " is #" << parent_subdomain_id.at(k) << "subdomain of parent element " << parent_element_id.at(k);
      stream << std::endl;
    }
  }
  
private:
  const mesh<parent_cell_type>& m;

  array<unsigned int> elements;
  array<unsigned int> parent_element_id;
  array<unsigned int> parent_subdomain_id;

private:
  submesh<parent_cell_type> submesh_from_selection(const std::vector<bool>& selected_elements) const {
    const std::size_t n_elements(std::count(selected_elements.begin(),
					    selected_elements.end(),
					    true));

    array<unsigned int> el_id{n_elements};
    array<unsigned int> sd_id{n_elements};
    array<unsigned int> el{n_elements, cell_type::n_vertex_per_element};

    for (std::size_t n(0), m(0); n < get_element_number(); ++n) {
      if (selected_elements[n]) {
	el_id.at(m) = parent_element_id.at(n);
	sd_id.at(m) = parent_subdomain_id.at(n);
	std::copy(&elements.at(n, 0),
		  &elements.at(n, 0) + cell_type::n_vertex_per_element,
		  &el.at(m, 0));
	++m;
      }
    }
    
    return submesh<parent_cell_type>(m, el, el_id, sd_id);
  }
  
};

template<typename cell>
class mesh {
public:
  typedef cell cell_type;

  mesh(const double* vertices,
       unsigned int n_vertices, unsigned int n_components,
       const unsigned int* elements, unsigned int n_elements)
    : vertices{n_vertices, n_components},
      elements{n_elements, cell_type::n_vertex_per_element},
      references{n_elements},
      h{n_elements} {
    (this->vertices).set_data(vertices);
    (this->elements).set_data(elements);
    references.fill(0);

    // We have to ensure that the vertices are ordered in the
    // elements list.
    if(not check_elements_admissibility())
      sort_elements();

    compute_element_diameter();
  }

  mesh(const double* vertices,
       unsigned int n_vertices, unsigned int n_components,
       const unsigned int* elements, unsigned int n_elements,
       const unsigned int* references)
    : mesh(vertices, n_vertices, n_components, elements, n_elements) {
    (this->references).set_data(references);
  }

  mesh(const submesh<cell_type, cell_type>& m)
    : vertices{0},
      elements{m.get_element_number(),
	       cell_type::n_vertex_per_element},
      references{m.get_element_number()},
      h{m.get_element_number()} {
    references.fill(0);

    std::set<unsigned int> selected_nodes;
    for (std::size_t k(0); k < get_element_number(); ++k) {
      for (std::size_t n(0); n < cell_type::n_vertex_per_element; ++n)
	selected_nodes.insert(m.get_elements().at(k, n));
    }

    {
      vertices = array<double>{selected_nodes.size(), cell_type::n_dimension};
      std::size_t new_node_id(0);
      for (const auto node_id: selected_nodes) {
	for (std::size_t k(0); k < cell_type::n_dimension; ++k)
	  vertices.at(new_node_id, k) = m.get_vertices().at(node_id, k);
	new_node_id += 1;
      }
    }

    {
      std::map<unsigned int, unsigned int> index_map;
      std::size_t new_node_id(0);
      for (const auto node_id: selected_nodes) {
	index_map[node_id] = new_node_id;
	new_node_id += 1;
      }

      for (std::size_t k(0); k < get_element_number(); ++k) {
	for (std::size_t n(0); n < cell_type::n_vertex_per_element; ++n) {
	  elements.at(k, n) = index_map[m.get_elements().at(k, n)];
	}
      }
    }
  }

  submesh<cell_type, cell_type> query_elements(const std::function<bool(const double*)>& f) const {
    std::vector<bool> selected_elements(get_element_number(), false);
    
    for (std::size_t k(0); k < get_element_number(); ++k) {
      for (std::size_t n(0); n < vertices.get_size(1); ++n)
	selected_elements.at(k) = selected_elements.at(k) || f(&vertices.at(elements.at(k, n), 0));
    }
    
    return submesh_from_selection(selected_elements);
  }
  
  std::size_t get_embedding_space_dimension() const { return vertices.get_size(1); }
  std::size_t get_element_number() const { return elements.get_size(0); }
  std::size_t get_vertex_number() const { return vertices.get_size(0); }
  
  const array<double>& get_vertices() const { return vertices; }
  const array<unsigned int>& get_elements() const { return elements; }

  double get_cell_volume(std::size_t k) const {
    return cell_type::get_cell_volume(vertices, elements, k);
  }

  array<double> get_jmt(std::size_t k) const {
    return cell_type::get_jmt(vertices, elements, k);
  }

  submesh<cell_type, cell_type> get_submesh_with_reference(std::size_t ref_id) {
    const std::size_t element_number(std::count(&references.at(0),
						&references.at(0) + references.get_size(0),
						ref_id));
    array<unsigned int> elements{element_number, cell_type::n_vertex_per_element};
    array<unsigned int> parent_element_id{element_number};
    array<unsigned int> parent_element_subdomain_id{element_number};
    parent_element_subdomain_id.fill(0);

    
    for (std::size_t k(0), pk(0); pk < get_element_number(); ++pk) {
      if (references.at(pk) == ref_id) {
	for (std::size_t n(0); n < cell_type::n_vertex_per_element; ++n)
	  elements.at(k, n) = this->elements.at(pk, n);
	  
	parent_element_id.at(k) = pk;
	
	k += 1;
      }
    }

    return submesh<cell_type, cell_type>(*this, elements, parent_element_id, parent_element_subdomain_id);
  }
  
  submesh<cell_type, ::cell::point> get_point_submesh(std::size_t k = 0) const {
    array<unsigned int> el_id{1}, sd_id{1}, el{1, 1};
    el_id.at(0) = k;
    sd_id.at(0) = 0;
    el.at(0, 0) = elements.at(k, 0);
    return  submesh<cell_type, ::cell::point>(*this, el, el_id, sd_id);
  }

  struct subdomain_info {
    std::size_t parent_element_id, parent_element_subdomain_id;
    ::cell::subdomain_type subdomain;

    bool operator<(const subdomain_info& op) const { return subdomain < op.subdomain; }
  };

  submesh<cell_type> get_boundary_submesh() const {
    using ::cell::subdomain_type;

    typedef std::map<subdomain_info, unsigned int> map_type;
    map_type subdomain_count;

    const std::size_t subdomain_id(cell_type::n_subdomain_type - 2);
    for (unsigned int k(0); k < get_element_number(); ++k) {
      for (unsigned int j(0); j < cell_type::n_subdomain(subdomain_id); ++j) {	
	subdomain_info current{k, j, cell::get_subdomain(elements, k, subdomain_id, j)};
	subdomain_count[current] += 1;
      }
    }
    const std::size_t n_elements(std::count_if(subdomain_count.begin(),
					       subdomain_count.end(),
					       [](const typename map_type::value_type& subdomain) {
						 return subdomain.second == 1;
					       }));

    array<unsigned int> el_id{n_elements};
    array<unsigned int> sd_id{n_elements};
    array<unsigned int> el{n_elements, cell_type::boundary_cell_type::n_vertex_per_element};

    unsigned int n(0);
    for (const auto& info: subdomain_count) {
      if (info.second == 1) {
	el_id.at(n) = info.first.parent_element_id;
	sd_id.at(n) = info.first.parent_element_subdomain_id;
	std::copy(info.first.subdomain.begin(),
		  info.first.subdomain.end(),
		  &el.at(n, 0));
	++n;
      }
    }
    
    return submesh<cell_type>(*this, el, el_id, sd_id);
  }

  void show(std::ostream& stream) const {
    for (std::size_t k(0); k < get_element_number(); ++k) {
      stream << "element " << k << ": ";
      for (std::size_t i(0); i < elements.get_size(1); ++i) {
	stream << elements.at(k, i) << " ";
      }
      stream << std::endl;
    }
  }

#warning "FIXME: implement a better than O(N) algorithm."
  std::size_t get_element_at(const double* x) const {
    for (std::size_t k(0); k < get_element_number(); ++k) {
      const array<double>
	bc_coord(cell_type::get_barycentric_coordinates(vertices, elements, k, x));

      const double* min(std::min_element(&bc_coord.at(0),
					 &bc_coord.at(0) + cell_type::n_vertex_per_element));
      if (*min >= 0.0)
	return k;
    }

    throw std::string("point out of bounds");
  }

  double get_h_max() const {
    return h_max;
  }

  double get_element_diameter(std::size_t k) const {
    return h.at(k);
  }
  
private:
  array<double> vertices;
  array<unsigned int> elements;
  array<unsigned int> references;
  array<double> h;
  double h_max;

  bool check_elements_admissibility() {
    for (unsigned int k(0); k < elements.get_size(0); ++k)
      for (unsigned int i(0); i < elements.get_size(1) - 1; ++i)
	if (elements.at(k, i) >= elements.at(k, i + 1))
	  return false;

    return true;
  }
  
  void sort_elements() {
    for (unsigned int k(0); k < elements.get_size(0); ++k)
      std::sort(&elements.at(k, 0),
		&elements.at(k, cell_type::n_vertex_per_element - 1));
  }

  void compute_element_diameter() {
    for (std::size_t k(0); k < get_element_number(); ++k)
      h.at(k) = cell_type::element_diameter(vertices, elements, k);


    h_max = *std::max_element(&h.at(0),
			      &h.at(0) + get_element_number());
  }

  submesh<cell_type, cell_type> submesh_from_selection(const std::vector<bool>& selected_elements) const {
    const std::size_t n_elements(std::count(selected_elements.begin(),
					    selected_elements.end(),
					    true));

    array<unsigned int> el_id{n_elements};
    array<unsigned int> sd_id{n_elements};
    array<unsigned int> el{n_elements, cell_type::n_vertex_per_element};

    for (std::size_t n(0), m(0); n < get_element_number(); ++n) {
      if (selected_elements[n]) {
	el_id.at(m) = n;
	sd_id.at(m) = 0;
	std::copy(&elements.at(n, 0),
		  &elements.at(n, 0) + cell_type::n_vertex_per_element,
		  &el.at(m, 0));
	++m;
      }
    }
    
    return submesh<cell_type, cell_type>(*this, el, el_id, sd_id);
  }
};

mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n);
mesh<cell::triangle> gen_square_mesh(double x_1, double x_2,
				     unsigned int n_1, unsigned int n_2);
mesh<cell::tetrahedron> gen_cube_mesh(double x_1, double x_2, double x_3,
				      unsigned int n_1, unsigned int n_2, unsigned int n_3);

#endif /* _MESH_H_ */
