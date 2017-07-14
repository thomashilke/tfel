#ifndef _MESH_H_
#define _MESH_H_

#include <algorithm>
#include <map>
#include <vector>
#include <ostream>

#include <spikes/array.hpp>

#include "cell.hpp"

template<typename cell>
class mesh;

template<typename cell>
class submesh {
public:
  typedef cell parent_cell_type;
  typedef typename cell::boundary_cell_type cell_type;

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
  
  template<typename F>
  submesh<parent_cell_type> query_elements(F f) const {
    std::vector<bool> selected_elements(get_element_number(), false);
    
    const array<double>& vertices(m.get_vertices());
    for (std::size_t k(0); k < get_element_number(); ++k) {
      for (std::size_t n(0); n < vertices.get_size(1); ++n)
	selected_elements.at(k) = selected_elements.at(k) || f(&vertices.at(elements.at(k, n), 0));
    }

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
};

template<typename cell>
class mesh {
public:
  typedef cell cell_type;

  mesh(const double* vertices,
       unsigned int n_vertices, unsigned int n_components,
       const unsigned int* elements, unsigned int n_elements)
    : vertices{n_vertices, n_components},
      elements{n_elements, cell_type::n_vertex_per_element} {
    (this->vertices).set_data(vertices);
    (this->elements).set_data(elements);

    // We have to ensure that the vertices are ordered in the
    // elements list.
    if(not check_elements_admissibility())
      sort_elements();
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

  unsigned int get_element_at(const double* x) const {
    // TODO
  }
  
private:
  array<double> vertices;
  array<unsigned int> elements;

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
		&elements.at(k, cell::n_vertex_per_element - 1));
  }
};

mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n);
mesh<cell::triangle> gen_square_mesh(double x_1, double x_2,
				     unsigned int n_1, unsigned int n_2);

#endif /* _MESH_H_ */
