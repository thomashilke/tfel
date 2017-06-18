#ifndef _MESH_H_
#define _MESH_H_

#include <algorithm>

#include <spikes/array.hpp>


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

#endif /* _MESH_H_ */
