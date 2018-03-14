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


template<typename value_t, typename mesh_t>
class mesh_data;

template<typename cell_t>
class fe_mesh;


template<typename p_cell_type, typename cell_t = typename p_cell_type::boundary_cell_type>
class submesh {
public:
  typedef p_cell_type parent_cell_type;
  typedef cell_t cell_type;

  submesh(const fe_mesh<parent_cell_type>& parent,
	  const array<unsigned int>& el,
	  const array<unsigned int>& el_id,
	  const array<unsigned int>& sd_id)
    : m(parent),
      cells(el),
      parent_cell_id(el_id),
      parent_subdomain_id(sd_id) {}
  
  submesh(const fe_mesh<parent_cell_type>& parent)
    : m(parent),
      cells{0, 0},
      parent_cell_id{0},
      parent_subdomain_id{0} {}

  std::size_t get_embedding_space_dimension() const { return m.get_embedding_space_dimension(); }
  double get_cell_volume(std::size_t k) const { return cell_type::get_cell_volume(m.get_vertices(), cells, k); }
  std::size_t get_cell_number() const { return cells.get_size(0); }
  std::size_t get_vertex_number() const { return m.get_vertices().get_size(0); }
  array<double> get_jmt(std::size_t k) const { return m.get_jmt(parent_cell_id.at(k)); }
  std::size_t get_subdomain_id(std::size_t k) const { return parent_subdomain_id.at(k); }
  std::size_t get_parent_cell_id(std::size_t k) const { return parent_cell_id.at(k); }
  const fe_mesh<parent_cell_type>& get_mesh() const {return m;}

  const array<double>& get_vertices() const { return m.get_vertices(); }
  const array<unsigned int>& get_cells() const { return cells; }
  
  submesh<parent_cell_type> query_cells(const std::function<bool(const double*)>& f) const {
    std::vector<bool> selected_cells(get_cell_number(), false);
    
    const array<double>& vertices(m.get_vertices());
    for (std::size_t k(0); k < get_cell_number(); ++k) {
      for (std::size_t n(0); n < vertices.get_size(1); ++n)
	selected_cells.at(k) = selected_cells.at(k) || f(&vertices.at(cells.at(k, n), 0));
    }
    
    return submesh_from_selection(selected_cells);
  }

  submesh<parent_cell_type> inflow_boundary(const std::function<array<double>(const double*)>& b) const {
    std::vector<bool> selected_cells(get_cell_number(), false);

    for (std::size_t k(0); k < get_cell_number(); ++k)
      for (std::size_t n(0); n < get_vertices().get_size(1); ++n) {
	const auto normal(cell_type::normal(m.get_vertices(), cells, k));
	const auto b_val(b(&m.get_vertices().at(cells.at(k, n), 0)));
	const double dot_product(dotp(normal, b_val));
	const bool is_inflow_cell(dot_product < 0.0);

	selected_cells.at(k) = selected_cells.at(k) || is_inflow_cell;
      }


    return submesh_from_selection(selected_cells);
  }
  

  void show(std::ostream& stream) const {
    for (std::size_t k(0); k < get_cell_number(); ++k) {
      stream << "cell " << k << ": ";
      for (std::size_t i(0); i < cells.get_size(1); ++i) {
	stream << cells.at(k, i) << " ";
      }
      stream << " is #" << parent_subdomain_id.at(k) << " subdomain of parent cell " << parent_cell_id.at(k);
      stream << std::endl;
    }
  }

  submesh<parent_cell_type, cell_type> submesh_from_selection(const mesh_data<bool, submesh<parent_cell_type, cell_type> >& cell_selection) const {
    std::vector<bool> selection(&cell_selection.value(0, 0), &cell_selection.value(0, 0) + get_cell_number());
    return submesh_from_selection(selection);
  }
  
private:
  const fe_mesh<parent_cell_type>& m;

  array<unsigned int> cells;
  array<unsigned int> parent_cell_id;
  array<unsigned int> parent_subdomain_id;

private:
  submesh<parent_cell_type> submesh_from_selection(const std::vector<bool>& selected_cells) const {
    const std::size_t n_cells(std::count(selected_cells.begin(),
                                         selected_cells.end(),
                                         true));

    array<unsigned int> el_id{n_cells};
    array<unsigned int> sd_id{n_cells};
    array<unsigned int> el{n_cells, cell_type::n_vertex_per_cell};

    for (std::size_t n(0), m(0); n < get_cell_number(); ++n) {
      if (selected_cells[n]) {
	el_id.at(m) = parent_cell_id.at(n);
	sd_id.at(m) = parent_subdomain_id.at(n);
	std::copy(&cells.at(n, 0),
		  &cells.at(n, 0) + cell_type::n_vertex_per_cell,
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
       const unsigned int* cells, unsigned int n_cells)
    : vertices{n_vertices, n_components},
      cells{n_cells, cell_type::n_vertex_per_cell},
      references{n_cells} {
        (this->vertices).set_data(vertices);
        (this->cells).set_data(cells);
        references.fill(0);

        // We have to ensure that the vertices are ordered in the
        // cells list.
        if(not check_cells_admissibility())
          sort_cells();
      }

  mesh(const double* vertices,
       unsigned int n_vertices, unsigned int n_components,
       const unsigned int* cells, unsigned int n_cells,
       const unsigned int* references)
    : mesh(vertices, n_vertices, n_components, cells, n_cells) {
    (this->references).set_data(references);
  }

  template<typename parent_cell_type>
  mesh(const submesh<parent_cell_type, cell_type>& m)
    : vertices{0},
      cells{m.get_cell_number(),
          cell_type::n_vertex_per_cell},
      references{m.get_cell_number()} {
        references.fill(0);

        std::set<unsigned int> selected_nodes;
        for (std::size_t k(0); k < get_cell_number(); ++k) {
          for (std::size_t n(0); n < cell_type::n_vertex_per_cell; ++n)
            selected_nodes.insert(m.get_cells().at(k, n));
        }

        {
          vertices = array<double>{selected_nodes.size(), parent_cell_type::n_dimension};
          std::size_t new_node_id(0);
          for (const auto node_id: selected_nodes) {
            for (std::size_t k(0); k < parent_cell_type::n_dimension; ++k)
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

          for (std::size_t k(0); k < get_cell_number(); ++k) {
            for (std::size_t n(0); n < cell_type::n_vertex_per_cell; ++n) {
              cells.at(k, n) = index_map[m.get_cells().at(k, n)];
            }
          }
        }
      }

  std::size_t get_embedding_space_dimension() const { return vertices.get_size(1); }
  std::size_t get_cell_number() const { return cells.get_size(0); }
  std::size_t get_vertex_number() const { return vertices.get_size(0); }
  
  const array<double>& get_vertices() const { return vertices; }
  const array<unsigned int>& get_cells() const { return cells; }


  struct subdomain_info {
    std::size_t parent_cell_id, parent_cell_subdomain_id;
    ::cell::subdomain_type subdomain;

    bool operator<(const subdomain_info& op) const { return subdomain < op.subdomain; }
  };

  void show(std::ostream& stream) const {
    for (std::size_t k(0); k < get_cell_number(); ++k) {
      stream << "cell " << k << ": ";
      for (std::size_t i(0); i < cells.get_size(1); ++i) {
	stream << cells.at(k, i) << " ";
      }
      stream << std::endl;
    }
  }

#warning "FIXME: implement a better than O(N) algorithm."
  std::size_t get_cell_at(const double* x) const {
    for (std::size_t k(0); k < get_cell_number(); ++k) {
      const array<double>
	bc_coord(cell_type::get_barycentric_coordinates(vertices, cells, k, x));

      const double* min(std::min_element(&bc_coord.at(0),
					 &bc_coord.at(0) + cell_type::n_vertex_per_cell));
      if (*min >= 0.0)
	return k;
    }

    throw std::string("point out of bounds");
  }

protected:
  mesh(): vertices{}, cells{}, references{} {}
  
  array<double> vertices;
  array<unsigned int> cells;
  array<unsigned int> references;

private:
  bool check_cells_admissibility() {
    for (unsigned int k(0); k < cells.get_size(0); ++k)
      for (unsigned int i(0); i < cells.get_size(1) - 1; ++i)
	if (cells.at(k, i) >= cells.at(k, i + 1))
	  return false;

    return true;
  }
  
  void sort_cells() {
    for (unsigned int k(0); k < cells.get_size(0); ++k)
      std::sort(&cells.at(k, 0),
		&cells.at(k, cell_type::n_vertex_per_cell - 1));
  }
};


template<typename cell>
class fe_mesh: public mesh<cell> {
  using mesh<cell>::vertices;
  using mesh<cell>::cells;
  using mesh<cell>::references;
  
public:
  typedef cell cell_type;

  fe_mesh(const double* vertices,
          unsigned int n_vertices, unsigned int n_components,
          const unsigned int* cells, unsigned int n_cells)
    : mesh<cell>(vertices, n_vertices, n_components, cells, n_cells),
      cell_volume{n_cells},
      h{n_cells} {
      compute_cell_diameter();
      compute_cell_volume();
      compute_jmt();
    }

  fe_mesh(const double* vertices,
          unsigned int n_vertices, unsigned int n_components,
          const unsigned int* cells, unsigned int n_cells,
          const unsigned int* references)
    : mesh<cell>(vertices, n_vertices, n_components, cells, n_cells, references),
      cell_volume{n_cells}, h{n_cells}, h_max(0.0) {
    compute_cell_diameter();
    compute_cell_volume();
    compute_jmt();    
  }

  template<typename parent_cell_type>
  fe_mesh(const submesh<parent_cell_type, cell_type>& m)
    : mesh<cell>(),
      cell_volume{m.get_cell_number()},
      h{m.get_cell_number()},
      h_max(0.0),
      jmt() {


      mesh<cell>::cells = array<unsigned int>{m.get_cell_number(),
                                              cell_type::n_vertex_per_cell};
      references = array<unsigned int>{m.get_cell_number()};
      references.fill(0);

      std::set<unsigned int> selected_nodes;
      for (std::size_t k(0); k < this->get_cell_number(); ++k) {
        for (std::size_t n(0); n < cell_type::n_vertex_per_cell; ++n)
          selected_nodes.insert(m.get_cells().at(k, n));
      }

      {
        vertices = array<double>{selected_nodes.size(), parent_cell_type::n_dimension};
        std::size_t new_node_id(0);
        for (const auto node_id: selected_nodes) {
          for (std::size_t k(0); k < parent_cell_type::n_dimension; ++k)
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

        for (std::size_t k(0); k < this->get_cell_number(); ++k) {
          for (std::size_t n(0); n < cell_type::n_vertex_per_cell; ++n) {
            mesh<cell>::cells.at(k, n) = index_map[m.get_cells().at(k, n)];
          }
        }
      }

      compute_cell_diameter();
      compute_cell_volume();
      compute_jmt();
    }

  submesh<cell_type, cell_type> query_cells(const std::function<bool(const double*)>& f) const {
    std::vector<bool> selected_cells(this->get_cell_number(), false);
    
    for (std::size_t k(0); k < this->get_cell_number(); ++k) {
      for (std::size_t n(0); n < vertices.get_size(1); ++n)
	selected_cells.at(k) = selected_cells.at(k) || f(&mesh<cell>::vertices.at(mesh<cell>::cells.at(k, n), 0));
    }
    
    return submesh_from_selection(selected_cells);
  }
  
  double get_cell_volume(std::size_t k) const {
    return cell_volume.at(k);
  }

  const array<double>& get_jmt(std::size_t k) const {
    return jmt[k];
  }

  submesh<cell_type, cell_type> get_submesh_with_reference(std::size_t ref_id) {
    const std::size_t cell_number(std::count(&references.at(0),
                                             &references.at(0) + references.get_size(0),
                                             ref_id));
    array<unsigned int> cells{cell_number, cell_type::n_vertex_per_cell};
    array<unsigned int> parent_cell_id{cell_number};
    array<unsigned int> parent_cell_subdomain_id{cell_number};
    parent_cell_subdomain_id.fill(0);

    
    for (std::size_t k(0), pk(0); pk < this->get_cell_number(); ++pk) {
      if (references.at(pk) == ref_id) {
	for (std::size_t n(0); n < cell_type::n_vertex_per_cell; ++n)
	  cells.at(k, n) = this->cells.at(pk, n);
	  
	parent_cell_id.at(k) = pk;
	
	k += 1;
      }
    }

    return submesh<cell_type, cell_type>(*this, cells, parent_cell_id, parent_cell_subdomain_id);
  }
  
  submesh<cell_type, ::cell::point> get_point_submesh(std::size_t k = 0) const {
    array<unsigned int> el_id{1}, sd_id{1}, el{1, 1};
    el_id.at(0) = k;
    sd_id.at(0) = 0;
    el.at(0, 0) = mesh<cell>::cells.at(k, 0);
    return  submesh<cell_type, ::cell::point>(*this, el, el_id, sd_id);
  }

  struct subdomain_info {
    std::size_t parent_cell_id, parent_cell_subdomain_id;
    ::cell::subdomain_type subdomain;

    bool operator<(const subdomain_info& op) const { return subdomain < op.subdomain; }
  };

  submesh<cell_type> get_boundary_submesh() const {
    using ::cell::subdomain_type;

    typedef std::map<subdomain_info, unsigned int> map_type;
    map_type subdomain_count;

    const std::size_t subdomain_id(cell_type::n_subdomain_type - 2);
    for (unsigned int k(0); k < this->get_cell_number(); ++k) {
      for (unsigned int j(0); j < cell_type::n_subdomain(subdomain_id); ++j) {	
	subdomain_info current{k, j, cell::get_subdomain(mesh<cell>::cells, k, subdomain_id, j)};
	subdomain_count[current] += 1;
      }
    }
    const std::size_t n_cells(std::count_if(subdomain_count.begin(),
                                            subdomain_count.end(),
                                            [](const typename map_type::value_type& subdomain) {
                                              return subdomain.second == 1;
                                            }));

    array<unsigned int> el_id{n_cells};
    array<unsigned int> sd_id{n_cells};
    array<unsigned int> el{n_cells, cell_type::boundary_cell_type::n_vertex_per_cell};

    unsigned int n(0);
    for (const auto& info: subdomain_count) {
      if (info.second == 1) {
	el_id.at(n) = info.first.parent_cell_id;
	sd_id.at(n) = info.first.parent_cell_subdomain_id;
	std::copy(info.first.subdomain.begin(),
		  info.first.subdomain.end(),
		  &el.at(n, 0));
	++n;
      }
    }
    
    return submesh<cell_type>(*this, el, el_id, sd_id);
  }

  void show(std::ostream& stream) const {
    for (std::size_t k(0); k < this->get_cell_number(); ++k) {
      stream << "cell " << k << ": ";
      for (std::size_t i(0); i < mesh<cell>::cells.get_size(1); ++i) {
	stream << mesh<cell>::cells.at(k, i) << " ";
      }
          stream << std::endl;
    }
  }

#warning "FIXME: implement a better than O(N) algorithm."
    std::size_t get_cell_at(const double* x) const {
      for (std::size_t k(0); k < this->get_cell_number(); ++k) {
        const array<double>
          bc_coord(cell_type::get_barycentric_coordinates(mesh<cell>::vertices,
                                                          mesh<cell>::cells, k, x));

        const double* min(std::min_element(&bc_coord.at(0),
                                           &bc_coord.at(0) + cell_type::n_vertex_per_cell));
        if (*min >= 0.0)
          return k;
      }

      throw std::string("point out of bounds");
    }

  double get_h_max() const {
    return h_max;
  }

  double get_cell_diameter(std::size_t k) const {
    return h.at(k);
  }
  
private:
  array<double> cell_volume;
  array<double> h;
  double h_max;
  std::vector<array<double> > jmt;
  
  void compute_cell_diameter() {
    for (std::size_t k(0); k < this->get_cell_number(); ++k)
      h.at(k) = cell_type::cell_diameter(mesh<cell>::vertices,
                                         mesh<cell>::cells, k);

    if (this->get_cell_number())
      h_max = *std::max_element(&h.at(0),
                                &h.at(0) + this->get_cell_number());
    else
      h_max = 0.0;
  }

  void compute_cell_volume() {
    for (std::size_t k(0); k < this->get_cell_number(); ++k)
      cell_volume.at(k) = cell_type::get_cell_volume(mesh<cell>::vertices,
                                                     mesh<cell>::cells, k);
  }

  void compute_jmt() {
    jmt.reserve(this->get_cell_number());
    for (std::size_t k(0); k < this->get_cell_number(); ++k)
      jmt.push_back(cell_type::get_jmt(mesh<cell>::vertices,
                                       mesh<cell>::cells, k));
  }

  submesh<cell_type, cell_type> submesh_from_selection(
                                                       const std::vector<bool>& selected_cells) const {
    
    const std::size_t n_cells(std::count(selected_cells.begin(),
                                         selected_cells.end(),
                                         true));

    array<unsigned int> el_id{n_cells};
    array<unsigned int> sd_id{n_cells};
    array<unsigned int> el{n_cells, cell_type::n_vertex_per_cell};

    for (std::size_t n(0), m(0); n < this->get_cell_number(); ++n) {
      if (selected_cells[n]) {
	el_id.at(m) = n;
	sd_id.at(m) = 0;
	std::copy(&mesh<cell>::cells.at(n, 0),
		  &mesh<cell>::cells.at(n, 0) + cell_type::n_vertex_per_cell,
		  &el.at(m, 0));
	++m;
      }
    }
    
    return submesh<cell_type, cell_type>(*this, el, el_id, sd_id);
  }
};


fe_mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n);
fe_mesh<cell::triangle> gen_square_mesh(double x_1, double x_2,
                                        unsigned int n_1, unsigned int n_2);
fe_mesh<cell::tetrahedron> gen_cube_mesh(double x_1, double x_2, double x_3,
                                         unsigned int n_1, unsigned int n_2, unsigned int n_3);

#endif /* _MESH_H_ */
