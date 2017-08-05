#ifndef _CELL_H_
#define _CELL_H_

#include <set>
#include <cmath>
#include <vector>

#include <spikes/array.hpp>

#include "linear_algebra.hpp"
#include "vector_operation.hpp"
#include "subdomain.hpp"

namespace cell {

  struct empty_cell {};


  struct point {
    static const std::size_t n_dimension = 0;
    static const std::size_t n_vertex_per_element = 1;
    static const std::size_t n_subdomain_type = 1;
    static constexpr std::size_t n_subdomain_of_type[1] = {1};
    static const bool is_simplicial = true;

    typedef empty_cell boundary_cell_type;
    
    static std::set<subdomain_type> get_vertex_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> vertex_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1); ++i) {
	  subdomain_type vertex;
	  vertex.insert(elements.at(k, i));
	  vertex_list.insert(vertex);
	}
      }
      
      return vertex_list;      
    }

    static std::set<subdomain_type> get_subdomain_list(const array<unsigned int>& elements,
						       std::size_t sd) {
      switch (sd) {
      case 0: return get_vertex_list(elements);
      }
      throw std::string("invalid subdomain id");
    }

    static array<double> map_points_to_subdomain(std::size_t subdomain_id, const array<double>& xs) {
      if (subdomain_id >= 1)
	throw std::string("point::map_points_on_subdomain: only one subdomains for a point.");
      
      return xs;
    }
    
    static double get_cell_volume(const array<double>& /*vertices*/,
				  const array<unsigned int>& /*elements*/,
				  unsigned int /*k*/) {
      return 1.0;
    }

    static array<double> map_points_to_space_coordinates(const array<double>& vertices,
							 const array<unsigned int>& elements,
							 std::size_t subdomain_id, const array<double>& xs) {
      array<double> hat_xs(xs);
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n)
	  hat_xs.at(i, n) = vertices.at(elements.at(subdomain_id, 0), n);
      }
      
      return hat_xs;
    }


    static void map_points_to_space_coordinates(array<double>& hat_xs,
						const array<double>& vertices,
						const array<unsigned int>& elements,
						std::size_t subdomain_id, const array<double>& xs) {
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n)
	  hat_xs.at(i, n) = vertices.at(elements.at(subdomain_id, 0), n);
      }
    }

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   std::size_t k) {
      return 0.0;
    }

    static array<double> barycenter(const array<double>& vertices,
				    const array<unsigned int>& elements,
				    std::size_t k) {
      array<double> bc{vertices.get_size(1)};
      std::copy(&vertices.at(elements.at(k, 0), 0),
		&vertices.at(elements.at(k, 0), 0) + vertices.get_size(1),
		&bc.at(0));
      return bc;
    }
  };

  
  struct edge {
    static const std::size_t n_dimension = 1;
    static const std::size_t n_vertex_per_element = 2;
    static const std::size_t n_subdomain_type = 2;
    static constexpr std::size_t n_subdomain_of_type[2] = {2, 1};
    static const bool is_simplicial = true;

    typedef point boundary_cell_type;
    
    static std::size_t n_subdomain(unsigned int i) {
      static const std::size_t n_sub[] = {2, 1};
      return n_sub[i];
    }

    static const subdomain_type get_subdomain(const array<unsigned int>& elements,
					      std::size_t k,
					      std::size_t sd, std::size_t j) {
      switch (sd) {
      case 0: return get_vertex(elements, k, j);
      case 1: return get_edge(elements, k, j);
      }
      throw std::string("invalid subdomain id");
    }

    static subdomain_type get_vertex(const array<unsigned int>& elements,
				     std::size_t k,
				     std::size_t j) {
      subdomain_type vertex;
      vertex.insert(elements.at(k, j));
      return vertex;
    }
    
    static subdomain_type get_edge(const array<unsigned int>& elements,
				   std::size_t k,
				   std::size_t j) {
      if (j != 0)
	throw std::string("invalid edge id");
      
      subdomain_type element;
      for (unsigned int i(0); i < 2; ++i)
	element.insert(elements.at(k, i));
      return element;
    }

    static std::set<subdomain_type> get_subdomain_list(const array<unsigned int>& elements,
						       std::size_t sd) {
      switch (sd) {
      case 0: return get_vertex_list(elements);
      case 1: return get_edge_list(elements);
      }
      throw std::string("invalid subdomain id");
    }

    static std::set<subdomain_type> get_vertex_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> vertex_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1); ++i) {
	  subdomain_type vertex;
	  vertex.insert(elements.at(k, i));
	  vertex_list.insert(vertex);
	}
      }
      
      return vertex_list;
    }

    static std::set<subdomain_type> get_edge_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> element_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	subdomain_type edge;
	for (unsigned int i(0); i < elements.get_size(1); ++i)
	  edge.insert(elements.at(k, i));
	element_list.insert(edge);
      }
      
      return element_list;
    }

    static double get_cell_volume(const array<double>& vertices,
				  const array<unsigned int>& elements,
				  unsigned int k) {
      return std::abs(vertices.at(elements.at(k, 0), 0)
		      - vertices.at(elements.at(k, 1), 0));
    }

    static array<double> get_jmt(const array<double>& vertices,
				 const array<unsigned int>& elements,
				 unsigned int k) {
      array<double> jmt{1,1};
      jmt.at(0,0) = 1.0 / (vertices.at(elements.at(k, 1), 0)
			   - vertices.at(elements.at(k, 0), 0));
      return jmt;
    }
    
    static array<double> map_points_to_subdomain(std::size_t subdomain_id, const array<double>& xs) {
      if (subdomain_id >= 2)
	throw std::string("edge::map_points_on_subdomain: only two subdomains for an edge.");
      
      array<double> hat_xs(xs);
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	if (subdomain_id == 0)
	  hat_xs.at(i, 0) = 0.0;
	else
	  hat_xs.at(i, 0) = 1.0;
      }

      return hat_xs;
    }

    static array<double> map_points_to_space_coordinates(const array<double>& vertices,
							 const array<unsigned int>& elements,
							 std::size_t subdomain_id, const array<double>& xs) {
      array<double> hat_xs{xs.get_size(0), vertices.get_size(1)};
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n)
	  hat_xs.at(i, n) = vertices.at(elements.at(subdomain_id, 0), n) +
	    xs.at(i, 0) * (vertices.at(elements.at(subdomain_id, 1), n)
			   - vertices.at(elements.at(subdomain_id, 0), n));
      }
      
      return hat_xs;
    }

    static void map_points_to_space_coordinates(array<double>& hat_xs,
						const array<double>& vertices,
						const array<unsigned int>& elements,
						std::size_t subdomain_id, const array<double>& xs) {
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n)
	  hat_xs.at(i, n) = vertices.at(elements.at(subdomain_id, 0), n) +
	    xs.at(i, 0) * (vertices.at(elements.at(subdomain_id, 1), n)
			   - vertices.at(elements.at(subdomain_id, 0), n));
      }
    }

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   std::size_t k) {
      return std::abs(vertices.at(elements.at(k, 0), 0)
		      - vertices.at(elements.at(k, 1), 0));
    }

    static array<double> normal(const array<double>& vertices,
				const array<unsigned int>& elements,
				std::size_t k) {
      if (vertices.get_size(1) != 2 or elements.get_size(1) != 2)
	throw std::string("cell::edge: normal is only defined for edges embedded in 2d space");

      array<double> n{2};
      n.at(0) =  - (vertices.at(elements.at(k, 1), 1) - vertices.at(elements.at(k, 0), 1));
      n.at(1) =    (vertices.at(elements.at(k, 1), 0) - vertices.at(elements.at(k, 0), 0));

      return n;
    }

    static array<double> get_barycentric_coordinate_map(const array<double>& vertices,
							const array<unsigned int>& elements,
							std::size_t k) {
      array<double> map{n_vertex_per_element, n_vertex_per_element};

      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	map.at(0, j) = 1.0;
      
      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	for (std::size_t i(0); i < n_dimension; ++i)
	  map.at(i + 1, j) = vertices.at(elements.at(k, j), i);

      matrix_inverse(&map.at(0, 0), n_vertex_per_element);
      
      return map;
    }

    static array<double> get_barycentric_coordinates(const array<double>& vertices,
						     const array<unsigned int>& elements,
						     std::size_t k,
						     const double* x) {
      const array<double> bc_map(get_barycentric_coordinate_map(vertices, elements, k));
      array<double>
	x_coord{n_vertex_per_element},
	bc_coord{n_vertex_per_element};

      x_coord.at(0) = 1.0;
      for (std::size_t i(0); i < n_dimension; ++i)
	x_coord.at(i + 1) = x[i];

      bc_coord.fill(0.0);
      for (std::size_t i(0); i < n_vertex_per_element; ++i)
	for (std::size_t j(0); j < n_vertex_per_element; ++j)
	  bc_coord.at(i) += bc_map.at(i, j) * x_coord.at(j);

      return bc_coord;
    }


    static array<double> barycenter(const array<double>& vertices,
				    const array<unsigned int>& elements,
				    std::size_t k) {
      array<double> bc{vertices.get_size(1)};
      const double weight(1.0 / elements.get_size(1));
      for (std::size_t m(0); m < vertices.get_size(1); ++m) {
	bc.at(m) = 0.0;
	for (std::size_t n(0); n < elements.get_size(1); ++n)
	  bc.at(m) += vertices.at(elements.at(k, n), m) * weight;
      }

      return bc;
    }

    struct fe {
      struct lagrange_p0;
      struct lagrange_p1;
      struct lagrange_p1_bubble;
    };
  };


  struct triangle {
    static const std::size_t n_dimension = 2;
    static const std::size_t n_vertex_per_element = 3;
    static const std::size_t n_subdomain_type = 3;
    static constexpr std::size_t n_subdomain_of_type[3] = {3, 3, 1};
    static const bool is_simplicial = true;

    typedef edge boundary_cell_type;

    static std::size_t n_subdomain(unsigned int i) {
      static const std::size_t n_sub[] = {3, 3, 1};
      return n_sub[i];
    }

    static subdomain_type get_subdomain(const array<unsigned int>& elements,
					std::size_t k,
					std::size_t sd, std::size_t j) {
      switch (sd) {
      case 0: return get_vertex(elements, k, j);
      case 1: return get_edge(elements, k, j);
      case 2: return get_element(elements, k, j);
      }
      throw std::string("invalid subdomain id");
    }
    
    static subdomain_type get_vertex(const array<unsigned int>& elements,
				     std::size_t k, std::size_t j) {
      subdomain_type vertex;
      vertex.insert(elements.at(k, j));
      return vertex;
    }
    
    static subdomain_type get_edge(const array<unsigned int>& elements,
				   std::size_t k, std::size_t j) {
      subdomain_type edge;
      if (j < 2) {
	edge.insert(elements.at(k, 0));
	edge.insert(elements.at(k, j + 1));
      } else {
	edge.insert(elements.at(k, 1));
	edge.insert(elements.at(k, 2));
      }
      return edge;
    }

    static subdomain_type get_element(const array<unsigned int>& elements,
				      std::size_t k, std::size_t j) {
      subdomain_type element;
      for (unsigned int i(0); i < 3; ++i)
	element.insert(elements.at(k, i));
      return element;
    }

    static std::set<subdomain_type> get_subdomain_list(const array<unsigned int>& elements,
						       std::size_t sd) {
      switch (sd) {
      case 0: return get_vertex_list(elements);
      case 1: return get_edge_list(elements);
      case 2: return get_element_list(elements);
      }
      throw std::string("invalid subdomain id");
    }
    
    static std::set<subdomain_type> get_vertex_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> vertex_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1); ++i) {
	  subdomain_type vertex;
	  vertex.insert(elements.at(k, i));
	  vertex_list.insert(vertex);
	}
      }
      
      return vertex_list;
    }

    static std::set<subdomain_type> get_edge_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> edge_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1) - 1; ++i) {
	  for (unsigned int j(i + 1); j < elements.get_size(1); ++j) {
	    subdomain_type edge;
	    edge.insert(elements.at(k, i));
	    edge.insert(elements.at(k, j));
	    edge_list.insert(edge);
	  }
	}
      }
      
      return edge_list;
    }

    static std::set<subdomain_type> get_element_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> element_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	subdomain_type triangle;
	for (unsigned int i(0); i < elements.get_size(1); ++i)
	  triangle.insert(elements.at(k, i));
	element_list.insert(triangle);
      }
      
      return element_list;
    }

    /*
     * Maps points defined on the reference triangle to space coordinates
     */
    static array<double> map_points_to_space_coordinates(const array<double>& vertices,
							 const array<unsigned int>& elements,
							 std::size_t subdomain_id, const array<double>& xs) {
      array<double> hat_xs{xs.get_size(0), vertices.get_size(1)};
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n)
	  hat_xs.at(i, n) = vertices.at(elements.at(subdomain_id, 0), n) 
	    + xs.at(i, 0) * (vertices.at(elements.at(subdomain_id, 1), n)
			     - vertices.at(elements.at(subdomain_id, 0), n))
	    + xs.at(i, 1) * (vertices.at(elements.at(subdomain_id, 2), n)
			   - vertices.at(elements.at(subdomain_id, 0), n));
      }
      
      return hat_xs;
    }

    static void map_points_to_space_coordinates(array<double>& hat_xs,
						const array<double>& vertices,
						const array<unsigned int>& elements,
						std::size_t subdomain_id, const array<double>& xs) {
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n)
	  hat_xs.at(i, n) = vertices.at(elements.at(subdomain_id, 0), n) 
	    + xs.at(i, 0) * (vertices.at(elements.at(subdomain_id, 1), n)
			     - vertices.at(elements.at(subdomain_id, 0), n))
	    + xs.at(i, 1) * (vertices.at(elements.at(subdomain_id, 2), n)
			   - vertices.at(elements.at(subdomain_id, 0), n));
      }
    }


    /*
     * Return the integration-wise cell volume
     */
    static double get_cell_volume(const array<double>& vertices,
				  const array<unsigned int>& elements,
				  unsigned int k) {
      const double
	a_1(vertices.at(elements.at(k, 1), 0) - vertices.at(elements.at(k, 0), 0)),
	a_2(vertices.at(elements.at(k, 1), 1) - vertices.at(elements.at(k, 0), 1)),
	b_1(vertices.at(elements.at(k, 2), 0) - vertices.at(elements.at(k, 0), 0)),
	b_2(vertices.at(elements.at(k, 2), 1) - vertices.at(elements.at(k, 0), 1));
      
      return 1.0 / 2.0 * std::abs(a_1 * b_2 - a_2 * b_1);
    }

    /*
     * Return the inverse transpose of the jacobian of the T_k(x) mapping
     */
    static array<double> get_jmt(const array<double>& vertices,
				 const array<unsigned int>& elements,
				 unsigned int k) {
      array<double> jmt{2,2};
      jmt.at(0,0) = vertices.at(elements.at(k, 1), 0) - vertices.at(elements.at(k, 0), 0);
      jmt.at(0,1) = vertices.at(elements.at(k, 2), 0) - vertices.at(elements.at(k, 0), 0);
      jmt.at(1,0) = vertices.at(elements.at(k, 1), 1) - vertices.at(elements.at(k, 0), 1);
      jmt.at(1,1) = vertices.at(elements.at(k, 2), 1) - vertices.at(elements.at(k, 0), 1);

      matrix_inverse(&jmt.at(0,0), 2);
      std::swap(jmt.at(0, 1), jmt.at(1, 0));
      
      return jmt;
    }

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   std::size_t k) {
      double longest_side(0.0);
      for (unsigned int i(0); i < 3 - 1; ++i)
	for (unsigned int j(i + 1); j < 3; ++j) {
	  double side_length(0.0);
	  for (unsigned int n(0); n < 2; ++n)
	    side_length += std::pow(  vertices.at(elements.at(k, i), n)
				    - vertices.at(elements.at(k, j), n), 2);
	  side_length = std::sqrt(side_length);

	  longest_side = std::max(longest_side, side_length);
	}
      return longest_side;
    }

    static array<double> get_barycentric_coordinate_map(const array<double>& vertices,
							const array<unsigned int>& elements,
							std::size_t k) {
      array<double> map{n_vertex_per_element, n_vertex_per_element};

      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	map.at(0, j) = 1.0;
      
      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	for (std::size_t i(0); i < n_dimension; ++i)
	  map.at(i + 1, j) = vertices.at(elements.at(k, j), i);

      matrix_inverse(&map.at(0, 0), n_vertex_per_element);
      
      return map;
    }

    static array<double> get_barycentric_coordinates(const array<double>& vertices,
						     const array<unsigned int>& elements,
						     std::size_t k,
						     const double* x) {
      const array<double> bc_map(get_barycentric_coordinate_map(vertices, elements, k));
      array<double>
	x_coord{n_vertex_per_element},
	bc_coord{n_vertex_per_element};

      x_coord.at(0) = 1.0;
      for (std::size_t i(0); i < n_dimension; ++i)
	x_coord.at(i + 1) = x[i];

      bc_coord.fill(0.0);
      for (std::size_t i(0); i < n_vertex_per_element; ++i)
	for (std::size_t j(0); j < n_vertex_per_element; ++j)
	  bc_coord.at(i) += bc_map.at(i, j) * x_coord.at(j);

      return bc_coord;
    }

    static array<double> barycenter(const array<double>& vertices,
				    const array<unsigned int>& elements,
				    std::size_t k) {
      array<double> bc{vertices.get_size(1)};
      const double weight(1.0 / elements.get_size(1));
      for (std::size_t m(0); m < vertices.get_size(1); ++m) {
	bc.at(m) = 0.0;
	for (std::size_t n(0); n < elements.get_size(1); ++n)
	  bc.at(m) += vertices.at(elements.at(k, n), m) * weight;
      }

      return bc;
    }

    static array<double> subdomain_normal(const array<double>& vertices,
					  const array<unsigned int>& elements,
					  std::size_t k,
					  std::size_t subdomain_id) {
      auto edge(get_edge(elements, k, subdomain_id));
      std::vector<unsigned int> vertices_id(edge.begin(), edge.end());
      std::size_t last_vertex(2 - subdomain_id);

      array<double> n{2};
      n.at(0) =  - (vertices.at(vertices_id[1], 1) - vertices.at(vertices_id[0], 1));
      n.at(1) =    (vertices.at(vertices_id[1], 0) - vertices.at(vertices_id[0], 0));

      array<double> m{2};
      m.at(0) =    (vertices.at(elements.at(k, last_vertex), 0) - vertices.at(vertices_id[0], 0));
      m.at(1) =    (vertices.at(elements.at(k, last_vertex), 1) - vertices.at(vertices_id[0], 1));

      if (dotp(n, m) > 0.0)
	for (std::size_t k(0); k < 2; ++k)
	  n.at(k) *= -1.0;

      return n;
    }

    struct fe {
      struct lagrange_p0;
      struct lagrange_p1;
      struct lagrange_p1_bubble;
    };
  };


  struct tetrahedron {
    static const std::size_t n_dimension = 3;
    static const std::size_t n_vertex_per_element = 4;
    static const std::size_t n_subdomain_type = 4;
    static constexpr const std::size_t n_subdomain_of_type[4] = {4, 6, 4, 1};
    static const bool is_simplicial = true;

    typedef triangle boundary_cell_type;

    static std::size_t n_subdomain(unsigned int i) {
      static const std::size_t n_sub[] = {4, 6, 4, 1};
      return n_sub[i];
    }

    static subdomain_type get_subdomain(const array<unsigned int>& elements,
					std::size_t k,
					std::size_t sd, std::size_t j) {
      switch (sd) {
      case 0: return get_vertex(elements, k, j);
      case 1: return get_edge(elements, k, j);
      case 2: return get_triangle(elements, k, j);
      case 3: return get_element(elements, k, j);
      }
      throw std::string("invalid subdomain id");
    }

    static subdomain_type get_vertex(const array<unsigned int>& elements,
				     std::size_t k,
				     std::size_t j) {
      subdomain_type vertex;
      vertex.insert(elements.at(k, j));
      return vertex;
    }
    
    static subdomain_type get_edge(const array<unsigned int>& elements,
				   std::size_t k,
				   std::size_t j) {
      subdomain_type edge;

      if (j < 3) {
	edge.insert(elements.at(k, 0));
	edge.insert(elements.at(k, 1 + j));
      } else if (j < 5) {
	edge.insert(elements.at(k, 1));
	edge.insert(elements.at(k, j - 1));
      } else {
	edge.insert(elements.at(k, 2));
	edge.insert(elements.at(k, 3));
      }
      
      return edge;
    }
    
    static subdomain_type get_triangle(const array<unsigned int>& elements,
				       std::size_t k,
				       std::size_t j) {
      subdomain_type triangle;

      if (j < 2) {
	triangle.insert(elements.at(k, 0));
	triangle.insert(elements.at(k, 1));
	triangle.insert(elements.at(k, 2 + j));
      } else {
	triangle.insert(elements.at(k, j - 2));
	triangle.insert(elements.at(k, 2));
	triangle.insert(elements.at(k, 3));
      }

      return triangle;
    }
    
    static subdomain_type get_element(const array<unsigned int>& elements,
				      std::size_t k,
				      std::size_t j) {
      subdomain_type element;
      for (std::size_t i(0); i < 4; ++i)
	element.insert(elements.at(k, i));
      return element;
    }

    static std::set<subdomain_type> get_subdomain_list(const array<unsigned int>& elements,
					     std::size_t sd) {
      switch (sd) {
      case 0: return get_vertex_list(elements);
      case 1: return get_edge_list(elements);
      case 2: return get_triangle_list(elements);
      case 3: return get_element_list(elements);
      }
      throw std::string("invalid subdomain id");
    }
    
    static std::set<subdomain_type> get_vertex_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> vertex_list;

      for (std::size_t k(0); k < elements.get_size(0); ++k)
	for (std::size_t i(0); i < elements.get_size(1); ++i) {
	  subdomain_type vertex;
	  vertex.insert(elements.at(k, i));
	  vertex_list.insert(vertex);
	}

      return vertex_list;
    }
    
    static std::set<subdomain_type> get_edge_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> edge_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1) - 1; ++i) {
	  for (unsigned int j(i + 1); j < elements.get_size(1); ++j) {
	    subdomain_type edge;
	    edge.insert(elements.at(k, i));
	    edge.insert(elements.at(k, j));
	    edge_list.insert(edge);
	  }
	}
      }
      
      return edge_list;
    }

    static std::set<subdomain_type> get_triangle_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> triangle_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1) - 2; ++i) {
	  for (unsigned int j(i + 1); j < elements.get_size(1) - 1; ++j) {
	    for (unsigned int l(j + 1); l < elements.get_size(1); ++l) {
	      subdomain_type triangle;
	      triangle.insert(elements.at(k, i));
	      triangle.insert(elements.at(k, j));
	      triangle.insert(elements.at(k, l));
	      triangle_list.insert(triangle);
	    }
	  }
	}
      }
      
      return triangle_list;
    }

    static std::set<subdomain_type> get_element_list(const array<unsigned int>& elements) {
      std::set<subdomain_type> element_list;

      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	subdomain_type tetrahedron;
	for (unsigned int i(0); i < elements.get_size(1); ++i)
	  tetrahedron.insert(elements.at(k, i));
	element_list.insert(tetrahedron);
      }
      
      return element_list;
    }

    static array<double> map_points_to_space_coordinates(const array<double>& vertices,
							 const array<unsigned int>& elements,
							 std::size_t k,
							 const array<double>& xs) {
      array<double> hat_xs{xs.get_size(0), vertices.get_size(1)};

      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n) {
	  hat_xs.at(i, n) = vertices.at(elements.at(k, 0), n);
	  for (std::size_t m(0); m < n_vertex_per_element - 1; ++m)
	    hat_xs.at(i, n) += xs.at(i, m) * (vertices.at(elements.at(k, m + 1), n)
					      - vertices.at(elements.at(k, 0), n));
	}
      }
      
      return hat_xs;
    }

    static void map_points_to_space_coordinates(array<double>& hat_xs,
						const array<double>& vertices,
						const array<unsigned int>& elements,
						std::size_t k,
						const array<double>& xs) {
      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n) {
	  hat_xs.at(i, n) = vertices.at(elements.at(k, 0), n);
	  for (std::size_t m(0); m < n_vertex_per_element - 1; ++m)
	    hat_xs.at(i, n) += xs.at(i, m) * (vertices.at(elements.at(k, m + 1), n)
					      - vertices.at(elements.at(k, 0), n));
	}
      }
    }

    static double get_cell_volume(const array<double>& vertices,
				  const array<unsigned int>& elements,
				  unsigned int k) {
      array<double> a{3}, b{3}, c{3};
      a.at(0) = vertices.at(elements.at(k, 1), 0) - vertices.at(elements.at(k, 0), 0);
      a.at(1) = vertices.at(elements.at(k, 1), 1) - vertices.at(elements.at(k, 0), 1);
      a.at(2) = vertices.at(elements.at(k, 1), 2) - vertices.at(elements.at(k, 0), 2);

      b.at(0) = vertices.at(elements.at(k, 2), 0) - vertices.at(elements.at(k, 0), 0);
      b.at(1) = vertices.at(elements.at(k, 2), 1) - vertices.at(elements.at(k, 0), 1);
      b.at(2) = vertices.at(elements.at(k, 2), 2) - vertices.at(elements.at(k, 0), 2);

      c.at(0) = vertices.at(elements.at(k, 3), 0) - vertices.at(elements.at(k, 0), 0);
      c.at(1) = vertices.at(elements.at(k, 3), 1) - vertices.at(elements.at(k, 0), 1);
      c.at(2) = vertices.at(elements.at(k, 3), 2) - vertices.at(elements.at(k, 0), 2);

      return std::abs(1.0 / 6.0 * dotp(crossp(a, b), c));
    }

    static array<double> get_jmt(const array<double>& vertices,
				 const array<unsigned int>& elements,
				 unsigned int k) {
      assert(vertices.get_size(1) == n_dimension);
      
      array<double> jmt{n_dimension, n_dimension};

      for (std::size_t i(0); i < n_dimension; ++i)
	for (std::size_t j(0); j < n_dimension; ++j) {
	  jmt.at(j, i) = vertices.at(elements.at(k, 1 + i), j) - vertices.at(elements.at(k, 0), j);
	}

      matrix_inverse(&jmt.at(0,0), n_dimension);

      for (std::size_t i(0); i < n_dimension - 1; ++i)
	for (std::size_t j(i + 1); j < n_dimension; ++j)
	  std::swap(jmt.at(i, j), jmt.at(j, i));

      return jmt;
    }

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   unsigned int k) {
      double longest_side(0.0);
      for (unsigned int i(0); i < n_vertex_per_element - 1; ++i)
	for (unsigned int j(i + 1); j < n_vertex_per_element; ++j) {
	  double side_length(0.0);
	  for (unsigned int n(0); n < vertices.get_size(1); ++n)
	    side_length += std::pow(vertices.at(elements.at(k, i), n) -
				    vertices.at(elements.at(k, j), n),
				    2);
	  side_length = std::sqrt(side_length);

	  longest_side = std::max(longest_side, side_length);
	}
      return longest_side;
    }

    static array<double> get_barycentric_coordinate_map(const array<double>& vertices,
							const array<unsigned int>& elements,
							std::size_t k) {
      array<double> map{n_vertex_per_element, n_vertex_per_element};

      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	map.at(0, j) = 1.0;
      
      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	for (std::size_t i(0); i < n_dimension; ++i)
	  map.at(i + 1, j) = vertices.at(elements.at(k, j), i);

      matrix_inverse(&map.at(0, 0), n_vertex_per_element);
      
      return map;
    }

    static array<double> get_barycentric_coordinates(const array<double>& vertices,
						     const array<unsigned int>& elements,
						     std::size_t k,
						     const double* x) {
      const array<double> bc_map(get_barycentric_coordinate_map(vertices, elements, k));
      array<double>
	x_coord{n_vertex_per_element},
	bc_coord{n_vertex_per_element};

      x_coord.at(0) = 1.0;
      for (std::size_t i(0); i < n_dimension; ++i)
	x_coord.at(i + 1) = x[i];

      bc_coord.fill(0.0);
      for (std::size_t i(0); i < n_vertex_per_element; ++i)
	for (std::size_t j(0); j < n_vertex_per_element; ++j)
	  bc_coord.at(i) += bc_map.at(i, j) * x_coord.at(j);

      return bc_coord;
    }
    
    static array<double> barycenter(const array<double>& vertices,
				    const array<unsigned int>& elements,
				    std::size_t k) {
      array<double> bc{vertices.get_size(1)};
      const double weight(1.0 / elements.get_size(1));
      for (std::size_t m(0); m < vertices.get_size(1); ++m) {
	bc.at(m) = 0.0;
	for (std::size_t n(0); n < elements.get_size(1); ++n)
	  bc.at(m) += vertices.at(elements.at(k, n), m) * weight;
      }

      return bc;
    }

    struct fe {
      struct lagrange_p0;
      struct lagrange_p1;
    };
  };

}

#include "fe.hpp"

#endif /* _CELL_H_ */
