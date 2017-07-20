#ifndef _CELL_H_
#define _CELL_H_

#include <set>
#include <cmath>

#include <spikes/array.hpp>

#include "linear_algebra.hpp"


namespace cell {
  namespace detail {
    template<std::size_t n>
    struct subdomain {
      subdomain(): fill(0) {}

      void insert(unsigned int node) {
	if (fill and nodes[fill - 1] >= node)
	  throw std::string("nodes must be ordered when inserted into a subdomain");
	
	nodes[fill] = node;
	fill += 1;
      }

      std::size_t size() const { return fill; }

      unsigned int* begin() { return nodes; }
      unsigned int* end() { return nodes + fill; }

      const unsigned int* begin() const { return nodes; }
      const unsigned int* end() const { return nodes + fill; }

      bool operator==(const subdomain<n>& s) const {
	if (fill != s.fill)
	  return false;
	
	for (unsigned int i(0); i < fill; ++i)
	  if (nodes[i] != s.nodes[i])
	    return false;
      }

      bool operator<(const subdomain<n>& s) const {
	if (fill != s.fill)
	  return fill < s.fill;
	for (unsigned int i(0); i < fill; ++i)
	  if (nodes[i] != s.nodes[i])
	    return nodes[i] < s.nodes[i];
	return false;
      }
      
      unsigned int nodes[n];
      unsigned short fill;
    };
  }
  
  //typedef std::set<std::size_t> subdomain_type;
  typedef detail::subdomain<4> subdomain_type;

  struct empty_cell {};
  
  struct point {
    static const std::size_t n_dimension = 0;
    static const std::size_t n_vertex_per_element = 1;
    static const std::size_t n_subdomain_type = 1;
    static const std::size_t n_subdomain_of_type[1];
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

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   std::size_t k) {
      return 0.0;
    }
  };
  
  struct edge {
    static const std::size_t n_dimension = 1;
    static const std::size_t n_vertex_per_element = 2;
    static const std::size_t n_subdomain_type = 2;
    static const std::size_t n_subdomain_of_type[2];
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

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   std::size_t k) {
      return std::abs(vertices.at(elements.at(k, 0), 0)
		      - vertices.at(elements.at(k, 1), 0));
    }

  };
    
  struct triangle {
    static const std::size_t n_dimension = 2;
    static const std::size_t n_vertex_per_element = 3;
    static const std::size_t n_subdomain_type = 3;
    static const std::size_t n_subdomain_of_type[3];
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
      
      return 0.5 * std::abs(a_1 * b_2 - a_2 * b_1);
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
      for (unsigned int i(0); i < 2 - 1; ++i)
	for (unsigned int j(i + 1); j < 2; ++j) {
	  double side_length(0.0);
	  for (unsigned int n(0); n < 2; ++n)
	    side_length += std::pow(  vertices.at(elements.at(k, i), n)
				    - vertices.at(elements.at(k, j), n), 2);
	  side_length = std::sqrt(side_length);

	  longest_side = std::max(longest_side, side_length);
	}
      return longest_side;
    }
  };

}

#endif /* _CELL_H_ */
