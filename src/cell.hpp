#ifndef _CELL_H_
#define _CELL_H_

#include <set>

#include "mesh.hpp"
#include "meta.hpp"

namespace cell {
  typedef std::set<std::size_t> subdomain_type;
  
  struct point {
    static const unsigned int dimension = 0;
    static const unsigned int node_number = 1;
    static const bool is_simplicial = true;

    /*typedef empty_type_list cell_type_of_submanifold;
      typedef empty_type_list disjoint_submanifolds;*/
  };
  
  struct edge {
    static const unsigned int n_dimension = 1;
    static const unsigned int node_number = 2;
    static const bool is_simplicial = true;

    /*typedef type_list<point> cell_type_of_submanifold;
      typedef type_list<holder<unsigned int, 2> > disjoint_submanifolds;*/
  };
    
  struct triangle {
    static const unsigned int n_dimension = 2;
    static const unsigned int n_vertex_per_element = 3;
    static const bool is_simplicial = true;
    
    /*typedef type_list<point, edge> cell_type_of_submanifold;
    typedef type_list<holder<unsigned int, 3>,
    holder<unsigned int, 3> > disjoint_submanifolds;*/

    std::set<subdomain_type> get_vertex_list(const mesh<triangle>& m) {
      std::set<subdomain_type> vertex_list;

      const array<unsigned int>& elements(m.get_elements());
      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	for (unsigned int i(0); i < elements.get_size(1); ++i) {
	  subdomain_type vertex;
	  vertex.insert(elements.at(k, i));
	  vertex_list.insert(vertex);
	}
      }
      
      return vertex_list;
    }

    std::set<subdomain_type> get_edge_list(const mesh<triangle>& m) {
      std::set<subdomain_type> edge_list;

      const array<unsigned int>& elements(m.get_elements());
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

    std::set<subdomain_type> get_element_list(const mesh<triangle>& m) {
      std::set<subdomain_type> element_list;

      const array<unsigned int>& elements(m.get_elements());
      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	subdomain_type triangle;
	for (unsigned int i(0); i < elements.get_size(1); ++i)
	  triangle.insert(elements.at(k, i));
	element_list.insert(triangle);
      }
      
      return element_list;
    }
    
  };

  //typedef type_list<point, edge, triangle> simplex_list;
  
  /*template<unsigned int dimension>
  struct simplex_of_dim {
    typedef 
    typename at<dimension, simplex_list>::type
    type;
    };*/
}

#endif /* _CELL_H_ */
