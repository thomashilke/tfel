#ifndef _CELL_H_
#define _CELL_H_

#include <set>

#include "mesh.hpp"
#include "meta.hpp"

namespace cell {
  typedef std::set<std::size_t> subdomain_type;
  
  struct point { };
  
  struct edge {
    static const std::size_t n_dimension = 1;
    static const std::size_t n_vertex_per_element = 2;
    static const std::size_t n_subdomain_type = 2;
    static const bool is_simplicial = true;

    static const std::size_t n_subdomain(unsigned int i) {
      static const std::size_t n_sub[] = {2, 1};
      return n_sub[i];
    }

    static const subdomain_type get_subdomain(const mesh<edge>& m,
					      std::size_t k,
					      std::size_t sd, std::size_t j) {
      switch (sd) {
      case 0: return get_vertex(m, k, j);
      case 1: return get_edge(m, k, j);
      }
      throw std::string("invalid subdomain id");
    }

    static subdomain_type get_vertex(const mesh<edge>& m,
				     std::size_t k,
				     std::size_t j) {
      subdomain_type vertex;
      vertex.insert(m.get_elements().at(k, j));
      return vertex;
    }
    
    static subdomain_type get_edge(const mesh<edge>& m,
				   std::size_t k,
				   std::size_t j) {
      subdomain_type element;
      for (unsigned int i(0); i < 2; ++i)
	element.insert(m.get_elements().at(k, i));
      return element;
    }


    static std::set<subdomain_type> get_subdomain_list(const mesh<edge>& m,
						       std::size_t sd) {
      switch (sd) {
      case 0: return get_vertex_list(m);
      case 1: return get_edge_list(m);
      }
      throw std::string("invalid subdomain id");
    }

    static std::set<subdomain_type> get_vertex_list(const mesh<edge>& m) {
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

    static std::set<subdomain_type> get_edge_list(const mesh<edge>& m) {
      std::set<subdomain_type> element_list;

      const array<unsigned int>& elements(m.get_elements());
      for (unsigned int k(0); k < elements.get_size(0); ++k) {
	subdomain_type edge;
	for (unsigned int i(0); i < elements.get_size(1); ++i)
	  edge.insert(elements.at(k, i));
	element_list.insert(edge);
      }
      
      return element_list;
    }
  };
    
  struct triangle {
    static const std::size_t n_dimension = 2;
    static const std::size_t n_vertex_per_element = 3;
    static const std::size_t n_subdomain_type = 3;
    static const bool is_simplicial = true;

    static const std::size_t n_subdomain(unsigned int i) {
      static const std::size_t n_sub[] = {3, 3, 1};
      return n_sub[i];
    }

    static subdomain_type get_subdomain(const mesh<triangle>& m,
					std::size_t k,
					std::size_t sd, std::size_t j) {
      switch (sd) {
      case 0: return get_vertex(m, k, j);
      case 1: return get_edge(m, k, j);
      case 2: return get_element(m, k, j);
      }
      throw std::string("invalid subdomain id");
    }
    
    static subdomain_type get_vertex(const mesh<triangle>& m, std::size_t k, std::size_t j) {
      subdomain_type vertex;
      vertex.insert(m.get_elements().at(k, j));
      return vertex;
    }
    
    static subdomain_type get_edge(const mesh<triangle>& m, std::size_t k, std::size_t j) {
      subdomain_type edge;
      if (j < 2) {
	edge.insert(m.get_elements().at(k, 0));
	edge.insert(m.get_elements().at(k, j + 1));
      } else {
	edge.insert(m.get_elements().at(k, 1));
	edge.insert(m.get_elements().at(k, 2));
      }
      return edge;
    }

    static subdomain_type get_element(const mesh<triangle>& m, std::size_t k, std::size_t j) {
      subdomain_type element;
      for (unsigned int i(0); i < 3; ++i)
	element.insert(m.get_elements().at(k, i));
      return element;
    }

    static std::set<subdomain_type> get_subdomain_list(const mesh<triangle>& m,
						       std::size_t sd) {
      switch (sd) {
      case 0: return get_vertex_list(m);
      case 1: return get_edge_list(m);
      case 2: return get_element_list(m);
      }
      throw std::string("invalid subdomain id");
    }
    
    static std::set<subdomain_type> get_vertex_list(const mesh<triangle>& m) {
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

    static std::set<subdomain_type> get_edge_list(const mesh<triangle>& m) {
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

    static std::set<subdomain_type> get_element_list(const mesh<triangle>& m) {
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

}

#endif /* _CELL_H_ */
