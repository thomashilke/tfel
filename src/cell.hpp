#ifndef _CELL_H_
#define _CELL_H_

#include "meta.hpp"

namespace cell {
  struct point {
    static const unsigned int dimension = 0;
    static const unsigned int node_number = 1;
    static const bool is_simplicial = true;

    typedef empty_type_list cell_type_of_submanifold;
    typedef empty_type_list disjoint_submanifolds;
  };
  
  struct edge {
    static const unsigned int dimension = 1;
    static const unsigned int node_number = 2;
    static const bool is_simplicial = true;

    typedef type_list<point> cell_type_of_submanifold;
    typedef type_list<holder<unsigned int, 2> > disjoint_submanifolds;
  };
    
  struct triangle {
    static const unsigned int dimension = 2;
    static const unsigned int node_number = 3;
    static const bool is_simplicial = true;
    
    typedef type_list<point, edge> cell_type_of_submanifold;
    typedef type_list<holder<unsigned int, 3>,
		      holder<unsigned int, 3> > disjoint_submanifolds;
  };

  typedef type_list<point, edge, triangle> simplex_list;
  
  template<unsigned int dimension>
  struct simplex_of_dim {
    typedef 
    typename at<dimension, simplex_list>::type
    type;
  };
}

#endif /* _CELL_H_ */
