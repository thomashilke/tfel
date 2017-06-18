#ifndef _FE_H_
#define _FE_H_

#include "cell.hpp"

namespace finite_element {

  struct edge_lagrange_p1 {
    typedef cell::edge cell_type;

    static const std::size_t n_node_per_element = 2;
 
    static double basis_function(unsigned int i,
				 unsigned int* derivatives,
				 double* x) {
      typedef double (*bf_type)(double*);
      static const bf_type bf[2][2] = {{bf_0_1, bf_0_2},
				       {bf_1_1, bf_1_2}};
      if (i >= 2)
	throw std::string("finite_element::edge_lagrange_p1::basis_function: "
			  "i out of range");

      if (derivatives[0] > 1)
	return 0.0;

      return bf[derivatives[0]][i](x);
    }

  protected:
    static double bf_0_1(double* x) { return 1.0 - x[0]; }
    static double bf_0_2(double* x) { return       x[0]; }
    static double bf_1_1(double* x) { return - 1.0; }
    static double bf_1_2(double* x) { return   1.0; }
  };

  
  struct edge_lagrange_p1_bubble: public edge_lagrange_p1 {
    using edge_lagrange_p1::cell_type;

    static double basis_function(unsigned int i,
				 unsigned int* derivatives,
				 double* x) {
      typedef double (*bf_type)(double*);
      static const bf_type bf[2][3] = {{edge_lagrange_p1::bf_0_1, edge_lagrange_p1::bf_0_2, bf_0_3},
				       {edge_lagrange_p1::bf_1_1, edge_lagrange_p1::bf_1_2, bf_1_3}};
      if (i >= 3)
	throw std::string("finite_element::edge_lagrange_p1::basis_function: "
			  "i out of range");

      if (derivatives[0] > 1) {
	if (i < 2)
	  return 0.0;
	else
	  return bf_2_3(x);
      }

      return bf[derivatives[0]][i](x);
    }
    
  protected:
    static double bf_0_3(double* x) {
      return edge_lagrange_p1::bf_0_1(x) * edge_lagrange_p1::bf_0_2(x);
    }
    
    static double bf_1_3(double* x) {
      return edge_lagrange_p1::bf_1_1(x) * edge_lagrange_p1::bf_0_2(x)
	+ edge_lagrange_p1::bf_0_1(x) * edge_lagrange_p1::bf_1_2(x);
    }
    
    static double bf_2_3(double* x) {
      return edge_lagrange_p1::bf_1_1(x) * edge_lagrange_p1::bf_1_2(x)
	+ edge_lagrange_p1::bf_1_1(x) * edge_lagrange_p1::bf_1_2(x);
    }
  };

  struct triangle_lagrange_p1 {
    typedef cell::triangle cell_type;

    static const std::size_t n_node_per_element = 3;

    static double basis_function(unsigned int i,
				 unsigned int* derivative,
				 double* x) {
      return 0.0;
    }

  protected:
    static double bf_0_1(double* x) { return 1.0 - x[0] - x[1]; }
    static double bf_0_2(double* x) { return x[0]; }
    static double bf_0_3(double* x) { return x[1]; }
    static double bf_10_1(double* x) { return -1.0; }
    static double bf_10_2(double* x) { return  0.0; }
    static double bf_10_3(double* x) { return  1.0; }
    static double bf_01_1(double* x) { return -1.0; }
    static double bf_01_2(double* x) { return  0.0; }
    static double bf_01_3(double* x) { return  1.0; }
  };
  
  /*template<unsigned int dimension>
  struct lagrange_p1 {
    typedef
    typename cell::simplex_of_dim<dimension>::type
    cell_type;
    };*/

}

#endif /* _FE_H_ */
