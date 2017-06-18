
#include <iostream>
#include <vector>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"
#include "mesh.hpp"
#include "fes.hpp"

class sparse_linear_system {
public:
  sparse_linear_system(std::size_t n_equation, std::size_t n_unknown) {}
  void accumulate(std::size_t i, std::size_t j, double v) {}
};

mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n) {
  std::vector<double> vertices(n + 1);
  std::vector<unsigned int> elements(2 * n);

  for (unsigned int i(0); i < n + 1; ++i)
    vertices[i] = x_1 + (x_2 - x_1) / n * i;

  for (unsigned int k(0); k < n; ++k) {
    elements[2 * k] = k;
    elements[2 * k + 1] = k + 1;
  }

  return mesh<cell::edge>(&vertices[0], n + 1, 1,
			  &elements[0], n);
}


void laplacian() {
  

}

int main(int argc, char *argv[]) {
  const std::size_t n(10);
  mesh<cell::edge> m(gen_segment_mesh(0.0, 1.0, n));

  typedef finite_element::edge_lagrange_p1 fe;
  finite_element_space<fe> fes(m);

  fes.show(std::cout);


  sparse_linear_system sys(fes.get_dof_number(),
			   fes.get_dof_number());
  for (unsigned int k(0); k < m.get_element_number(); ++k) {
    for (unsigned int i(0); i < fe::n_dof_per_element; ++i) {
      for (unsigned int j(0); j < fe::n_dof_per_element; ++j) {

	double a_el(0.0);
	const std::size_t n_q(0);
	for (unsigned int q(0); q < n_q; ++q) {
	  // 1. Get the jacobian at x_q: the mesh should provide it,
	  // 2. Get the jmt at x_q: the mesh should provide it,
	  // 3. Evaluate the basis functions \phi_i at x_q: we need to define quadrature rules,
	  //                                                and the finite element provide \phi_i,
	  // 4. Evaluate the derivatives \partial_k \phi_i at x_q: the finite element provide \partial_k\phi_i
	  // 5. Get Q and \omega_q: the quadrature provides it.

	  const std::size_t dim(0);
	  const std::size_t n_dof(0);
	  
	  const double jac(0.0);
	  const array<double> xq{n_q};
	  const array<double> jmt{dim, dim};
	  const array<double> phi{n_dof}; // That can be cached before the loops
	  const array<double> dphi{n_dof}; // And that too
	  const double Q(0.0);
	  const array<double> omega{n_q};

	  for (unsigned int n(0); n < dim; ++n)
	    a_el += Q * omega.at(q) * jac * dphi.at(n, i) * dphi.at(n, j);
	}
	sys.accumulate(fes.get_dof(k, j),
		       fes.get_dof(k, i), a_el);
      }
    }
  }

  
  return 0;
}

