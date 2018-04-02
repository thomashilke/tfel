#include "../src/core/fes.hpp"
#include "../src/core/projector.hpp"

double f(const double* x) {
  return 0.5 * x[1] * x[1];
  //return x[1];
  //return x[0];
}

int main() {

  using cell_type = cell::triangle;
  using fe_type = cell_type::fe::lagrange_p2;
  using fes_type = finite_element_space<fe_type>;
  using fes_elem_type = fes_type::element;
  using quad_type = quad::triangle::qf5pT;
  
  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 10, 10));
  fes_type fes(m);

  fes_elem_type fh(projector::lagrange(f, fes));

  std::cout << "int f = "
            << integrate<quad_type>(make_expr<fe_type>(fh), m) << std::endl;
  std::cout << "int df/dx = "
            << integrate<quad_type>(d<2>(make_expr<fe_type>(fh)), m) << std::endl;

  
  return 0;
}
