
#include "../src/mesh.hpp"
#include "../src/fe.hpp"
#include "../src/fes.hpp"
#include "../src/quadrature.hpp"
#include "../src/projector.hpp"

double func(const double* x) {
  return std::sin(x[0] * M_PI) * std::sin(x[1] * M_PI);
}

int main(int , char**) {

  mesh<cell::triangle> m(gen_square_mesh(1.0, 1.0, 3, 3));
  typedef finite_element::triangle_lagrange_p1 fe_type;
  finite_element_space<fe_type> fes(m);

  const auto func_h(projector::l2<fe_type, quad::triangle::qf1pT>(func, fes));
  
  return 0;
}
