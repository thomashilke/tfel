#include "../../formulations/stationary_advection_diffusion_2d.hpp"

double b_0(const double* x) {
  return 0.0;
}

double b_1(const double* x) {
  return 1.0;
}

double u_bc(const double* x) {
  return x[0] > 0.5 ? 1.0 : 0.0;
}

int main(int argc, char *argv[]) {
  const double diffusivity(1e-3);
  
  using cell_type = cell::triangle;
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 30, 30));
  
  using fe_type = finite_element::triangle_lagrange_p1;
  stationary_advection_diffusion<fe_type> sad(m, diffusivity);

  sad.set_boundary_value(u_bc);
  sad.set_advection_velocity(b_0, b_1);
  sad.solve();

  const auto solution(sad.get_solution());

  exporter::ensight6("stationary_advection_diffusion_2d",
		     solution, "u");
  
  return 0;
}
