#include "../src/formulations/stokes_2d.hpp"

double f_0(const double* x) { return 0.0; }
double f_1(const double* x) { return 0.0; }
double u_0_bc(const double* x) {
    return x[1] > 0.999 and x[0] != 0.0 and x[0] != 1.0? 1.0 : 0.0;
}
double u_1_bc(const double* x) {
    return 0.0;
}


int main(int argc, char *argv[]) {
  using cell_type = cell::triangle;
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 50, 50));

  stokes_2d<> s2d(m, 0.1);
  s2d.set_velocity_condition(u_0_bc, u_1_bc);
  s2d.set_force(f_0, f_1);
  s2d.solve();

  const auto solution(s2d.get_solution());
  exporter::ensight6("stokes_p",
		     solution.get_component<0>(), "u0",
		     solution.get_component<1>(), "u1",
		     solution.get_component<2>(), "p");

  {
    finite_element_space<finite_element::triangle_lagrange_p0> fes(m);
    const auto h(build_element_diameter_function(m, fes));
    exporter::ensight6("diam", h, "h");
  }
  
  return 0;
}
