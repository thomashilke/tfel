#include "../../formulations/stokes_2d.hpp"


/*
 *  Force terms
 */
double f_0(const double* x) { return 0.0; }

double f_1(const double* x) { return 0.0; }


/*
 *  Boundary values
 */
double u_0_bc(const double* x) {
    return x[1] > 0.999 and x[0] != 0.0 and x[0] != 1.0? 1.0 : 0.0;
}

double u_1_bc(const double* x) {
    return 0.0;
}


/*
 *  Compute the Stokes flow in a square cavity, zero force
 *  terms, zero velocity on the boundaries, exept on the top
 *  boundary where we take the velocity uniformly u = (1, 0).
 *
 *  The solution is exported in the file
 *  stokes_2d_driven_cavity.case .
 */
int main(int argc, char *argv[]) {
  try {
    using cell_type = cell::triangle;
    mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 50, 50));

    using vfe = finite_element::triangle_lagrange_p1_bubble;
    //using vfe = finite_element::triangle_lagrange_p1;
    using pfe = finite_element::triangle_lagrange_p1;

    stokes_2d<vfe, pfe> s2d(m, 1.0);
    s2d.set_velocity_condition(u_0_bc, u_1_bc);
    s2d.set_force(f_0, f_1);
    s2d.solve();

    const auto solution(s2d.get_solution());
    exporter::ensight6("stokes_2d_driven_cavity",
		       solution.get_component<0>(), "u0",
		       solution.get_component<1>(), "u1",
		       solution.get_component<2>(), "p");
  }
  catch (const std::string& exept) {
    std::cout << exept << std::endl;
  }
  
  return 0;
}
