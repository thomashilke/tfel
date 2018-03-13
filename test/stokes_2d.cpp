#include "../src/formulations/stokes_2d.hpp"

/*
 *  Regularized driven cavity
 */
/*
double f0(const double* x) { return 0.0; }
double f1(const double* x) { return 0.0; }

double u0_bv(const double* x) {
  return ((x[1] > 0.999
           and x[0] != 0.0
           and x[0] != 1.0)
          ? x[0] * (1.0 - x[0]) : 0.0);
}

double u1_bv(const double* x) {
  return 0.0;
}
*/

/*
 *  Exact polynomial solution
 */
double f0(const double* x) { return 0.0; }
double f1(const double* x) { return 0.0; }

double u0_bv(const double* x) {
  return x[0];
}

double u1_bv(const double* x) {
  return -x[1];
}


int main(int argc, char *argv[]) {
  try {
    using cell_type = cell::triangle;
    using mesh_type = fe_mesh<cell_type>;
    using velocity_fe_type = cell_type::fe::lagrange_p1_bubble;
    //using velocity_fe_type = cell_type::fe::lagrange_p1;
    using pressure_fe_type = cell_type::fe::lagrange_p1;
    
    mesh_type m(gen_square_mesh(1.0, 1.0, 100, 100));

    stokes_2d<velocity_fe_type, pressure_fe_type> s2d(m, 0.1);
    s2d.set_velocity_condition(u0_bv, u1_bv);
    s2d.set_force(f0, f1);
    s2d.solve();
  
    const auto solution(s2d.get_solution());
    exporter::ensight6("stokes",
  		       mesh_data<double, mesh_type>(m, mesh_data_kind::vertex,
  		       to_mesh_vertex_data<velocity_fe_type>(solution.get_component<0>()),
                       to_mesh_vertex_data<velocity_fe_type>(solution.get_component<1>())), "u",
  		       to_mesh_vertex_data<pressure_fe_type>(solution.get_component<2>()),  "p");
  
    {
      finite_element_space<cell::triangle::fe::lagrange_p0> fes(m);
      const auto h(build_element_diameter_function(m, fes));
      exporter::ensight6("diam", to_mesh_cell_data<cell::triangle::fe::lagrange_p0>(h), "h");
    }
    
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
