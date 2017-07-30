
#include "../../formulations/steady_advection_diffusion_1d.hpp"

int main(
  int argc,
  char *argv[])
{
  using cell_type = cell::edge;
  using fe_type = finite_element::edge_lagrange_p1;
  mesh<cell_type> m(gen_segment_mesh(0.0, 1.0, 40));
  submesh<cell_type> dm(m.get_boundary_submesh());

  submesh<cell_type>
    left_boundary(dm.query_elements([](const double* x) { return x[0] < 0.5;} )),
    right_boundary(dm.query_elements([](const double* x) { return x[0] > 0.5;} ));

  steady_advection_diffusion<fe_type> sad(m, dm, 1.0/50.0);
  sad.set_boundary_value([](const double* x) { return x[0]; });
  sad.set_advection_velocity([](const double* x) { return 1.0; });
  sad.solve();

  exporter::ascii<fe_type>("sad.dat", sad.get_solution());
  
  return 0;
}
