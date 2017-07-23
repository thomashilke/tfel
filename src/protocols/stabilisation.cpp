
#include "../formulations/steady_advection_1d.hpp"

int main(
  int argc,
  char *argv[])
{
  using cell_type = cell::edge;
  using fe_type = finite_element::edge_lagrange_p1;
  mesh<cell_type> m(gen_segment_mesh(0.0, 1.0, 64));
  submesh<cell_type> dm(m.get_boundary_submesh());

  sumbesh<cell_typ>
    left_boundary(dm.query_elements([](const double* x) { return x[0] < 0.5;} )),
    right_boundary(dm.query_elements([](const double* x) { return x[0] > 0.5;} ));

  steady_advection_diffusion<fe_type> sad(m, left_boundary, 1.0);
  sad.set_boundary_value([](const double* x) { return x[0]; });
  sad.set_advection_velocity([](const double* x) { return 1.0; });
  sad.solve();

  exporter::ascii("sad.dat", sad.get_solution());
  
  return 0;
}
