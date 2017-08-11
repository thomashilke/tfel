
#include "../../formulations/steady_advection_diffusion_1d.hpp"

int main(
  int argc,
  char *argv[])
{
  using cell_type = cell::edge;
  using mesh_type = fe_mesh<cell_type>;
  using fe_type = cell::edge::fe::lagrange_p1;
  fe_mesh<cell_type> m(gen_segment_mesh(0.0, 1.0, 40));
  submesh<cell_type> dm(m.get_boundary_submesh());

  submesh<cell_type>
    left_boundary(dm.query_elements([](const double* x) { return x[0] < 0.5;} )),
    right_boundary(dm.query_elements([](const double* x) { return x[0] > 0.5;} ));

  steady_advection_diffusion<fe_type> sad(m, dm, 1.0/50.0);
  sad.set_boundary_value([](const double* x) { return x[0]; });
  sad.set_advection_velocity([](const double* x) { return 1.0; });
  sad.solve();

  exporter::ascii<mesh_type>("sad.dat", to_mesh_vertex_data<fe_type>(sad.get_solution()));
  
  return 0;
}
