
#include <iostream>

#include "../src/core/mesh.hpp"
#include "../src/core/mesh_data.hpp"
#include "../src/core/export.hpp"

double b_0(const double* x) {
  //return 1.0;
  return - (x[1] - 0.5);
  return std::sin(M_PI * x[0]) * std::cos(M_PI * x[1]);
}

double b_1(const double* x) {
  //return 1.0;
  return (x[0] - 0.5);
  return - std::sin(M_PI * x[1]) * std::cos(M_PI * x[0]);
}

int main(int argc, char *argv[]) {
  mesh<cell::triangle> m(gen_square_mesh(1.0, 1.0, 4, 4));
  m.show(std::cout);

  submesh<cell::triangle> dm(m.get_boundary_submesh());
  dm.show(std::cout);

  mesh_data<double, submesh<cell::triangle> > n(compute_boundary_normals(dm));
  mesh_data<double, submesh<cell::triangle> > b_bc(evaluate_on_cells(dm, b_0, b_1));
  mesh_data<double, submesh<cell::triangle> > b_v(evaluate_on_vertices(dm, b_0, b_1));

  mesh_data<bool, submesh<cell::triangle> > inflow_boundary_cells(
    dm, mesh_data<bool, submesh<cell::triangle> >::data_type::cell,
    (n[0] * b_bc[0] + n[1] * b_bc[1]) < -1.0e-6 );

  //submesh<cell::triangle, cell::edge> inflow_m(dm.submesh_from_selection(inflow_boundary_cells));
  mesh<cell::edge> inflow_m(dm.submesh_from_selection(inflow_boundary_cells));
  inflow_m.show(std::cout);

  exporter::ensight6_geometry("square", m);
  exporter::ensight6_geometry("inflow", inflow_m);
  
  return 0;
}
