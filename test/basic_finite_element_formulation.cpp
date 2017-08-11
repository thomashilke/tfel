
#include <functional>

#include "../src/formulations/unsteady_diffusion_2d.hpp"


double bell(double x) {
  const double epsilon(1e-5);
  if (std::abs(x) < 1.0 - epsilon) {
    return std::exp(1.0 - 1.0 / (1.0 - x * x));
  } else {
    return 0.0;
  }
}

double shift_scale_bell(double x, double width, double x_0) {
  return bell((x - x_0) / width);
}

double initial_condition(const double* x) { return 0.0; }
double heat_source(const double* x) {
  const double x_0(0.5), x_1(0.75);
  const double r(std::sqrt(std::pow(x[0] - x_0, 2) + std::pow(x[1] - x_1, 2)));
  return 10.0 * shift_scale_bell(r, 0.15, 0.0);
}



int main(int argc, char *argv[]) {

  typedef cell::triangle::fe::lagrange_p1 fe_type;
  typedef fe_mesh<fe_type::cell_type> mesh_type;
  typedef submesh<fe_type::cell_type> submesh_type;

  unsigned int M(200);
  double delta_t(0.1);
  double t_end(M * delta_t);
  double D(0.001);
  std::size_t n(50);

  double time(0.0);
  auto source_function([&](const double* x) -> double {
      if (time > t_end / 2)
	return 0.0;
      else
	return heat_source(x);
    });
  
  mesh_type m(gen_square_mesh(1.0, 1.0, n, n));
  submesh_type dm(m.get_boundary_submesh());
  submesh_type dirichlet_boundary(dm.query_elements([&](const double* x) -> bool {
	  return not ((x[1] < 1.0 / (2.0 * n))
		      and x[0] > 1.0 / 3.0
		      and x[0] < 2.0 / 3.0);
		      }));

  unsteady_diffusion_2d<fe_type> h(m, dirichlet_boundary, delta_t, D);
  h.set_initial_condition(initial_condition);

  exporter::ensight6_transient<mesh_type> ens("heat", m, "solution");
  ens.export_time_step(time, to_mesh_vertex_data<fe_type>(h.get_solution()));

  for (std::size_t i(0); i < M; ++i) {
    std::cout << "step" << i << std::endl;

    time += delta_t;
    h.set_source(source_function);
    h.step();

    ens.export_time_step(time, to_mesh_vertex_data<fe_type>(h.get_solution()));
  }

  return 0;
}


