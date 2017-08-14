#include "../../formulations/unsteady_advection_diffusion_2d.hpp"

double b_0(const double* x) { return 1.0 / 2.0; }
double b_1(const double* x) { return 0.0; }


double u_bc(const double* x) {
  return 1.0;
}

double ic(const double* x) {
  return std::tanh(- (x[0] - 0.25) * 20.0);
}

double exact_solution(double time, const double* x) {
  const double b(1.0 / 2.0);
  const double xp[] = {x[0] - b * time, x[1]};
  return ic(xp);
}

double sqr(double x) { return x * x; }

double error(std::size_t M, std::size_t N);

int main(int argc, char *argv[]) {
  std::vector<std::size_t>
    Ms{64, 181, 512, 1449, 4096},
    Ns{16,  32,  64,  128,  256};

  std::cout << "# delta_t h l2_error" << std::endl;
  for (std::size_t i(0); i < Ms.size(); ++i) {
    std::size_t N(Ns[i]), M(Ms[i]);
    const double err(error(M, N));
    std::cout << 1.0 / M << " " << 1.0 / N << " " << err << std::endl;
  }
}

double error(std::size_t M, std::size_t N) {
  const double diffusivity(0.0);
  const double t_end(1.0);
  const double delta_t(t_end / M);
    
  
  using cell_type = cell::triangle;
  using mesh_type = fe_mesh<cell_type>;
  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, N, N));
  submesh<cell_type> dm(m.get_boundary_submesh());

  /*
    submesh<cell_type> inflow_boundary(dm.inflow_boundary([&](const double* x) -> array<double> {
      array<double> b{2};
      b.at(0) = vortex::b_0(x);
      b.at(1) = vortex::b_1(x);
      return b;
    }));
  */

  submesh<cell_type> inflow_boundary(dm.query_cells([](const double* x) {
	return x[0] < 0.0001;
      }));
  
  using fe_type = cell::triangle::fe::lagrange_p1;
  unsteady_advection_diffusion<fe_type> tad(m, inflow_boundary, delta_t, diffusivity);
  tad.set_boundary_value(u_bc);
  tad.set_advection_velocity(b_0, b_1);
  tad.set_initial_condition(ic);


  exporter::ensight6_transient<mesh_type>
    ens("unsteady_advection_diffusion_rotating_hill",
    m, "solution");

  double time(0.0);
  ens.export_time_step(time, to_mesh_vertex_data<fe_type>(tad.get_solution()));
  for (std::size_t k(0); k < M; ++k) {
    std::cerr << "step " << k << std::endl;
    time += delta_t;
    tad.step();
    ens.export_time_step(time, to_mesh_vertex_data<fe_type>(tad.get_solution()));
  }


  finite_element_space<fe_type> fes(m);
  exporter::ensight6("exact_sol"
		     ,
		     evaluate_on_vertices<fe_mesh<cell_type> >(m, [=](const double* x) -> double { return exact_solution(time, x) ;})
		     //projector::lagrange<fe_type>(, fes)
		     , "solution");
  
  return std::sqrt(integrate<quad::triangle::qf5pT>(compose(sqr,
							    (make_expr(std::bind(exact_solution, time, std::placeholders::_1))
							     - make_expr<fe_type>(tad.get_solution()))), m));
}
