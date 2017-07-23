#include "../../formulations/transient_advection_diffusion_2d.hpp"


double bell(double r_sqr) {
  const double epsilon(1e-5);
  if (std::abs(r_sqr) < 1.0 - epsilon) {
    return std::exp(1.0 - 1.0 / (1.0 - r_sqr));
  } else {
    return 0.0;
  }
}

double shift_scale_bell(const double* x, double width, double x_0, double x_1) {
  return bell((std::pow(x[0] - x_0, 2) +
	       std::pow(x[1] - x_1, 2)) / (width * width));
}


struct vortex {
  static double b_0(const double* x) {
    return - (x[1] - 0.5) * M_PI;
  }

  static double b_1(const double* x) {
    return (x[0] - 0.5) * M_PI;
  }
};


struct kim_moin {
  static double b_0(const double* x) {
    return std::sin(x[0] * M_PI) * std::cos(x[1] * M_PI);
  }

  static double b_1(const double* x) {
    return - std::cos(x[0] * M_PI) * std::sin(x[1] * M_PI);
  }
};


double u_bc(const double* x) {
  return 0.0;
}

double exact_solution(double time, const double* x) {
  double xp[] = {std::cos(-M_PI * time) * x[0] + std::sin(-M_PI * time) * x[1],
		 -std::sin(-M_PI * time) * x[0] + std::cos(-M_PI * time) * x[1]};
  return shift_scale_bell(xp, 0.10, 0.5, 0.75);
}


int main(int argc, char *argv[]) {
  const double diffusivity(1e-7);
  const double t_end(0.5);
  const std::size_t M(16), N(16);
  const double delta_t(t_end / M);
    
  
  using cell_type = cell::triangle;
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, N, N));
  
  using fe_type = finite_element::triangle_lagrange_p1;
  transient_advection_diffusion<fe_type> tad(m, delta_t, diffusivity);
  tad.set_advection_velocity(vortex::b_0, vortex::b_1);
  tad.set_initial_condition([](const double* x) {
      return shift_scale_bell(x, 0.10, 0.5, 0.75);
    });

  tad.set_boundary_value(u_bc);

  exporter::ensight6_transient<fe_type>
    ens("transient_advection_diffusion_rotating_hill",
	m, "solution");

  double time(0.0);
  ens.export_time_step(time, tad.get_solution());
  for (std::size_t k(0); k < M; ++k) {
    std::cout << "step " << k << std::endl;
    time += delta_t;
    tad.step();
    ens.export_time_step(time, tad.get_solution());
  }
  
  return 0;
}
