
#include <cmath>

#include "../src/core/mesh.hpp"
#include "../src/core/fe.hpp"
#include "../src/core/fes.hpp"
#include "../src/core/quadrature.hpp"
#include "../src/core/projector.hpp"
#include "../src/core/export.hpp"
#include "../src/core/timer.hpp"
#include "../src/core/expression.hpp"

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


double velocity_x(const double* x) { return -x[1]; }
double velocity_y(const double* x) { return  x[0]; }


double initial_condition(const double* x) {
  //return 0.0;
  return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  //const double x_0(0.5), x_1(0.5);
  //const double r(std::sqrt(std::pow(x[0] - x_0, 2) + std::pow(x[1] - x_1, 2)));
  //return shift_scale_bell(r, 0.25, 0.0);
}

double source(const double* x) {
  //return (-1.0 + 2.0 * M_PI * M_PI) * initial_condition(x);
  return 0.0;
}

double bc(const double* x) {
  return 0.0;
  return x[0] * x[0] + x[1] * x[1];
}

int main(int argc, char *argv[]) {
  typedef cell::triangle cell_type;
  typedef cell::triangle::fe::lagrange_p1 fe_type;
  typedef finite_element_space<fe_type> fes_type;
  typedef typename finite_element_space<fe_type>::element element_type;
  typedef quad::triangle::qf5pT q_type;
  using mesh_type = mesh<cell_type>;
  
  const std::size_t n(100);
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());

  fes_type fes(m, dm, bc);
  const element_type c_init(projector::l2<fe_type,
			                  q_type>(initial_condition,
						  fes));
  exporter::ensight6("initial_condition", to_mesh_vertex_data<fe_type>(c_init), "c_init");


  const std::size_t M(100);
  const double delta_t(2.0 / M);
  const double diffusivity(1.0);

  bilinear_form<fes_type, fes_type> a(fes, fes); {
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());
    a += integrate<q_type>(diffusivity * (d<1>(u) * d<1>(v) + d<2>(u) * d<2>(v)), m);
    a += integrate<q_type>((1.0 / delta_t) * u * v, m);
  }

  element_type c(c_init);

  exporter::ensight6_transient<mesh_type> ens("solution", m, "c");
  ens.export_time_step(0.0, to_mesh_vertex_data<fe_type>(c));
  
  for (std::size_t k(0); k < M; ++k) {
    std::cout << "step " << k << std::endl;
    linear_form<fes_type> f(fes); {
      const auto v(f.get_test_function());
      f += integrate<q_type>((1.0 / delta_t) * make_expr<fe_type>(c) * v
			     - initial_condition * v
			     + (1.0 - (k + 1) * delta_t) * 2 * M_PI * M_PI * (initial_condition * v), m);
    }
    c = a.solve(f);
    ens.export_time_step((k + 1) * delta_t, to_mesh_vertex_data<fe_type>(c));
  }
  
  return 0;
}
