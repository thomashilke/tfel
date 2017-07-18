
#include <functional>

#include "../src/basic_fe_formulation.hpp"


template<typename fe>
class heat: public basic_fe_formulation<fe> {
public:
  using typename basic_fe_formulation<fe>::fe_type;
  using typename basic_fe_formulation<fe>::fes_type;
  using typename basic_fe_formulation<fe>::cell_type;
  using typename basic_fe_formulation<fe>::element_type;
  using typename basic_fe_formulation<fe>::volume_quadrature_type;
  using typename basic_fe_formulation<fe>::boundary_quadrature_type;
  
  heat(const mesh<cell_type>& m,
       double delta_t, double diffusivity)
    : delta_t(delta_t), diffusivity(diffusivity),
      m(m), dm(m.get_boundary_submesh()),
      fes(m, dm), a(fes, fes), f(fes),
      source(fes), solution(fes) {
    assemble_bilinear_form();
  }

  heat(const mesh<cell_type>& m,
       const submesh<cell_type>& dm,
       double delta_t, double diffusivity)
    : delta_t(delta_t), diffusivity(diffusivity),
      m(m), dm(dm),
      fes(m, dm), a(fes, fes), f(fes),
      source(fes), solution(fes) {
    assemble_bilinear_form();
  }

  void set_initial_condition(double (*ic)(const double* )) {
    solution = projector::l2<fe_type, volume_quadrature_type>(ic, fes);
  }

  void set_source(double (*src)(const double*)) {
    //source = projector::l2<fe_type, volume_quadrature_type>(src, fes);
    set_source(std::function<double(const double*)>(src));
  }

  void set_source(const std::function<double(const double*)>& src) {
    source = projector::l2<fe_type, volume_quadrature_type>(src, fes);
  }

  void step() {
    f.clear();

    const auto v(f.get_test_function());

    if (assemble_time_derivative)
      f += integrate<volume_quadrature_type>((1.0 / delta_t) * make_expr<fe_type>(solution) * v, m);

    if (assemble_source)
      f += integrate<volume_quadrature_type>(make_expr<fe_type>(source) * v, m);
    
    solution = a.solve(f);
  }

  const element_type& get_solution() const {
    return solution;
  }
  
private:
  const double delta_t, diffusivity;
  
  const mesh<cell_type>& m;
  const submesh<cell_type>& dm;
  fes_type fes;

  bilinear_form<fes_type, fes_type> a;
  linear_form<fes_type> f;

  static const bool assemble_time_derivative = true;
  static const bool assemble_laplacian = true;
  static const bool assemble_source = true;
  static const bool assemble_transport = true;
  
  element_type source;
  element_type solution;

private:
  static double velocity_x(const double* x) {
    return - std::sin(x[0] * M_PI) * std::cos(x[1] * M_PI);
    return - x[1] / M_PI / 0.5;
  }
  static double velocity_y(const double* x) {
    return std::cos(x[0] * M_PI) * std::sin(x[1] * M_PI);
    return  x[0] / M_PI / 0.5;
  }
  
  void assemble_bilinear_form() {
    //a.clear();
    
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());

    if (assemble_time_derivative)
      a += integrate<volume_quadrature_type>((1.0 / delta_t) * u * v, m);

    if (assemble_laplacian)
      a += integrate<volume_quadrature_type>(diffusivity * (d<1>(u) * d<1>(v)
							    + d<2>(u) * d<2>(v)),
					     m);
    if (assemble_transport)
      a += integrate<volume_quadrature_type>((velocity_x * d<1>(u) + velocity_y * d<2>(u)) * v, m);
  }
};

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

  typedef finite_element::triangle_lagrange_p1 fe_type;
  typedef mesh<fe_type::cell_type> mesh_type;
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

  heat<fe_type> h(m, dirichlet_boundary, delta_t, D);
  h.set_initial_condition(initial_condition);

  exporter::ensight6_transient<fe_type> ens("heat", m, "solution");
  ens.export_time_step(time, h.get_solution());

  for (std::size_t i(0); i < M; ++i) {
    std::cout << "step" << i << std::endl;

    time += delta_t;
    h.set_source(source_function);
    h.step();

    ens.export_time_step(time, h.get_solution());
  }

  return 0;
}


