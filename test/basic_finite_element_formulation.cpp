
#include "../src/mesh.hpp"
#include "../src/fe.hpp"
#include "../src/fes.hpp"
#include "../src/quadrature.hpp"
#include "../src/projector.hpp"
#include "../src/export.hpp"
#include "../src/timer.hpp"
#include "../src/expression.hpp"

template<typename fe>
class basic_fe_formulation {
public:
  typedef fe fe_type;
  typedef typename fe_type::cell_type cell_type;
  typedef finite_element_space<fe_type> fes_type;
  typedef typename fes_type::element element_type;
  typedef typename default_quadrature<cell_type>::type volume_quadrature_type;
  typedef typename default_quadrature<typename cell_type::boundary_cell_type>::type boundary_quadrature_type;

  virtual ~basic_fe_formulation() {}
  
private:
};


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
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());

    a += integrate<volume_quadrature_type>((1.0 / delta_t) * u * v, m);
    a += integrate<volume_quadrature_type>(diffusivity * (d<1>(u) * d<1>(v)
							  + d<2>(u) * d<2>(v)),
					   m);
  }

  void set_initial_condition(double (*ic)(const double* )) {
    solution = projector::l2<fe_type, volume_quadrature_type>(ic, fes);
  }

  void set_source(double (*src)(const double*)) {
    source = projector::l2<fe_type, volume_quadrature_type>(src, fes);
  }

  void step() {
    f.clear();

    const auto v(f.get_test_function());
    f += integrate<volume_quadrature_type>((1.0 / delta_t) * make_expr<fe_type>(solution) * v, m);
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
  
  element_type solution;
  element_type source;
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
  const double x_0(0.5), x_1(0.5);
  const double r(std::sqrt(std::pow(x[0] - x_0, 2) + std::pow(x[1] - x_1, 2)));
  return 10.0 * shift_scale_bell(r, 0.15, 0.0);
}

int main(int argc, char *argv[]) {

  typedef finite_element::triangle_lagrange_p1 fe_type;
  typedef mesh<fe_type::cell_type> mesh_type;

  unsigned int M(200);
  double delta_t(0.1);
  double D(0.1);
  std::size_t n(50);

  mesh_type m(gen_square_mesh(1.0, 1.0, n, n));

  heat<fe_type> h(m, delta_t, D);
  h.set_initial_condition(initial_condition);
  h.set_source(heat_source);

  double time(0.0);
  exporter::ensight6_transient<fe_type> ens("heat", m, "solution");
  ens.export_time_step(time, h.get_solution());

  for (std::size_t i(0); i < M; ++i) {
    std::cout << "step" << i << std::endl;
    h.step();

    time += delta_t;
    ens.export_time_step(time, h.get_solution());
  }

  return 0;
}


