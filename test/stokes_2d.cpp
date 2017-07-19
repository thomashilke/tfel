#include "../src/basic_fe_formulation.hpp"

using stokes_2d_fe = composite_finite_element<finite_element::triangle_lagrange_p1,
					      finite_element::triangle_lagrange_p1,
					      finite_element::triangle_lagrange_p1>;

class stokes_2d
  : public basic_fe_formulation<stokes_2d_fe> {
public:
  using basic_fe_formulation<stokes_2d_fe>::fe_type;
  using basic_fe_formulation<stokes_2d_fe>::fes_type;
  using basic_fe_formulation<stokes_2d_fe>::cell_type;
  using basic_fe_formulation<stokes_2d_fe>::element_type;
  using basic_fe_formulation<stokes_2d_fe>::volume_quadrature_type;
  using basic_fe_formulation<stokes_2d_fe>::boundary_quadrature_type;

  using velocity_fe_type = fe_type::fe_type<0>;
  using pressure_fe_type = fe_type::fe_type<2>;

  using velocity_fes_type = finite_element_space<velocity_fe_type>;
  using pressure_fes_type = finite_element_space<pressure_fe_type>;
  
  stokes_2d(const mesh<cell_type>& m, double mu)
    : viscosity(mu),
      m(m), dm(m.get_boundary_submesh()),
      pressure_point_m(m.get_point_submesh(m.get_element_number() / 2)),
      fes(m),
      v_fes(fes.get_finite_element_space<0>()),
      p_fes(fes.get_finite_element_space<2>()),
      a(fes, fes), f(fes),
      source(fes.get_finite_element_space<0>()),
      solution(fes) {
    
    fes.set_dirichlet_boundary_condition<0, cell::edge>(dm, u_0_bc);
    fes.set_dirichlet_boundary_condition<1, cell::edge>(dm, u_1_bc);
    fes.set_dirichlet_boundary_condition<2, cell::point>(pressure_point_m);
    
    assemble_bilinear_form();
    assemble_linear_form();
  }

  void set_source(const std::function<double(const double*)>& src) {
    source = projector::l2<velocity_fe_type, volume_quadrature_type>(src, v_fes);
  }

  void solve() {
    solution = a.solve(f);
  }
  
  const element_type& get_solution() const {
    return solution;
  }

private:
  double viscosity;

  const mesh<cell_type>& m;
  const submesh<cell_type> dm;
  const submesh<cell_type, cell::point> pressure_point_m;
  
  fes_type fes;
  const velocity_fes_type& v_fes;
  const pressure_fes_type& p_fes;
  
  bilinear_form<fes_type, fes_type> a;
  linear_form<fes_type> f;

  velocity_fes_type::element source;
  fes_type::element solution;

private:
  void assemble_bilinear_form() {

    const double h(1.0 / 100.0);
    const double stab_coefficient(10.0);
    
    const auto u_0(a.get_trial_function<0>());
    const auto u_1(a.get_trial_function<1>());
    const auto p(a.get_trial_function<2>());
    
    const auto v_0(a.get_test_function<0>());
    const auto v_1(a.get_test_function<1>());
    const auto q(a.get_test_function<2>());

    a += integrate<volume_quadrature_type>(  viscosity * (  d<1>(u_0)*d<1>(v_0) + d<2>(u_0)*d<2>(v_0)
							  + d<1>(u_1)*d<1>(v_1) + d<2>(u_1)*d<2>(v_1))
					     - p * (d<1>(v_0) + d<2>(v_1))
					     - q * (d<1>(u_0) + d<2>(u_1))
					     - stab_coefficient * (h * h) * (d<1>(p) * d<1>(q) + d<2>(p) * d<2>(q)),
					   m);
  }

  static double f_0(const double* x) { return 0.0; }
  static double f_1(const double* x) { return 0.0; }
  static double f_2(const double* x) { return 0.0; }

  static double u_0_bc(const double* x) { return x[0]; }
  static double u_1_bc(const double* x) { return -x[1]; }
  
  void assemble_linear_form() {
    const auto v_0(f.get_test_function<0>());
    const auto v_1(f.get_test_function<1>());
    const auto q(f.get_test_function<2>());

    f += integrate<volume_quadrature_type>(  f_0 * v_0
					   + f_1 * v_1
					   , m);
  }
};



int main(int argc, char *argv[]) {
  using cell_type = cell::triangle;
  
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 100, 100));

  stokes_2d s2d(m, 1.0);
  //s2d.solve();
  return 0;
  const auto solution(s2d.get_solution());
  exporter::ensight6<stokes_2d::velocity_fe_type>("stokes_u_0", solution.get_component<0>(), "u_0");
  exporter::ensight6<stokes_2d::velocity_fe_type>("stokes_u_1", solution.get_component<1>(), "u_1");
  exporter::ensight6<stokes_2d::velocity_fe_type>("stokes_p", solution.get_component<2>(), "p");
  
  return 0;
}
