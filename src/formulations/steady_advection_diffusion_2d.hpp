#ifndef _STEADY_ADVECTION_DIFFUSION_2D_H_
#define _STEADY_ADVECTION_DIFFUSION_2D_H_

#include "../core/basic_fe_formulation.hpp"


/*
 *  The advection field b is supposed to be in H(div).
 */
template<typename fe>
class steady_advection_diffusion: public basic_fe_formulation<fe> {
public:
  using fe_type = typename basic_fe_formulation<fe>::fe_type;
  using fes_type = typename basic_fe_formulation<fe>::fes_type;
  using cell_type = typename basic_fe_formulation<fe>::cell_type;
  using element_type = typename basic_fe_formulation<fe>::element_type;
  using volume_quadrature_type = typename basic_fe_formulation<fe>::volume_quadrature_type;
  using boundary_quadrature_type = typename basic_fe_formulation<fe>::boundary_quadrature_type;


  steady_advection_diffusion(const fe_mesh<cell_type>& m, double diffusivity)
    : diffusivity(diffusivity),
      m(m), dm(m.get_boundary_submesh()), fes(m, dm),
      a(fes, fes), f(fes),
      solution(fes),
      b_0(null_function), b_1(null_function), src(null_function) {}

  void set_boundary_value(const std::function<double(const double*)>& u_bc) {
    fes.set_dirichlet_condition(u_bc);
  }
  
  void set_advection_velocity(const std::function<double(const double*)>& b_0,
			      const std::function<double(const double*)>& b_1) {
    this->b_0 = b_0;
    this->b_1 = b_1;
    assemble_bilinear_form();
  }

  void set_source_term(const std::function<double(const double*)>& src) {
    this->src = src;
  }
  
  void solve() {
    assemble_linear_form();
    solution = a.solve(f);
  }

  element_type get_solution() const { return solution; }
      
private:
  const double diffusivity;
  
  const fe_mesh<cell_type>& m;
  const submesh<cell_type> dm;
  
  fes_type fes;

  bilinear_form<fes_type, fes_type> a;
  linear_form<fes_type> f;

  element_type solution;
  
  std::function<double(const double*)> b_0, b_1, src;

  double diffusion_stabilisation = false;
  double supg_stabilisation = true;

private:
  static double null_function(const double* x) { return 0.0; }
  
  void assemble_bilinear_form() {
    a.clear();
    
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());

    a += integrate<volume_quadrature_type>(diffusivity * (d<1>(u) * d<1>(v) + d<2>(u) * d<2>(v)) +
					   make_expr(b_0) * d<1>(u) * v +
					   make_expr(b_1) * d<2>(u) * v
					   , m);
    if (diffusion_stabilisation) {

      const double h(1.0/3.0);
      const double b_norm(1.0);
      a += integrate<volume_quadrature_type>(h * b_norm * (d<1>(u) * d<1>(v) + d<2>(u) * d<2>(v))
					     , m);
    }

    if (supg_stabilisation) {
      if (not std::is_same<fe_type, cell::triangle::fe::lagrange_p1>::value)
	throw std::string("steady advection diffusion 2d:"
			  " supg stabilisation is not available for non piece wise linear finite elements");
      
      const double delta(1.0);
      const double b_norm(1.0);
      
      finite_element_space<cell::triangle::fe::lagrange_p0> p0_fes(m);
      const auto h(build_element_diameter_function(m, p0_fes));

      a += integrate<volume_quadrature_type>(delta / b_norm * make_expr<cell::triangle::fe::lagrange_p0>(h) *
					     (make_expr(b_0) * d<1>(u) + make_expr(b_1) * d<2>(u)) *
					     (make_expr(b_0) * d<1>(v) + make_expr(b_1) * d<2>(v))
					     , m);
    }
  }

  void assemble_linear_form() {
    f.clear();

    const auto v(f.get_test_function());

    f += integrate<volume_quadrature_type>(make_expr(src) * v
					   , m);

    if (supg_stabilisation) {
      const double delta(1.0);
      const double b_norm(1.0);
      
      finite_element_space<cell::triangle::fe::lagrange_p0> p0_fes(m);
      const auto h(build_element_diameter_function(m, p0_fes));

      f += integrate<volume_quadrature_type>(delta / b_norm * make_expr<cell::triangle::fe::lagrange_p0>(h) *
					     make_expr(src) *
					     (make_expr(b_0) * d<1>(v) + make_expr(b_1) * d<2>(v))
					     , m);
    }
  }
};

#endif /* _STEADY_ADVECTION_DIFFUSION_2D_H_ */
