#ifndef _STEADY_ADVECTION_DIFFUSION_1D_H_
#define _STEADY_ADVECTION_DIFFUSION_1D_H_

#include "../core/basic_fe_formulation.hpp"


/*
 *  The advection field b is supposed to be in H(div), that is constant in 1d.
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


  steady_advection_diffusion(const mesh<cell_type>& m,
			     const submesh<cell_type>& dm,
			     double diffusivity)
    : diffusivity(diffusivity),
      m(m), dm(dm),
      fes(m, dm), p0_fes(m),
      a(fes, fes), f(fes),
      solution(fes),
      h(build_element_diameter_function(m, p0_fes)), bk_norm(p0_fes),
      b(null_function), src(null_function) {

    static_assert(cell_type::n_dimension == 1, "This is a 1d formulation");
  }

  steady_advection_diffusion(const mesh<cell_type>& m,
			     double diffusivity)
    : steady_advection_diffusion(m, m.get_boundary_submesh(), diffusivity) {}

  void set_boundary_value(const std::function<double(const double*)>& u_bc) {
    fes.set_dirichlet_condition(u_bc);
  }
  
  void set_advection_velocity(const std::function<double(const double*)>& b) {
    this->b = b;

    const auto bk = projector::l2<cell::edge::fe::lagrange_p0,
				  volume_quadrature_type>(b, p0_fes);

    using p0_fe = cell::edge::fe::lagrange_p0;
    bk_norm = projector::l2<cell::edge::fe::lagrange_p0,
			    volume_quadrature_type>(compose(std::abs, make_expr<p0_fe>(bk)), p0_fes);
    
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

  const mesh<cell_type>& m;
  const submesh<cell_type> dm;

  fes_type fes;
  finite_element_space<cell::edge::fe::lagrange_p0> p0_fes;

  bilinear_form<fes_type, fes_type> a;
  linear_form<fes_type> f;

  element_type solution;
  finite_element_space<cell::edge::fe::lagrange_p0>::element h;
  finite_element_space<cell::edge::fe::lagrange_p0>::element bk_norm;

  std::function<double(const double*)> b, src;

  double diffusion_stabilisation = false;
  double supg_stabilisation = false;

private:
  static double null_function(const double* x) { return 0.0; }
  static double inv(double x) { return 1.0 / x; }
  
  void assemble_bilinear_form() {
    a.clear();
    
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());

    a += integrate<volume_quadrature_type>(diffusivity * d<1>(u) * d<1>(v) +
					   make_expr(b) * d<1>(u) * v
					   , m);
    if (diffusion_stabilisation) {

      const double h(m.get_h_max());
      const double b_norm(1.0);
      a += integrate<volume_quadrature_type>(h * b_norm * d<1>(u) * d<1>(v)
					     , m);
    }

    if (supg_stabilisation) {
      if (not std::is_same<fe_type, cell::edge::fe::lagrange_p1>::value)
	throw std::string("steady advection diffusion 1d:"
			  " supg stabilisation is not available for non piece wise linear finite elements");
      
      const double delta(0.5);
      a += integrate<volume_quadrature_type>(delta *
					     compose(inv, make_expr<cell::edge::fe::lagrange_p0>(bk_norm)) *
					     make_expr<cell::edge::fe::lagrange_p0>(h) *
					     (make_expr(b) * d<1>(u)) *
					     (make_expr(b) * d<1>(v))
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
      
      finite_element_space<cell::edge::fe::lagrange_p0> p0_fes(m);
      const auto h(build_element_diameter_function(m, p0_fes));

      f += integrate<volume_quadrature_type>(delta *
					     compose(inv, make_expr<cell::edge::fe::lagrange_p0>(bk_norm)) *
					     make_expr<cell::edge::fe::lagrange_p0>(h) *
					     make_expr(src) *
					     (make_expr(b) * d<1>(v))
					     , m);
    }
  }
};


#endif /* _STEADY_ADVECTION_DIFFUSION_1D_H_ */
