#ifndef _UNSTEADY_ADVECTION_DIFFUSION_3D_H_
#define _UNSTEADY_ADVECTION_DIFFUSION_3D_H_

#include "../core/basic_fe_formulation.hpp"


/*
 *  The advection field b is supposed to be in H(div).
 */
template<typename fe>
class unsteady_advection_diffusion: public basic_fe_formulation<fe> {
public:
  using fe_type = typename basic_fe_formulation<fe>::fe_type;
  using fes_type = typename basic_fe_formulation<fe>::fes_type;
  using cell_type = typename basic_fe_formulation<fe>::cell_type;
  using element_type = typename basic_fe_formulation<fe>::element_type;
  using v_quad_type = typename basic_fe_formulation<fe>::volume_quadrature_type;
  using s_quad_type = typename basic_fe_formulation<fe>::boundary_quadrature_type;
  using p0_fe_type = finite_element::tetrahedron_lagrange_p0;
  using p0_element_type = finite_element_space<p0_fe_type>::element;

  unsteady_advection_diffusion(const mesh<cell_type>& m,
			       const submesh<cell_type>& dm,
			       double delta_t, double diffusivity, double supg_delta)
    : delta_t(delta_t), diffusivity(diffusivity), supg_delta(supg_delta),
      m(m), dm(dm), fes(m, dm), p0_fes(m),
      a(fes, fes), f(fes),
      solution(fes),
      bk_norm(p0_fes), h(build_element_diameter_function<cell_type>(m, p0_fes)),
      b_0(null_function), b_1(null_function), b_2(null_function),
      src(null_function) {}

  unsteady_advection_diffusion(const mesh<cell_type>& m, double delta_t, double diffusivity)
    : unsteady_advection_diffusion(m, submesh<cell_type>(m), delta_t, diffusivity) {}

  void set_boundary_value(const std::function<double(const double*)>& u_bc) {
    fes.set_dirichlet_condition(u_bc);
  }

  void set_initial_condition(const std::function<double(const double*)>& ic) {
    solution = projector::lagrange<fe_type>(ic, fes);
  }
  
  void set_advection_velocity(const std::function<double(const double*)>& b_0,
			      const std::function<double(const double*)>& b_1,
			      const std::function<double(const double*)>& b_2) {
    this->b_0 = b_0;
    this->b_1 = b_1;
    this->b_1 = b_2;

    const auto bk_0 = projector::l2<p0_fe_type, v_quad_type>(b_0, p0_fes);
    const auto bk_1 = projector::l2<p0_fe_type, v_quad_type>(b_1, p0_fes);
    const auto bk_2 = projector::l2<p0_fe_type, v_quad_type>(b_2, p0_fes);

    bk_norm = projector::l2<p0_fe_type,
			    v_quad_type>(compose(std::sqrt,
						 make_expr<p0_fe_type>(bk_0) * make_expr<p0_fe_type>(bk_0) + 
						 make_expr<p0_fe_type>(bk_1) * make_expr<p0_fe_type>(bk_1) +
						 make_expr<p0_fe_type>(bk_2) * make_expr<p0_fe_type>(bk_2)
						 )
					 , p0_fes);

    exporter::ensight6("supg_stab",
		       , bk_0, "bk_0"
		       , bk_1, "bk_1"
		       , bk_2, "bk_2"
		       , bk_norm, "bk_norm"
		       , h, "h");
    
    assemble_bilinear_form();
  }

  void set_source_term(const std::function<double(const double*)>& src) {
    this->src = src;
  }
  
  void step() {
    assemble_linear_form();
    
    dictionary p(dictionary()
                 .set("maxits",  2000u)
                 .set("restart", 1000u)
                 .set("rtol",    1.e-8)
                 .set("abstol",  1.e-50)
                 .set("dtol",    1.e20)
                 .set("ilufill", 2u));
    solver::petsc::gmres_ilu s(p);
    
    solution = a.solve(f, s);
  }

  element_type get_solution() const {
    return solution;
  }
      
private:
  const double delta_t;
  const double diffusivity;
  const double supg_delta;

  
  const mesh<cell_type>& m;
  const submesh<cell_type> dm;
  
  fes_type fes;
  finite_element_space<p0_fe_type> p0_fes;

  bilinear_form<fes_type, fes_type> a;
  linear_form<fes_type> f;

  element_type solution;
  p0_element_type  bk_norm, h;
  
  std::function<double(const double*)> b_0, b_1, b_2, src;

  bool assemble_time_derivative = true;
  bool assemble_advection_diffusion = true;
  bool assemble_source = false;
  bool use_diffusion_stabilisation = false;
  bool use_supg_stabilisation = true;

private:
  static double null_function(const double* x) { return 0.0; }
  static double inv(double x) { return 1.0 / x; }
  
  void assemble_bilinear_form() {
    a.clear();
    
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());


    if (assemble_time_derivative)
      a += integrate<v_quad_type>((1.0 / delta_t) * u * v, m);
    
    if (assemble_advection_diffusion)
      a += integrate<v_quad_type>(diffusivity * (d<1>(u) * d<1>(v) + d<2>(u) * d<2>(v) + d<3>(u) * d<3>(v)) +
				  make_expr(b_0) * d<1>(u) * v +
				  make_expr(b_1) * d<2>(u) * v +
				  make_expr(b_2) * d<3>(u) * v
				  , m);
    
    if (use_diffusion_stabilisation) {
      a += integrate<v_quad_type>(make_expr<p0_fe_type>(h) *
				  make_expr<p0_fe_type>(bk_norm) * 
				  (d<1>(u) * d<1>(v) + d<2>(u) * d<2>(v) + d<3>(u) * d<3>(v))
				  , m);
    }

    if (use_supg_stabilisation) {
      if (not std::is_same<fe_type, finite_element::tetrahedron_lagrange_p1>::value)
	throw std::string("unsteady advection diffusion 3d:"
			  " supg stabilisation is not available for non piece wise linear finite elements");
      
      a += integrate<v_quad_type>(//delta_t *
				  supg_delta *
				  compose(inv, make_expr<p0_fe_type>(bk_norm)) *
				  make_expr<p0_fe_type>(h) *
				  ( (1.0 / delta_t) * u +
				    make_expr(b_0) * d<1>(u) +
				    make_expr(b_1) * d<2>(u) +
				    make_expr(b_2) * d<3>(u) ) *
				  ( make_expr(b_0) * d<1>(v) +
				    make_expr(b_1) * d<2>(v) +
				    make_expr(b_2) * d<3>(v) )
				  , m);
    }
  }

  void assemble_linear_form() {
    f.clear();

    const auto v(f.get_test_function());

    
    if (assemble_time_derivative)
      f += integrate<v_quad_type>((1.0 / delta_t) * make_expr<fe_type>(solution) * v, m);

    if (assemble_source)
      f += integrate<v_quad_type>(make_expr(src) * v, m);

    if (supg_stabilisation) {
      f += integrate<v_quad_type>(//delta_t *
				  supg_delta *
				  compose(inv, make_expr<p0_fe_type>(bk_norm)) *
				  make_expr<p0_fe_type>(h) *
				  ( (1.0 / delta_t) * make_expr<fe_type>(solution) +
				    make_expr(src) ) *
				  ( make_expr(b_0) * d<1>(v) +
				    make_expr(b_1) * d<2>(v) +
				    make_expr(b_2) * d<3>(v) )
				  , m);
    }
  }
};

#endif /* _UNSTEADY_ADVECTION_DIFFUSION_3D_H_ */
