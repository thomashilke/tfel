#ifndef _FRM_STOKES_2D_H_
#define _FRM_STOKES_2D_H_

#include "../core/basic_fe_formulation.hpp"
#include "../core/element_diameter.hpp"


/*
 *  Stokes problem weak formulation
 *
 *  This formulation is valid at least for the following
 *  finite element choices on triangles:
 *   - p1_bubble - p1,
 *   - p1 - p1.
 *  Proper stabilisation term is added in the latter case.
 *  (See Quarteroni-Valli, 2008, Chap. 9.4, eq. (9.4.2), pp. 312)
 *
 *  The value of the pressure is pinned in one point (the
 *  first node of the first element of the mesh)
 *
 *  The dirichlet boundary is taken as the entire boundary
 *  of the domain m.
 */


template<typename v_fe_type = cell::triangle::fe::lagrange_p1_bubble,
	 typename p_fe_type = cell::triangle::fe::lagrange_p1>
class stokes_2d
  : public basic_fe_formulation<composite_finite_element<v_fe_type,
							 v_fe_type,
							 p_fe_type> > {
public:
  using stokes_2d_fe = composite_finite_element<v_fe_type, v_fe_type, p_fe_type>;
  using typename basic_fe_formulation<stokes_2d_fe>::fe_type;
  using typename basic_fe_formulation<stokes_2d_fe>::fes_type;
  using typename basic_fe_formulation<stokes_2d_fe>::cell_type;
  using typename basic_fe_formulation<stokes_2d_fe>::element_type;
  using typename basic_fe_formulation<stokes_2d_fe>::volume_quadrature_type;
  using typename basic_fe_formulation<stokes_2d_fe>::boundary_quadrature_type;

  using velocity_fe_type = typename fe_type::template fe_type<0>;
  using pressure_fe_type = typename fe_type::template fe_type<2>;

  using velocity_fes_type = finite_element_space<velocity_fe_type>;
  using pressure_fes_type = finite_element_space<pressure_fe_type>;
  
  stokes_2d(const fe_mesh<cell_type>& m, double mu)
    : viscosity(mu),
      m(m), dm(m.get_boundary_submesh()),
      pressure_point_m(m.get_point_submesh(0)),
      fes(m),
      a(fes, fes), f(fes),
      solution(fes) {
    
    fes.template add_dirichlet_boundary<2>(pressure_point_m, 0.0);
    a.clear();
    
    // The dirichlet dof must be known when assembling the system.
    // The values are relevent only when solving, though.
    assemble_bilinear_form();
  }

  const element_type& solve() {
    dictionary p(dictionary()
                 .set("maxits",  2000u)
                 .set("restart", 1000u)
                 .set("rtol",    1.e-8)
                 .set("abstol",  1.e-50)
                 .set("dtol",    1.e20)
                 .set("ilufill", 2u));
    solver::petsc::gmres_ilu s(p);

    solution = a.solve(f, s);
    return solution;
  }
  
  const element_type& get_solution() const { return solution; }

  void set_velocity_condition(const std::function<double(const double*)>& u_0_bc,
			      const std::function<double(const double*)>& u_1_bc) {
    fes.template add_dirichlet_boundary<0>(dm, u_0_bc);
    fes.template add_dirichlet_boundary<1>(dm, u_1_bc);
  }

  void set_force(const std::function<double(const double*)>& f_0,
		 const std::function<double(const double*)>& f_1) {
    f.clear();
    
    const auto v_0(f.template get_test_function<0>());
    const auto v_1(f.template get_test_function<1>());

    f += integrate<quad::triangle::qf5pT>(make_expr(f_0) * v_0 +
					  make_expr(f_1) * v_1
					  , m);
  }

private:
  double viscosity;

  const fe_mesh<cell_type>& m;
  const submesh<cell_type> dm;
  const submesh<cell_type, cell::point> pressure_point_m;
  
  fes_type fes;

  bilinear_form<fes_type, fes_type> a;
  linear_form<fes_type> f;

  typename fes_type::element solution;

private:
  void assemble_bilinear_form() {
    a.clear();
    
    const auto u_0(a.template get_trial_function<0>());
    const auto u_1(a.template get_trial_function<1>());
    const auto p(a.template get_trial_function<2>());
    
    const auto v_0(a.template get_test_function<0>());
    const auto v_1(a.template get_test_function<1>());
    const auto q(a.template get_test_function<2>());
    
    // Stokes's weak formulation
    a += integrate<quad::triangle::qf5pT>(viscosity * (d<1>(u_0)*d<1>(v_0) + d<2>(u_0)*d<2>(v_0) +
						       d<1>(u_1)*d<1>(v_1) + d<2>(u_1)*d<2>(v_1)) + 
					  p * (d<1>(v_0) + d<2>(v_1)) + 
					  q * (d<1>(u_0) + d<2>(u_1))
					  , m);

    // Stabilisation for the P_1 - P_1 finite element choice
    if (not std::is_same<velocity_fe_type, cell::triangle::fe::lagrange_p1_bubble>::value) {
      const auto h(cell_diameter(m));
      const double stab_coefficient(- 1.0);

      a += integrate<quad::triangle::qf5pT>(stab_coefficient *
                                            make_expr(h, 0) * make_expr(h, 0) * 
					    (d<1>(p)*d<1>(q) + d<2>(p)*d<2>(q))
					    , m);
    }
  }
};

#endif /* _FRM_STOKES_2D_H_ */
