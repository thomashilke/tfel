#ifndef _TRANSIENT_DIFFUSION_2D_H_
#define _TRANSIENT_DIFFUSION_2D_H_

#include "../basic_fe_formulation.hpp"


template<typename fe>
class transient_diffusion_2d: public basic_fe_formulation<fe> {
public:
  using typename basic_fe_formulation<fe>::fe_type;
  using typename basic_fe_formulation<fe>::fes_type;
  using typename basic_fe_formulation<fe>::cell_type;
  using typename basic_fe_formulation<fe>::element_type;
  using typename basic_fe_formulation<fe>::volume_quadrature_type;
  using typename basic_fe_formulation<fe>::boundary_quadrature_type;
  
  transient_diffusion_2d(const mesh<cell_type>& m,
			 double delta_t, double diffusivity)
    : delta_t(delta_t), diffusivity(diffusivity),
      m(m), dm(m.get_boundary_submesh()),
      fes(m, dm), a(fes, fes), f(fes),
      source(fes), solution(fes) {
    assemble_bilinear_form();
  }

  transient_diffusion_2d(const mesh<cell_type>& m,
			 const submesh<cell_type>& dm,
			 double delta_t, double diffusivity)
    : delta_t(delta_t), diffusivity(diffusivity),
      m(m), dm(dm),
      fes(m, dm), a(fes, fes), f(fes),
      source(fes), solution(fes) {
    assemble_bilinear_form();
  }

  void set_boundary_value(const std::function<double(const double*)>& u_bc) {
    fes.set_dirichlet_condition(u_bc);
  }

  void set_initial_condition(double (*ic)(const double* )) {
    solution = projector::l2<fe_type, volume_quadrature_type>(ic, fes);
  }

  void set_source(const std::function<double(const double*)>& src) {
    source = projector::l2<fe_type, volume_quadrature_type>(src, fes);
  }

  void step() {
    assemble_linear_form();
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
  void assemble_bilinear_form() {
    a.clear();
    
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());

    if (assemble_time_derivative)
      a += integrate<volume_quadrature_type>((1.0 / delta_t) * u * v, m);

    if (assemble_laplacian)
      a += integrate<volume_quadrature_type>(diffusivity * (d<1>(u) * d<1>(v)
							    + d<2>(u) * d<2>(v)),
					     m);
  }

  void assemble_linear_form() {
    f.clear();

    const auto v(f.get_test_function());

    if (assemble_time_derivative)
      f += integrate<volume_quadrature_type>((1.0 / delta_t) * make_expr<fe_type>(solution) * v, m);

    if (assemble_source)
      f += integrate<volume_quadrature_type>(make_expr<fe_type>(source) * v, m);
  }
};


#endif /* _TRANSIENT_DIFFUSION_2D_H_ */
