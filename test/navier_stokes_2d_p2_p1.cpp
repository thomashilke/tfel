
#include <spikes/timer.hpp>

#include "../src/core/composite_fe.hpp"
#include "../src/core/composite_fes.hpp"
#include "../src/core/composite_form.hpp"
#include "../src/core/quadrature.hpp"
#include "../src/core/export.hpp"

/*
 *  Regularized driven south and east boundary condition for u
 */

double f0(const double* x) {
  return 0.0;
}

double f1(const double* x) {
  return 0.0;
}

double u0_bv(const double* x) {
  return x[1] < 0.00001 ? 4.0 * x[0] * (1.0 - x[0]) : 0.0;
}

double u1_bv(const double* x) {
  return x[0] < 0.00001 ? 4.0 * x[1] * (1.0 - x[1]) : 0.0;
}

double p_v(const double* x) {
  return 0.0;
}


void benchmark(std::size_t n) {
  timer t;

  double
    mesh_timestamp(0.0),
    fes_timestamp(0.0),
    ic_timestamp(0.0),
    lin_timestamp(0.0),
    bilin_timestamp(0.0);

  
  std::cout << std::setw(10) << std::right << t.tic()
            << "  benchmark start" << std::endl;
  
  using cell_type = cell::triangle;
  using u_fe_type = cell_type::fe::lagrange_p2;
  using p_fe_type = cell_type::fe::lagrange_p1;
  using quad_type = quad::triangle::qf5pT;

  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());
  submesh<cell_type, cell::point> pinned_pressure_point(m.get_point_submesh(n * n / 2 + n/2));

  mesh_timestamp = t.tic();
  std::cout << std::setw(10) << std::right << mesh_timestamp
            << "  mesh, boundary submesh and pressure pin point" << std::endl;

  using fe_type = composite_finite_element<u_fe_type, u_fe_type, p_fe_type>;
  using fes_type = composite_finite_element_space<fe_type>;
  fes_type fes(m);

  fes.set_dirichlet_boundary_condition<0>(dm, u0_bv);
  fes.set_dirichlet_boundary_condition<1>(dm, u1_bv);
  fes.set_dirichlet_boundary_condition<2>(pinned_pressure_point, p_v);

  fes_timestamp = t.tic();
  std::cout << std::setw(10) << std::right << fes_timestamp
            << "  finite element space setup" << std::endl;

  /*
   *  velocity at successive timesteps (u_n and u_{n+1})
   */
  array<double> u_comp{fes.get_dof_number<0>()};
  array<double> p_comp{fes.get_dof_number<2>()};
  u_comp.fill(0.0);
  p_comp.fill(0.0);
  finite_element_space<u_fe_type>::element
    u0(fes.get_finite_element_space<0>(), u_comp),
    u1(fes.get_finite_element_space<1>(), u_comp);
  finite_element_space<p_fe_type>::element
    p(fes.get_finite_element_space<2>(), p_comp);
    
  fes_type::element x(fes, u0, u1, p);
  fes_type::element xp(x);

  ic_timestamp = t.tic();
  std::cout << std::setw(10) << std::right << ic_timestamp
            << "  initial conditions" << std::endl;

  /*
   *  model and numerical parameters
   */
  const double reynolds(5000.0);
  const double time_step(0.5);

  auto v0_h(xp.template get_component<0>());
  auto v1_h(xp.template get_component<1>());

  auto v0(make_expr<u_fe_type>(v0_h));
  auto v1(make_expr<u_fe_type>(v1_h)); 

        
  bilinear_form<fes_type, fes_type> a(fes, fes); {
    auto w0(a.get_test_function<0>());
    auto w1(a.get_test_function<1>());
    auto q (a.get_test_function<2>());
      
    auto vn0(a.get_trial_function<0>());
    auto vn1(a.get_trial_function<1>());
    auto pn (a.get_trial_function<2>());
    
    a += integrate<quad_type>(vn0 * w0 + vn1 * w1 +
          
                              time_step * (vn0 * d<1>(v0) * w0 +
                                           vn0 * d<1>(v1) * w1 +
                                           vn1 * d<2>(v0) * w0 +
                                           vn1 * d<2>(v1) * w1) +

                              time_step * (v0 * d<1>(vn0) * w0 +
                                           v0 * d<1>(vn1) * w1 +
                                           v1 * d<2>(vn0) * w0 +
                                           v1 * d<2>(vn1) * w1) +

                              (time_step / reynolds) * (d<1>(vn0) * d<1>(w0) + d<2>(vn0) * d<2>(w0) +
                                                        d<1>(vn1) * d<1>(w1) + d<2>(vn1) * d<2>(w1)) +
                                    
                              time_step * pn * (d<1>(w0) + d<2>(w1)) +
                                    
                              time_step * q * (d<1>(vn0) + d<2>(vn1))
                              , m);
    /*
    // \int v^{k+1}w
    a += integrate<quad_type>(vn0 * w0 + vn1 * w1
    , m);
          
    // dt \int v^{k+1}\nabla v^k w
    a += integrate<quad_type>( time_step * (vn0 * d<1>(v0) * w0 +
    vn0 * d<1>(v1) * w1 +
    vn1 * d<2>(v0) * w0 +
    vn1 * d<2>(v1) * w1)
    , m);

    // dt \int v^k \nabla v^{k+1} w
    a += integrate<quad_type>( time_step * (v0 * d<1>(vn0) * w0 +
    v0 * d<1>(vn1) * w1 +
    v1 * d<2>(vn0) * w0 +
    v1 * d<2>(vn1) * w1)
    , m);

    // \frac{dt}{Re} \int \nabla v^{k+1}\nabla w
    a += integrate<quad_type>( (time_step / reynolds) * (d<1>(vn0) * d<1>(w0) + d<2>(vn0) * d<2>(w0) +
    d<1>(vn1) * d<1>(w1) + d<2>(vn1) * d<2>(w1))
    , m);

    // dt \int p \div w
    a += integrate<quad_type>(time_step * pn * (d<1>(w0) + d<2>(w1))
    , m);

    // dt \int q \div v^{k+1}
    a += integrate<quad_type>(time_step * q * (d<1>(vn0) + d<2>(vn1))
    , m);
    */
  }
  bilin_timestamp = t.tic();
  std::cout << std::setw(10) << std::right << bilin_timestamp
            << "  bilinear form ("
            << (bilin_timestamp - ic_timestamp) / m.get_cell_number()
            << " millisec/cell)" << std::endl;

  linear_form<fes_type> f(fes); {
    auto w0(f.get_test_function<0>());
    auto w1(f.get_test_function<1>());

    auto u0_h(x.template get_component<0>());
    auto u1_h(x.template get_component<1>());

    auto u0(make_expr<u_fe_type>(u0_h));
    auto u1(make_expr<u_fe_type>(u1_h));

    f += integrate<quad_type>( u0 * w0 + u1 * w1 +

                               time_step * (f0 * w0 +
                                            f1 * w1) + 

                               time_step * (v0 * d<1>(v0) * w0 +
                                            v0 * d<1>(v1) * w1 +
                                            v1 * d<2>(v0) * w0 +
                                            v1 * d<2>(v1) * w1)
                               , m);
    /*
      f += integrate<quad_type>( u0 * w0 + u1 * w1
      , m);
      
      f += integrate<quad_type>( time_step * (f0 * w0 +
      f1 * w1)
      , m);

      f += integrate<quad_type>( time_step * (v0 * d<1>(v0) * w0 +
      v0 * d<1>(v1) * w1 +
      v1 * d<2>(v0) * w0 +
      v1 * d<2>(v1) * w1)
      , m);
    */
  }
  lin_timestamp = t.tic();
  std::cout << std::setw(10) << std::right << lin_timestamp
            << "  linear form ("
            << (lin_timestamp - bilin_timestamp) / m.get_cell_number()
            << " millisec/cell)" << std::endl;
}

void navier_stokes(std::size_t n);

int main(int argc, char* argv[]) {
  try {
    if (argc != 2)
      throw std::string("wrong number of arguments. "
                        "An single integer describing the "
                        "number of subdivisions is expected.");
    
    const std::size_t n(std::stoi(argv[1]));
    //navier_stokes(n);
    benchmark(n);
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  
  return 0;
}

void navier_stokes(std::size_t n) {
  using cell_type = cell::triangle;
  using u_fe_type = cell_type::fe::lagrange_p2;
  using p_fe_type = cell_type::fe::lagrange_p1;
  using quad_type = quad::triangle::qf5pT;

  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());
  submesh<cell_type, cell::point> pinned_pressure_point(m.get_point_submesh(n * n / 2 + n/2));

  using fe_type = composite_finite_element<u_fe_type, u_fe_type, p_fe_type>;
  using fes_type = composite_finite_element_space<fe_type>;
  fes_type fes(m);

  fes.set_dirichlet_boundary_condition<0>(dm, u0_bv);
  fes.set_dirichlet_boundary_condition<1>(dm, u1_bv);
  fes.set_dirichlet_boundary_condition<2>(pinned_pressure_point, p_v);

  /*
   *  velocity at successive timesteps (u_n and u_{n+1})
   */
  array<double> u_comp{fes.get_dof_number<0>()};
  array<double> p_comp{fes.get_dof_number<2>()};
  u_comp.fill(0.0);
  p_comp.fill(0.0);
  finite_element_space<u_fe_type>::element
    u0(fes.get_finite_element_space<0>(), u_comp),
    u1(fes.get_finite_element_space<1>(), u_comp);
  finite_element_space<p_fe_type>::element
    p(fes.get_finite_element_space<2>(), p_comp);
    
  fes_type::element x(fes, u0, u1, p);
  fes_type::element xp(x);

  /*
   *  model and numerical parameters
   */
  const double reynolds(5000.0);
  const double end_time(20.0);
  const double time_step(0.5);
  const std::size_t end_time_step(std::ceil(end_time / time_step));
  const std::size_t newton_max_step(10);
  const double newton_rtol(1.0e-8);

  /*
   *  storage for timing statistics 
   */
  double
    bilinear_form_elapsed_time(0.0),
    linear_form_elapsed_time(0.0),
    linear_solve_elapsed_time(0.0);

  /*
   *  time iteration loop
   */
  double time(0.0);
  for (std::size_t k(0); k < end_time_step; ++k) {
    std::cout << "time step #" << k << "(time " << time << ")" << std::endl;
      
    time += time_step;


    /*
     *  newton iteration loop
     */
    bool newton_done(false);
    for (std::size_t j(0); j < newton_max_step and not newton_done; ++j) {
      std::cout << "newton step #" << j << " of " << newton_max_step << std::endl;

      auto v0_h(xp.template get_component<0>());
      auto v1_h(xp.template get_component<1>());

      auto v0(make_expr<u_fe_type>(v0_h));
      auto v1(make_expr<u_fe_type>(v1_h)); 

        
      bilinear_form<fes_type, fes_type> a(fes, fes); {
        timer t;
          
        auto w0(a.get_test_function<0>());
        auto w1(a.get_test_function<1>());
        auto q (a.get_test_function<2>());
      
        auto vn0(a.get_trial_function<0>());
        auto vn1(a.get_trial_function<1>());
        auto pn (a.get_trial_function<2>());

        a += integrate<quad_type>(vn0 * w0 + vn1 * w1 +
          
                                  time_step * (vn0 * d<1>(v0) * w0 +
                                               vn0 * d<1>(v1) * w1 +
                                               vn1 * d<2>(v0) * w0 +
                                               vn1 * d<2>(v1) * w1) +

                                  time_step * (v0 * d<1>(vn0) * w0 +
                                               v0 * d<1>(vn1) * w1 +
                                               v1 * d<2>(vn0) * w0 +
                                               v1 * d<2>(vn1) * w1) +

                                  (time_step / reynolds) * (d<1>(vn0) * d<1>(w0) + d<2>(vn0) * d<2>(w0) +
                                                            d<1>(vn1) * d<1>(w1) + d<2>(vn1) * d<2>(w1)) +
                                    
                                  time_step * pn * (d<1>(w0) + d<2>(w1)) +
                                    
                                  time_step * q * (d<1>(vn0) + d<2>(vn1))
                                  , m);

        bilinear_form_elapsed_time = t.tic();
      }

      linear_form<fes_type> f(fes); {
        timer t;
          
        auto w0(f.get_test_function<0>());
        auto w1(f.get_test_function<1>());

        auto u0_h(x.template get_component<0>());
        auto u1_h(x.template get_component<1>());

        auto u0(make_expr<u_fe_type>(u0_h));
        auto u1(make_expr<u_fe_type>(u1_h));

        f += integrate<quad_type>( u0 * w0 + u1 * w1 +

                                   time_step * (f0 * w0 +
                                                f1 * w1) + 

                                   time_step * (v0 * d<1>(v0) * w0 +
                                                v0 * d<1>(v1) * w1 +
                                                v1 * d<2>(v0) * w0 +
                                                v1 * d<2>(v1) * w1)
                                   , m);
        /*
          f += integrate<quad_type>( u0 * w0 + u1 * w1
          , m);
      
          f += integrate<quad_type>( time_step * (f0 * w0 +
          f1 * w1)
          , m);

          f += integrate<quad_type>( time_step * (v0 * d<1>(v0) * w0 +
          v0 * d<1>(v1) * w1 +
          v1 * d<2>(v0) * w0 +
          v1 * d<2>(v1) * w1)
          , m);
        */
        linear_form_elapsed_time = t.tic();
      }

      timer t;
      auto xpp(a.solve(f));
      linear_solve_elapsed_time = t.tic();

        
      /*
       *  compute the velocity step norm, velocity norm and relative step norm
       */
      auto vn0_h(xpp.template get_component<0>());
      auto vn1_h(xpp.template get_component<1>());

      const double
        newton_step_0_norm(std::sqrt(integrate<quad_type>((make_expr<u_fe_type>(vn0_h) - make_expr<u_fe_type>(v0_h)) *
                                                          (make_expr<u_fe_type>(vn0_h) - make_expr<u_fe_type>(v0_h))
                                                          , m))),
        newton_step_1_norm(std::sqrt(integrate<quad_type>((make_expr<u_fe_type>(vn1_h) - make_expr<u_fe_type>(v1_h)) *
                                                          (make_expr<u_fe_type>(vn1_h) - make_expr<u_fe_type>(v1_h))
                                                          , m)));
      const double
        velocity_0_norm(std::sqrt(integrate<quad_type>(make_expr<u_fe_type>(vn0_h) * make_expr<u_fe_type>(vn0_h), m))),
        velocity_1_norm(std::sqrt(integrate<quad_type>(make_expr<u_fe_type>(vn1_h) * make_expr<u_fe_type>(vn1_h), m)));

      std::cout << "newton step norm (" << newton_step_0_norm << ", " << newton_step_1_norm << ")" << std::endl;
      std::cout << "velocity norm (" << velocity_0_norm << ", " << velocity_1_norm << ")" << std::endl;
      std::cout << "newton relative step norm ("
                << newton_step_0_norm / velocity_0_norm << ", "
                << newton_step_1_norm / velocity_1_norm << ")"
                << std::endl;

      const double
        newton_relative_step_0_norm(newton_step_0_norm / velocity_0_norm),
        newton_relative_step_1_norm(newton_step_1_norm / velocity_1_norm);
        
      newton_done = (newton_relative_step_0_norm < newton_rtol and
                     newton_relative_step_1_norm < newton_rtol);
        
      xp = xpp;
      std::cout << "elapsed time (bilinear, linear, solve): ("
                << bilinear_form_elapsed_time << ", "
                << linear_form_elapsed_time << ", "
                << linear_solve_elapsed_time << ")" << std::endl;
    }

    x = xp;
  }

  /*
   *  export the solution at final time
   */
  auto u0_h(xp.template get_component<0>());
  auto u1_h(xp.template get_component<1>());
  exporter::ensight6("navier_stokes_2d_p2_p1",
                     to_mesh_vertex_data<u_fe_type>(u0_h), "u0",
                     to_mesh_vertex_data<u_fe_type>(u1_h), "u1");
}
