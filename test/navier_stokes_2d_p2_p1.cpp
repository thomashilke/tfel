
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


int main(int argc, char* argv[]) {
  try {
    if (argc != 2)
      throw std::string("wrong number of arguments. "
                        "An single integer describing the "
                        "number of subdivisions is expected.");
    
    const std::size_t n(std::stoi(argv[1]));
    
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

    //array<double> tmp{fes.get_dof_number()};
    //tmp.fill(0.0);
    /*
     *  velocity at successive timesteps (u_n and u_{n+1})
     */
    array<double> u_comp{fes.get_dof_number<0>()};
    array<double> p_comp{fes.get_dof_number<2>()};
    u_comp.fill(0.0);
    p_comp.fill(0.0);
    finite_element_space<u_fe_type>::element u0(fes.get_finite_element_space<0>(), u_comp),
      u1(fes.get_finite_element_space<1>(), u_comp);
    finite_element_space<p_fe_type>::element p(fes.get_finite_element_space<2>(), p_comp);
    
    fes_type::element u(fes, u0, u1, p);
    fes_type::element up(u);


    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto v0(a.get_test_function<0>());
      auto v1(a.get_test_function<1>());
      auto q (a.get_test_function<2>());
      
      auto u0(a.get_trial_function<0>());
      auto u1(a.get_trial_function<1>());
      auto p (a.get_trial_function<2>());
      
      a += integrate<quad_type>(    d<1>(u0) * d<1>(v0) + d<2>(u0) * d<2>(v0)
                                  + d<1>(u1) * d<1>(v1) + d<2>(u1) * d<2>(v1)
                                  + p * (d<1>(v0) + d<2>(v1)) 
                                  + q * (d<1>(u0) + d<2>(u1))
                                , m);
    }

    linear_form<fes_type> f(fes); {
      auto v0(f.get_test_function<0>());
      auto v1(f.get_test_function<1>());
      
      f += integrate<quad_type>(    f0 * v0
                                  + f1 * v1
                                , m);
    }

    const fes_type::element x(a.solve(f));

    {
    const auto u0(x.get_component<0>());
    const auto u1(x.get_component<1>());
    const auto p (x.get_component<2>());

    const double u0_error(std::sqrt(integrate<quad_type>((make_expr<u_fe_type>(u0) - make_expr(u0_bv)) *
                                                         (make_expr<u_fe_type>(u0) - make_expr(u0_bv))
                                                         , m)));
    const double u1_error(std::sqrt(integrate<quad_type>((make_expr<u_fe_type>(u1) - make_expr(u1_bv)) *
                                                         (make_expr<u_fe_type>(u1) - make_expr(u1_bv))
                                                         , m)));
    const double p_error(std::sqrt(integrate<quad_type>((make_expr<p_fe_type>(p)) *
                                                        (make_expr<p_fe_type>(p))
                                                        , m)));
    const double u_error(std::sqrt(u0_error * u0_error + u1_error * u1_error));

    std::cout << "u L2 error: " << u_error << std::endl
              << "p L2 error: " << p_error << std::endl;
    
    exporter::ensight6("stokes",
                       to_mesh_vertex_data<u_fe_type>(u0), "u0",
                       to_mesh_vertex_data<u_fe_type>(u1), "u1",
                       to_mesh_vertex_data<p_fe_type>(p),  "p");
    }
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  
  return 0;
}
