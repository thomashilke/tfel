
#include "../src/core/composite_fe.hpp"
#include "../src/core/composite_fes.hpp"
#include "../src/core/composite_form.hpp"
#include "../src/core/quadrature.hpp"
#include "../src/core/export.hpp"

/*
 *  Kim & Moin solution of the Stokes formulation on a unit square
 */
/*
double f0(const double* x) {
  return - 2.0 * std::cos((x[0] - 0.5) * M_PI) * std::sin((x[1] - 0.5) * M_PI) * M_PI * M_PI;
}

double f1(const double* x) {
  return   2.0 * std::sin((x[0] - 0.5) * M_PI) * std::cos((x[1] - 0.5) * M_PI) * M_PI * M_PI;
}

double u0_bv(const double* x) {
  return - std::cos((x[0] - 0.5) * M_PI) * std::sin((x[1] - 0.5) * M_PI);
}

double u1_bv(const double* x) {
  return std::sin((x[0] - 0.5) * M_PI) * std::cos((x[1] - 0.5) * M_PI);
}

double p_v(const double* x) {
  return 0.0;
}
*/

/*
 *  Two sided driven cavity on a cube
 */
double f0(const double* x) {
  return 0.0;
}

double f1(const double* x) {
  return 0.0;
}

double u0_bv(const double* x) {
  return x[1] < 0.0001 ? x[0] * (1.0 - x[0]) : 0.0;
}

double u1_bv(const double* x) {
  return x[0] < 0.0001 ? x[1] * (1.0 - x[1]) : 0.0;
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
    //using u_fe_type = cell_type::fe::lagrange_p1_bubble;
    using p_fe_type = cell_type::fe::lagrange_p1;
    using quad_type = quad::triangle::qf5pT;

    fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
    submesh<cell_type> dm(m.get_boundary_submesh());
    submesh<cell_type, cell::point> pinned_pressure_point(m.get_point_submesh(n * n / 2 + n/2));

    using fe_type = composite_finite_element<u_fe_type, u_fe_type, p_fe_type>;
    using fes_type = composite_finite_element_space<fe_type>;
    fes_type fes(m);

    fes.add_dirichlet_boundary<0>(dm, u0_bv);
    fes.add_dirichlet_boundary<1>(dm, u1_bv);
    fes.add_dirichlet_boundary<2>(pinned_pressure_point, p_v);
    
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
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  
  return 0;
}
