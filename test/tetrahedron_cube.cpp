
#include <iostream>

#include "../src/core/mesh.hpp"
#include "../src/core/export.hpp"
#include "../src/core/form.hpp"

double u_bc(const double* x) {
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

double f(const double* x) {
  return -6.0;
}

int main(int argc, char *argv[]) {
  try {
    using cell_type = cell::tetrahedron;
    using fe_type = cell::tetrahedron::fe::lagrange_p1_bubble;
    using fes_type = finite_element_space<fe_type>;
    using quad_type = quad::tetrahedron::qf5pTet;
    
    fe_mesh<cell_type> m(gen_cube_mesh(1.0, 1.0, 1.0, 10, 10, 10));
    //m.show(std::cout);

    submesh<cell_type> dm(m.get_boundary_submesh());
    //dm.show(std::cout);

    fes_type fes(m, dm);
    //fes.show(std::cout);
    
    fes.set_dirichlet_condition(u_bc);

    bilinear_form<fes_type, fes_type> a(fes, fes); {
      const auto u(a.get_trial_function());
      const auto v(a.get_test_function());

      a += integrate<quad_type>(
	d<1>(u) * d<1>(v) +
	d<2>(u) * d<2>(v) +
	d<3>(u) * d<3>(v)
	, m);

      //a.show(std::cout);
    }

    linear_form<fes_type> f(fes); {
      const auto v(f.get_test_function());

      f += integrate<quad_type>(
	::f * v,
	m);
    }

    std::cout << "\\int_\\Omega x_1^2 dx = "
              << integrate<quad::tetrahedron::qf5pTet>(make_expr(std::function<double(const double*)>([](const double* x){return x[1]; })), m) - 0.5
              << std::endl;
    
    exporter::ensight6("cube",
		       to_mesh_vertex_data<fe_type>(a.solve(f)), "solution");
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  
  return 0;
}
