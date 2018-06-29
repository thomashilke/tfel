
#include "../src/tfel.hpp"

double src(const double* x) {
  return (x[0] - 0.5) + (x[1] - 0.5);
}

int main(int argc, char *argv[]) {
  const std::size_t n(100);

  using cell_type = cell::triangle;
  using quad_type = quad::triangle::qf5pT;
  using fe_type = cell::triangle::fe::lagrange_p1;
  using fes_type = finite_element_space<fe_type>;

  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type, cell::point> dm(m.get_point_submesh(0));
  fes_type fes(m);

  bilinear_form<fes_type, fes_type> a(fes, fes, 1, 1); {
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());

    a += integrate<quad_type>(d<1>(u) * d<1>(v) +
			      d<2>(u) * d<2>(v)
			      , m);

    a.algebraic_trial_block(0) += integrate<quad_type>(u, m);
    a.algebraic_test_block(0) += integrate<quad_type>(v, m);
    a.algebraic_block(0, 0) = 0.0;
  }

  linear_form<fes_type> f(fes, 1); {
    const auto v(f.get_test_function());

    f += integrate<quad_type>(src * v
			      , m);
    f.algebraic_equation_value(0) = 1.0;
  }

  dictionary p(dictionary()
               .set("maxits",  2000u)
               .set("restart", 1000u)
               .set("rtol",    1.e-8)
               .set("abstol",  1.e-50)
               .set("dtol",    1.e20)
               .set("ilufill", 2u));
  solver::petsc::gmres_ilu s(p);

  exporter::ensight6("poisson_constraint"
		     , to_mesh_vertex_data<fe_type>(a.solve(f, s)), "solution");
  
  return 0;
}
