
#include <tfel/tfel.hpp>

double gaussian(double x, double sigma) {
  return std::exp(- x * x / (2.0 * sigma * sigma));
}

double src(const double* x) {
  const double sigma(0.1);
  const double center[] = {0.5, 0.5};
  const double r(std::sqrt(std::pow(x[0] - center[0], 2) +
			   std::pow(x[1] - center[1], 2)));
  return gaussian(r, sigma);
}

int main(int argc, char *argv[]) {
  using cell_type    = cell::triangle;
  using fe_type      = cell_type::fe::lagrange_p1;
  using fes_type     = finite_element_space<fe_type>;
  using element_type = fes_type::element;
  using v_quad_type  = quad::triangle::qf5pT;
  
  const std::size_t n(100);
  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());

  fes_type fes(m, dm);

  bilinear_form<fes_type, fes_type> a(fes, fes); {
    const auto u(a.get_trial_function());
    const auto v(a.get_test_function());
    
    a += integrate<v_quad_type>(d<1>(u) * d<1>(v) +
				d<2>(u) * d<2>(v)
				, m);
  }

  linear_form<fes_type> f(fes); {
    const auto v(f.get_test_function());

    f += integrate<v_quad_type>(make_expr(src) * v
				, m);
  }

  element_type solution(a.solve(f));

  /*
   *  Compute the a posteriori error estimate
   */
  array<double> du{m.get_cell_number(), 2};
  double barycenter[] = {1.0 / 3.0, 1.0 / 3.0};
  for (std::unsigned int k(0); k < m.get_cell_number(); ++k) {
    const unsigned int d1[] = {1, 0};
    const unsigned int d2[] = {0, 1};
    du.at(k, 0) = solution.evaluate(k, barycenter, d1);
    du.at(k, 1) = solution.evaluate(k, barycenter, d1);
  }

  const std::size_t subdomain_id(cell_type::n_subdomain_type - 2);
  array<double> nu{m.get_cell_number()};
  nu.fill(0.0);
  for (std::unsigned int k(0); k < m.get_cell_number(); ++k) {
    for (std::unsigned int j(0); j < cell_type::n_subdomain(subdomain_id)) {
      const array<double>& vertices(m.get_vertices());
      const array<unsigned int>& cells(m.get_cells());
      
      
      double normal[] = {- (vertices.at(cells.at(k, 1), 1) - vertices.at(cells.at(k, 0), 1)),
                           (vertices.at(cells.at(k, 1), 0) - vertices.at(cells.at(k, 0), 0))};
      
      nu.at(k) += std::pow((du.at(k, 0) - du.at(, 0)), 2);
    }
  }
  
  exporter::ensight6("poisson"
		     , to_mesh_vertex_data<fe_type>(solution), "solution"
		     , to_mesh_vertex_data<fe_type>(projector::lagrange<fe_type>(src, fes)), "source");

  
  return 0;
}
