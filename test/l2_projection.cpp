#include "../src/core/projector.hpp"
#include "../src/core/export.hpp"

double func(const double* x) {
  return std::sin(x[0] * M_PI);
}

int main(int argc, char *argv[]) {
  typedef cell::edge cell_type;
  typedef finite_element::edge_lagrange_p1 fe_type;
  typedef finite_element_space<fe_type> fes_type;

  const bool export_to_stdout(true);
  const std::size_t n(4);
  
  try {
    mesh<cell_type> m(gen_segment_mesh(0.0, 1.0, n));
    fes_type fes(m);
    const auto func_h(projector::l2<fe_type, quad::edge::gauss3>(func, fes));

    if (export_to_stdout) {
      std::size_t n_k(n), n_i(20);
      
      const double h_k(1.0 / n_k);
      const double h_i(h_k / n_i);
      std::cerr << "# x, func, func_h \n";
      for (std::size_t k(0); k < n_k; ++k)
	for (std::size_t i(0); i <= n_i; ++i) {
	  const double x_hat(1.0 / n_i * i);
	  const double x(k * h_k + i * h_i);
	  std::cerr << k * h_k + i * h_i << " "
		    << func(&x) << " "
		    << func_h.evaluate(k, &x_hat) << std::endl;
	}
    }
  } catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
