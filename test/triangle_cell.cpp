
#include "../src/mesh.hpp"
#include "../src/fe.hpp"
#include "../src/fes.hpp"
#include "../src/quadrature.hpp"
#include "../src/projector.hpp"
#include "../src/export.hpp"
#include "../src/timer.hpp"

double func(const double* x) {
  //return x[0] + x[1];
  //return 1.0;
  return std::sin(x[0] * M_PI) * std::sin(x[1] * M_PI);
}

void projection(std::size_t n) {
  timer t;
  
  const bool enable_stdout_export(false);
  const std::size_t n_x(n), n_y(n);
  
  try {
    mesh<cell::triangle> m(gen_square_mesh(1.0, 1.0, n_x, n_y));
    std::cout << std::setw(40) << "mesh: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
    
    typedef finite_element::triangle_lagrange_p1 fe_type;
    finite_element_space<fe_type> fes(m);
    std::cout << std::setw(40) << "finite element space: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;

    const auto func_h(projector::l2<fe_type, quad::triangle::qf5pT>(func, fes));
    std::cout << std::setw(40) << "projection: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
    
    const double h_x(1.0 / n_x), h_y(1.0 / n_y);
    if (enable_stdout_export) {
      std::cerr.precision(12);
      for (unsigned int i(0); i < n_x; ++i)
	for (unsigned int j(0); j < n_y; ++j) {
	  const double x[] = {i * h_x, j * h_y};
	  std::cerr << x[0] << " " << x[1] << " "
		    << func(x) << " "
		    << func_h.evaluate(x) << '\n';
	}
    }

    exporter::ensight6("l2projection", func_h, "func_h");
    std::cout << std::setw(40) << "ensight export: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
    
  } catch (const std::string& e) {
    std::cout << e << std::endl;
  }
}

int main(int , char**) {
  std::vector<std::size_t> ns{20, 20, 40, 80, 160, 320, 640, 1280};

  for (auto n: ns) {
    projection(n);
    std::cout << std::endl;
  }

  return 0;
}
