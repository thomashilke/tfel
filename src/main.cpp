
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"
#include "mesh.hpp"
#include "fes.hpp"
#include "quadrature.hpp"
#include "expression.hpp"
#include "timer.hpp"
#include "sparse_linear_system.hpp"
#include "form.hpp"
#include "linear_solver.hpp"
#include "export.hpp"

double bell(double x) {
  const double epsilon(1e-6);
  if (std::abs(x) < 1.0 - epsilon) {
    return std::exp(1.0 - 1.0 / (1.0 - x * x));
  } else {
    return 0.0;
  }
}

double shift_scale_bell(double x, double width, double x_0) {
  return bell((x - x_0) / width);
}

double force(const double* x) {
  return shift_scale_bell(x[0], 0.25, 0.5);
}

double g(const double* x) {
  return 0.0;
}

int main(int argc, char *argv[]) {
  timer t;
  try {
    const std::size_t n(10000);
    
    mesh<cell::edge> m(gen_segment_mesh(0.0, 1.0, n));
    submesh<cell::edge> dm(m.get_boundary_submesh());
    submesh<cell::edge> left_boundary(dm.query_elements([](const double* x) -> bool {
	  return *x < 1.0 / 2.0;
	}));
    submesh<cell::edge> right_boundary(dm.query_elements([](const double* x) -> bool {
	  return *x > 1.0 / 2.0;
	}));
    std::cout << std::setw(40) << "mesh: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
 
    typedef finite_element::edge_lagrange_p1 fe;
    finite_element_space<fe> fes(m, left_boundary);
    std::cout << std::setw(40) << "finite element system: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
    
    bilinear_form<finite_element_space<fe>, finite_element_space<fe> >
      a(fes, fes); {
      const auto& phi(a.get_trial_function());
      const auto& psi(a.get_test_function());
    
      a += integrate<quad::edge::gauss3>(d<1>(phi) * d<1>(psi), m);
      std::cout << std::setw(40) << "bilinear form volume integration: "
		<< std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
    }

    linear_form<finite_element_space<fe> > f(fes); {
      const auto& psi(f.get_test_function());

      f += integrate<quad::edge::gauss3>(force * psi, m);
      std::cout << std::setw(40) << "linear form volume integration: "
		<< std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
      
      f += integrate<quad::point::eval>(g * psi, right_boundary);
      std::cout << std::setw(40) << "linear form boundary integration: "
		<< std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
    }

    typedef finite_element_space<fe>::element element_type;
    const element_type u(a.solve(f));
    std::cout << std::setw(40) << "solve linear system: "
		<< std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;


    exporter::ascii("u.dat", u);
    std::cout << std::setw(40) << "export the solution: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
    
  return 0;
}

