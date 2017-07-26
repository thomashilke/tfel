
#include <cmath>
#include <iostream>

#include "../src/core/mesh.hpp"
#include "../src/core/quadrature.hpp"
#include "../src/core/cell.hpp"
#include "../src/core/expression.hpp"
#include "../src/core/form.hpp"

double f_2(const double* x) {
  return std::sin(x[0]) + std::sin(x[1]);
}

template<typename Q>
double error(unsigned int n) {
  const mesh<cell::triangle> m(gen_square_mesh(M_PI, M_PI, n, n));
  return std::abs(integrate<Q>(make_expr(f_2), m) - 4.0 * M_PI); 
}


int main(int argc, char *argv[]) {
  for (unsigned int i(1); i < 12; ++i) {
    const std::size_t n(std::pow(2, i + 1));
    const double h(M_PI / static_cast<double>(n));
    std::cout.precision(16);

    std::cout << h << " "
	      << error<quad::triangle::qf1pT>(n) << " "
      	      << error<quad::triangle::qf2pT>(n) << " "
      	      << error<quad::triangle::qf5pT>(n) << " "
      	      << error<quad::triangle::qf1pTlump>(n) << std::endl;
  }
  
  return 0;
}

