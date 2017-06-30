
#include <cmath>
#include <iostream>

#include "../src/quadrature.hpp"

template<typename Q>
double error(unsigned int n) {
  // integrate sin(x) on [0,pi]
  return std::abs(integrate<Q, double(*)(double)>(std::sin, 0, M_PI, n)
		  - (- std::cos(M_PI) + std::cos(0.0))) / std::abs((- std::cos(M_PI) + std::cos(0.0)));
}

int main(int argc, char *argv[]) {
  const double length(M_PI);
  for (unsigned int i(1); i < 16; ++i) {
    const std::size_t n(std::pow(2, i + 1));
    const double h(length / static_cast<double>(n));
    std::cout.precision(16);
    std::cout << h << " "
	      << error<quad::gauss1>(n) << " "
      	      << error<quad::gauss2>(n) << " "
      	      << error<quad::gauss3>(n) << " "
      	      << error<quad::gauss4>(n) << " "
      	      << error<quad::gauss5>(n) << std::endl;

  }
  
  return 0;
}

