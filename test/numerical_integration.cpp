
#include <iostream>

#include "../src/core/quadrature.hpp"

double f_1(double x) {
  return x;
}


int main(int argc, char *argv[]) {
  std::cout << integrate<quad::gauss1>(f_1, 0.0, 1.0, 10) << std::endl;
  std::cout << integrate<quad::gauss2>(f_1, 0.0, 1.0, 10) << std::endl;
  std::cout << integrate<quad::gauss3>(f_1, 0.0, 1.0, 10) << std::endl;
  std::cout << integrate<quad::gauss4>(f_1, 0.0, 1.0, 10) << std::endl;
  std::cout << integrate<quad::gauss5>(f_1, 0.0, 1.0, 10) << std::endl;

  return 0;
}
