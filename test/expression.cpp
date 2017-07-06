
#include <iostream>
#include <typeinfo>

#include "../src/expression.hpp"

template<typename T>
void print_type() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

template<typename T>
void print_type(const T& e) {
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

double g(const double* x) {
  return 1.0;
}
  
int main(int argc, char *argv[]) {
  test_function_t phi((test_function<0,0>()));
  trial_function_t psi((trial_function<0,0>()));
  constant_t c(constant(10));

  auto mass(phi * psi);
  auto mass2(phi * psi + phi * psi);
  auto mass21(10. * phi * psi + 0.5 * phi * psi - g * phi * psi);

  auto dmass(d<1>(phi * psi));
  
  print_type<test_function<0, 0>>();
  print_type<trial_function<0, 0>>();
  print_type(mass);
  print_type(dmass);
  
  return 0;
}
