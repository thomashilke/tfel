
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <cmath>

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
  return std::exp(x[0]);
}


void test_expression() {
  const std::size_t k(0);
  const double x(0.5), x_hat(0.5);
  const double phi[2] = {std::sin(x), std::cos(x)};
  const double psi[2] = {x * x, 2.0 * x};

  test_function_t _1((form<0,0,0>()));
  trial_function_t _2((form<1,0,0>()));
  
  assert(_1(k, &x, &x_hat, phi, psi) == phi[0]);
  assert(_2(k, &x, &x_hat, phi, psi) == psi[0]);
  assert(d<1>(_1)(k, &x, &x_hat, phi, psi) == phi[1]);
  assert(d<1>(_2)(k, &x, &x_hat, phi, psi) == psi[1]);
  assert((_1 + _2)(k, &x, &x_hat, phi, psi) == phi[0] + psi[0]);
  assert((_1 - _2)(k, &x, &x_hat, phi, psi) == phi[0] - psi[0]);
  assert((_1 * _2)(k, &x, &x_hat, phi, psi) == phi[0] * psi[0]);

  assert(d<1>(_1 + _2)(k, &x, &x_hat, phi, psi) == phi[1] + psi[1]);
  assert(d<1>(_1 - _2)(k, &x, &x_hat, phi, psi) == phi[1] - psi[1]);
  assert(d<1>(_1 * _2)(k, &x, &x_hat, phi, psi) == phi[1] * psi[0] + phi[0] * psi[1]);

  assert((5.0 * _1)(k, &x, &x_hat, phi, psi) == 5.0 * phi[0]);
  assert((5.0 * d<1>(_1))(k, &x, &x_hat, phi, psi) == 5.0 * phi[1]);
  assert((5.0 * _2)(k, &x, &x_hat, phi, psi) == 5.0 * psi[0]);
  assert((5.0 * d<1>(_2))(k, &x, &x_hat, phi, psi) == 5.0 * psi[1]);

  assert((g * _1)(k, &x, &x_hat, phi, psi) == g(&x) * phi[0]);
  assert((g * _1 * _2)(k, &x, &x_hat, phi, psi) == g(&x) * phi[0] * psi[0]);

  assert((_1 * _1 * _1)(k, &x, &x_hat, phi, psi) == phi[0] * phi[0] * phi[0]);
}





int main(int argc, char *argv[]) {
  trial_function_t psi((form<1,0,0>()));
  test_function_t phi((form<0,0,0>()));

  auto mass(phi * psi);
  auto mass2(phi * psi + phi * psi);
  auto mass21(10. * phi * psi + 0.5 * phi * psi - g * phi * psi);

  auto dmass(d<1>(phi * psi));
  
  print_type<form<0, 0, 0>>();
  print_type<form<1, 0, 0>>();
  print_type(mass);
  print_type(dmass);

  argument<2>(1.0, 'c', 67ul);
  
  return 0;
}
