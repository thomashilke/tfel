
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <cmath>

#include "../src/expression.hpp"
#include "../src/meta.hpp"


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

void test_expression_call_wrapper() {
  std::size_t k(0);
  const double x(0.5), x_hat(0.5);
  const double phi[2] = {std::sin(x), std::cos(x)};
  const double psi[2] = {x * x, 2.0 * x};

  const double* psi_phi[2] = {phi, psi};
  
  test_function_t _1((form<0,0,0>()));
  trial_function_t _2((form<1,0,0>()));

  std::cout << "wrapped call: " << expression_call_wrapper<0, 2>::call(_1, psi_phi, k, &x, &x_hat) << std::endl;
  std::cout << "std call: " << _1(k, &x, &x_hat, phi, psi) << std::endl;
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
  print_type(mass2);
  print_type(mass21);
  print_type(dmass);
  
  argument<2>(1.0, 'c', 67ul);

  test_expression_call_wrapper();
  
  return 0;
}
