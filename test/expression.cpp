
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <cmath>

#include "../src/expression.hpp"
#include "../src/meta.hpp"
#include "../src/fe_value_manager.hpp"
#include "../src/fe.hpp"
#include "../src/composite_form.hpp"

typedef expression<form<1, 0, 0> > trial_function_t;
typedef expression<form<0, 0, 0> > test_function_t;

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


void test_fe_value_manager() {
  std::size_t nq(1);
  
  using fe = finite_element::triangle_lagrange_p1;
  using fe_list = type_list<fe, fe, fe>;
  
  fe_value_manager<fe_list> values(nq), zvalues(nq);
  values.clear();
  zvalues.clear();

  const double
    *phi_1(&values.get_values<0>().at(0,0,0)),
    *phi_2(&values.get_values<1>().at(0,0,0)),
    *phi_3(&values.get_values<2>().at(0,0,0));

  const double
    *zphi_1(&zvalues.get_values<0>().at(0,0,0)),
    *zphi_2(&zvalues.get_values<1>().at(0,0,0)),
    *zphi_3(&zvalues.get_values<2>().at(0,0,0));

  assert(phi_1 == phi_2);
  assert(phi_2 == phi_3);

  assert(zphi_1 == zphi_2);
  assert(zphi_2 == zphi_3);

  array<double> jmt{2, 2};
  jmt.at(0, 0) = 1.0;
  jmt.at(1, 0) = 0.0;
  jmt.at(0, 1) = 0.0;
  jmt.at(1, 1) = 1.0;

  array<double> xq{1, 2};
  xq.at(0, 0) = 1.0/3.0;
  xq.at(0, 1) = 1.0/3.0;

  values.prepare(jmt, xq);

  assert(phi_1[0] == fe::phi(0, &xq.at(0,0)));
  assert(phi_1[1] == fe::phi(1, &xq.at(0,0)));
  assert(phi_1[2] == fe::phi(2, &xq.at(0,0)));

  assert(phi_1[3] == fe::dphi(0, 0, &xq.at(0,0)));
  assert(phi_1[4] == fe::dphi(0, 1, &xq.at(0,0)));
  assert(phi_1[5] == fe::dphi(0, 2, &xq.at(0,0)));

  assert(phi_1[6] == fe::dphi(1, 0, &xq.at(0,0)));
  assert(phi_1[7] == fe::dphi(1, 1, &xq.at(0,0)));
  assert(phi_1[8] == fe::dphi(1, 2, &xq.at(0,0)));
}

void test_valuation_selection() {
  std::size_t nq(1);
  
  using fe_0 = finite_element::triangle_lagrange_p0;
  using fe_1 = finite_element::triangle_lagrange_p1;
  using fe_2 = finite_element::triangle_lagrange_p1_bubble;
  using fe_list = type_list<fe_0, fe_1, fe_2>;
  
  fe_value_manager<fe_list> values(nq), zvalues(nq);
  values.clear();
  zvalues.clear();

  const double
    *phi_1(&values.get_values<0>().at(0,0,0)),
    *phi_2(&values.get_values<1>().at(0,0,0)),
    *phi_3(&values.get_values<2>().at(0,0,0));

  const double
    *zphi_1(&zvalues.get_values<0>().at(0,0,0)),
    *zphi_2(&zvalues.get_values<1>().at(0,0,0)),
    *zphi_3(&zvalues.get_values<2>().at(0,0,0));

  assert(phi_1 != phi_2);
  assert(phi_2 != phi_3);

  assert(zphi_1 != zphi_2);
  assert(zphi_2 != zphi_3);

  array<double> jmt{2, 2};
  jmt.at(0, 0) = 1.0;
  jmt.at(1, 0) = 0.0;
  jmt.at(0, 1) = 0.0;
  jmt.at(1, 1) = 1.0;

  array<double> xq{1, 2};
  xq.at(0, 0) = 1.0/3.0;
  xq.at(0, 1) = 1.0/3.0;

  values.prepare(jmt, xq);

  const double* phi_psi[6] = {};
  select_function_valuation<fe_list, 0>(phi_psi, 0, 0, values, zvalues);
  assert(phi_psi[0] == phi_1);
  assert(phi_psi[1] == zphi_2);
  assert(phi_psi[2] == zphi_3);
  
  select_function_valuation<fe_list, 1>(phi_psi, 0, 0, values, zvalues);
  assert(phi_psi[0] == zphi_1);
  assert(phi_psi[1] == phi_2);
  assert(phi_psi[2] == zphi_3);
  
  select_function_valuation<fe_list, 2>(phi_psi, 0, 0, values, zvalues);
  assert(phi_psi[0] == zphi_1);
  assert(phi_psi[1] == zphi_2);
  assert(phi_psi[2] == phi_3);


  for (int i(0); i < 3; ++i) {
    select_function_valuation<fe_list, 1>(phi_psi, 0, i, values, zvalues);
    assert(&zvalues.get_values<0>().at(0,0,0) == phi_psi[0]);
    assert(&values.get_values<1>().at(0,i,0) == phi_psi[1]);
    assert(&zvalues.get_values<2>().at(0,0,0) == phi_psi[2]);
  }

  for (int i(0); i < 4; ++i) {
    select_function_valuation<fe_list, 2>(phi_psi, 0, i, values, zvalues);
    assert(&zvalues.get_values<0>().at(0,0,0) == phi_psi[0]);
    assert(&zvalues.get_values<1>().at(0,0,0) == phi_psi[1]);
    assert(&values.get_values<2>().at(0,i,0) == phi_psi[2]);
  }
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

void test_basis_function() {
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
}

int main(int argc, char *argv[]) {
  //test_basis_function();
  test_expression();
  test_valuation_selection();
  //  test_expression_call_wrapper();
  
  return 0;
}
