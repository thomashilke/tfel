#include "quadrature.hpp"

#include <cmath>
using namespace std;

const double quad::gauss1::x[1] = {0.0};
const double quad::gauss1::w[1] = {2.0};

const double quad::gauss2::x[2] = {-sqrt(1.0/3.0), sqrt(1.0/3.0)};
const double quad::gauss2::w[2] = {1.0, 1.0};


const double quad::gauss3::x[3] = {-sqrt(3.0/5.0), 0.0 , sqrt(3.0/5.0)};
const double quad::gauss3::w[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};


const double quad::gauss4::x[4] = {-sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)),
				    sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)),
				   -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0)),
				    sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0))};
const double quad::gauss4::w[4] = {(18.0 + sqrt(30.0))/36.0,
				   (18.0 + sqrt(30.0))/36.0,
				   (18.0 - sqrt(30.0))/36.0,
				   (18.0 - sqrt(30.0))/36.0};


const double quad::gauss5::x[5] = {-1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0)),
				    1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0)),
				    0.0,
				   -1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0)),
				    1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0))};
const double quad::gauss5::w[5] = { (322.0 + 13.0*sqrt(70.0)) / 900.0,
				    (322.0 + 13.0*sqrt(70.0)) / 900.0,
				    128.0/225.0,
				    (322.0 - 13.0*sqrt(70.0)) / 900.0,
				    (322.0 - 13.0*sqrt(70.0)) / 900.0};

