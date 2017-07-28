#include "quadrature.hpp"

#include <cmath>
using namespace std;

const double quad::point::eval::x[1][1] = {{0.0}};
const double quad::point::eval::w[1] = {1.0};


const double quad::edge::gauss1::x[1][1] = {{1.0/2.0}};
const double quad::edge::gauss1::w[1] = {1.0};


const double quad::edge::gauss2::x[2][1] = {{1.0/2.0 - sqrt(1.0/3.0) / 2.0},
					    {1.0/2.0 + sqrt(1.0/3.0) / 2.0}};
const double quad::edge::gauss2::w[2] = {1.0/2.0, 1.0/2.0};


const double quad::edge::gauss3::x[3][1] = {{1.0/2.0 - sqrt(3.0/5.0) / 2.0},
					    {1.0/2.0 + 0.0 / 2.0},
					    {1.0/2.0 + sqrt(3.0/5.0) / 2.0}};
const double quad::edge::gauss3::w[3] = {5.0/9.0 / 2.0, 8.0/9.0 / 2.0, 5.0/9.0 / 2.0};


const double quad::edge::gauss4::x[4][1] = { {1.0/2.0 - sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)) / 2.0},
					     {1.0/2.0 + sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)) / 2.0},
					     {1.0/2.0 - sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0)) / 2.0},
					     {1.0/2.0 + sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0)) / 2.0}};
const double quad::edge::gauss4::w[4] = {(18.0 + sqrt(30.0))/36.0 / 2.0,
					 (18.0 + sqrt(30.0))/36.0 / 2.0,
					 (18.0 - sqrt(30.0))/36.0 / 2.0,
					 (18.0 - sqrt(30.0))/36.0 / 2.0};


const double quad::edge::gauss5::x[5][1] = { {1.0/2.0 - 1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0)) / 2.0},
					     {1.0/2.0 + 1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0)) / 2.0},
					     {1.0/2.0 + 0.0 / 2.0},
					     {1.0/2.0 - 1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0)) / 2.0},
					     {1.0/2.0 + 1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0)) / 2.0}};
const double quad::edge::gauss5::w[5] = { (322.0 + 13.0*sqrt(70.0)) / 900.0 / 2.0,
					  (322.0 + 13.0*sqrt(70.0)) / 900.0 / 2.0,
					  128.0/225.0 / 2.0,
					  (322.0 - 13.0*sqrt(70.0)) / 900.0 / 2.0,
					  (322.0 - 13.0*sqrt(70.0)) / 900.0 / 2.0};


const double quad::triangle::qf1pT::x[1][2] = {{1.0/3.0, 1.0/3.0}};
const double quad::triangle::qf1pT::w[1] = {1.0};


const double quad::triangle::qf2pT::x[3][2] = {{1.0/2.0, 1.0/2.0},
					       {1.0/2.0, 0.0},
					       {0.0,     1.0/2.0}};
const double quad::triangle::qf2pT::w[3] = {1.0/3.0, 1.0/3.0, 1.0/3.0};


const double quad::triangle::qf5pT::x[7][2] = {{1.0/3.0, 1.0/3.0},
					       {(6.0       - sqrt(15.0)) / 21.0, (6.0       - sqrt(15.0)) / 21.0},
					       {(6.0       - sqrt(15.0)) / 21.0, (9.0 + 2.0 * sqrt(15.0)) / 21.0},
					       {(9.0 + 2.0 * sqrt(15.0)) / 21.0, (6.0       - sqrt(15.0)) / 21.0},
					       {(6.0       + sqrt(15.0)) / 21.0, (6.0       + sqrt(15.0)) / 21.0},
					       {(6.0       + sqrt(15.0)) / 21.0, (9.0 - 2.0 * sqrt(15.0)) / 21.0},
					       {(9.0 - 2.0 * sqrt(15.0)) / 21.0, (6.0       + sqrt(15.0)) / 21.0}};
const double quad::triangle::qf5pT::w[7] = {0.225,
					    (155.0 - sqrt(15.0)) / 1200.0,
					    (155.0 - sqrt(15.0)) / 1200.0,
					    (155.0 - sqrt(15.0)) / 1200.0,
					    (155.0 + sqrt(15.0)) / 1200.0,
					    (155.0 + sqrt(15.0)) / 1200.0,
					    (155.0 + sqrt(15.0)) / 1200.0};


const double quad::triangle::qf1pTlump::x[3][2] = {{0.0, 0.0},
						   {1.0, 0.0},
						   {0.0, 1.0}};
const double quad::triangle::qf1pTlump::w[3] = {1.0/3.0, 1.0/3.0, 1.0/3.0};


const double quad::tetrahedron::qf1pTet::x[1][3] = {{0.25, 0.25, 0.25}};
const double quad::tetrahedron::qf1pTet::w[1] = {1.0 / 2.0};

  

const double quad::tetrahedron::qf4pTet::x[4][3] = {{0.138196601125011, 0.138196601125011, 0.138196601125011},
						    {0.585410196624969, 0.138196601125011, 0.138196601125011},
						    {0.138196601125011, 0.585410196624969, 0.138196601125011},
						    {0.138196601125011, 0.138196601125011, 0.585410196624969}};
const double quad::tetrahedron::qf4pTet::w[4] = {0.25 / 2.0, 0.25 / 2.0, 0.25 / 2.0, 0.25 / 2.0};

/*const double quad::tetrahedron::qf5pTet::x[5][3] = {{},
						    {},
						    {},
						    {},
						    {}};
						    const double quad::tetrahedron::qf5pTet::w[5] = {/ 2.0, / 2.0, / 2.0, / 2.0, / 2.0};*/
