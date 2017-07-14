#include "fe.hpp"

constexpr std::size_t finite_element::edge_lagrange_p1::n_dof[2];
constexpr std::size_t finite_element::edge_lagrange_p1_bubble::n_dof[2];

constexpr double finite_element::edge_lagrange_p0::x[1][1];
constexpr double finite_element::edge_lagrange_p1::x[2][1];
constexpr double finite_element::edge_lagrange_p1_bubble::x[3][1];
constexpr double finite_element::triangle_lagrange_p0::x[1][2];
constexpr double finite_element::triangle_lagrange_p1::x[3][2];