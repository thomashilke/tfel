#include "fe.hpp"

constexpr std::size_t finite_element::edge_lagrange_p0::n_dof[2];
constexpr std::size_t finite_element::edge_lagrange_p1::n_dof[2];
constexpr std::size_t finite_element::edge_lagrange_p1_bubble::n_dof[2];

constexpr std::size_t finite_element::triangle_lagrange_p0::n_dof[3];
constexpr std::size_t finite_element::triangle_lagrange_p1::n_dof[3];
constexpr std::size_t finite_element::triangle_lagrange_p1_bubble::n_dof[3];

constexpr std::size_t finite_element::tetrahedron_lagrange_p0::n_dof[4];
constexpr std::size_t finite_element::tetrahedron_lagrange_p1::n_dof[4];


constexpr double finite_element::edge_lagrange_p0::x[1][1];
constexpr double finite_element::edge_lagrange_p1::x[2][1];
constexpr double finite_element::edge_lagrange_p1_bubble::x[3][1];

constexpr double finite_element::triangle_lagrange_p0::x[1][2];
constexpr double finite_element::triangle_lagrange_p1::x[3][2];
constexpr double finite_element::triangle_lagrange_p1_bubble::x[4][2];

constexpr double finite_element::tetrahedron_lagrange_p0::x[1][3];
constexpr double finite_element::tetrahedron_lagrange_p1::x[4][3];
