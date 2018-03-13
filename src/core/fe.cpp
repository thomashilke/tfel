#include "fe.hpp"

constexpr std::size_t cell::edge::fe::lagrange_p0::n_dof[2];
constexpr std::size_t cell::edge::fe::lagrange_p1::n_dof[2];
constexpr std::size_t cell::edge::fe::lagrange_p1_bubble::n_dof[2];
constexpr std::size_t cell::edge::fe::lagrange_p2::n_dof[2];

constexpr double cell::edge::fe::lagrange_p0::x[1][1];
constexpr double cell::edge::fe::lagrange_p1::x[2][1];
constexpr double cell::edge::fe::lagrange_p1_bubble::x[3][1];
constexpr double cell::edge::fe::lagrange_p2::x[3][1];


constexpr std::size_t cell::triangle::fe::lagrange_p0::n_dof[3];
constexpr std::size_t cell::triangle::fe::lagrange_p1::n_dof[3];
constexpr std::size_t cell::triangle::fe::lagrange_p1_bubble::n_dof[3];
constexpr std::size_t cell::triangle::fe::lagrange_p2::n_dof[3];

constexpr double cell::triangle::fe::lagrange_p0::x[1][2];
constexpr double cell::triangle::fe::lagrange_p1::x[3][2];
constexpr double cell::triangle::fe::lagrange_p1_bubble::x[4][2];
constexpr double cell::triangle::fe::lagrange_p2::x[6][2];


constexpr std::size_t cell::tetrahedron::fe::lagrange_p0::n_dof[4];
constexpr std::size_t cell::tetrahedron::fe::lagrange_p1::n_dof[4];

constexpr double cell::tetrahedron::fe::lagrange_p0::x[1][3];
constexpr double cell::tetrahedron::fe::lagrange_p1::x[4][3];
