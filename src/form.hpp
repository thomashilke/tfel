#ifndef FORM_H
#define FORM_H

#include <cstddef>
#include <ostream>
#include <algorithm>

#include <spikes/array.hpp>

#include "mesh.hpp"
#include "expression.hpp"
#include "sparse_linear_system.hpp"
#include "linear_solver.hpp"


template<typename cell_t, typename quadrature_t, typename form_t>
struct mesh_integration_proxy {
  mesh_integration_proxy(const form_t& f, const mesh<cell_t>& m)
    : f(f), m(m), xq{quadrature_t::n_point, quadrature_type::n_point_space_dimension} {
    xq.set_data(&quadrature_t::x[0][0]);
  }

  typedef quadrature_t quadrature_type;
  typedef form_t form_type;
  typedef cell_t cell_type;
  
  const form_type f;
  const mesh<cell_type>& m;
  array<double> xq;

  std::size_t get_global_element_id(std::size_t k) const {
    return k;
  }

  array<double> get_quadrature_points(std::size_t /*k*/) const {
    return xq;
  }
};

template<typename cell_t, typename quadrature_t, typename form_t>
struct submesh_integration_proxy {
  submesh_integration_proxy(const form_t& f, const submesh<cell_t>& m)
    : f(f), m(m) {}

  typedef quadrature_t quadrature_type;
  typedef form_t form_type;
  typedef typename submesh<cell_t>::cell_type cell_type;
  
  const form_type f;
  const submesh<cell_t>& m;
  
  std::size_t get_global_element_id(std::size_t k) const {
    return m.get_parent_element_id(k);
  }
  
  array<double> get_quadrature_points(std::size_t k) const {
    const std::size_t n_q(quadrature_type::n_point);
    array<double> xq{n_q, quadrature_type::n_point_space_dimension};
    xq.set_data(&quadrature_type::x[0][0]);
    array<double> hat_xq(cell_t::map_points_to_subdomain(m.get_subdomain_id(k), xq));
    return hat_xq;
  }
};


template<typename quadrature_type, typename cell_type, typename form_type>
typename std::enable_if<(form_type::rank > 0),
		 mesh_integration_proxy<cell_type, quadrature_type, expression<form_type> >
		 >::type
integrate(const expression<form_type>& f, const mesh<cell_type>& m) {
  return mesh_integration_proxy<cell_type, quadrature_type, expression<form_type> >(f, m);
}

template<typename quadrature_type, typename cell_type, typename form_type>
typename std::enable_if<(form_type::rank > 0),
		 submesh_integration_proxy<cell_type, quadrature_type, expression<form_type> >
		 >::type

integrate(const expression<form_type>& f, const submesh<cell_type>& m) {
  return submesh_integration_proxy<cell_type, quadrature_type, expression<form_type> >(f, m);
}

template<typename T>
double integrate_with_proxy(const T& integration_proxy) {
  static_assert(T::form_type::rank == 0, "integrate_with_proxy expects rank-0 expression.");

  typedef typename T::quadrature_type quadrature_type;
  typedef typename T::cell_type cell_type;

  const auto& m(integration_proxy.m);

  // prepare the quadrature weights
  const std::size_t n_q(quadrature_type::n_point);
  array<double> omega{n_q};
  omega.set_data(&quadrature_type::w[0]);

  double result(0.0);
  for (unsigned int k(0); k < m.get_element_number(); ++k) {
    const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
    const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
								      m.get_elements(),
								      k, xq_hat));
    // evaluate the expression
    const double volume(m.get_cell_volume(k));
    double rhs_el(0.0);
    for (unsigned int q(0); q < n_q; ++q) {
      rhs_el += volume * omega.at(q)
	* integration_proxy.f(k, &xq.at(q, 0), &xq_hat.at(q, 0));
    }
    result += rhs_el;
  }
  return result;
}

template<typename quadrature_type, typename cell_type, typename form_type>
typename std::enable_if<form_type::rank == 0, double>::type
integrate(const expression<form_type>& f, const mesh<cell_type>& m) {
  mesh_integration_proxy<cell_type, quadrature_type, expression<form_type> > proxy(f, m);
  return integrate_with_proxy(proxy);
}

template<typename quadrature_type, typename cell_type, typename form_type>
typename std::enable_if<form_type::rank == 0, double>::type
integrate(const expression<form_type>& f, const submesh<cell_type>& m) {
  submesh_integration_proxy<cell_type, quadrature_type, expression<form_type> > proxy(f, m);
  return integrate_with_proxy(proxy);
}


#include "linear_form.hpp"
#include "bilinear_form.hpp"

#endif /* FORM_H */
