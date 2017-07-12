#ifndef FORM_H
#define FORM_H

#include <cstddef>
#include <ostream>

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
  
  const mesh<cell_type>& m;
  const form_type f;
  array<double> xq;

  std::size_t get_global_element_id(std::size_t k) const {
    return k;
  }

  array<double> get_quadrature_points(std::size_t k) const {
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
  
  const submesh<cell_t>& m;
  const form_type f;
  
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
  const std::size_t dim(m.get_embedding_space_dimension());

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

template<typename test_fes_type>
class linear_form {
public:
  linear_form(const test_fes_type& te_fes)
    : test_fes(te_fes), f{te_fes.get_dof_number()} {}

  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 1, "linear_form expects rank-1 expression.");
    
    typedef typename test_fes_type::fe_type test_fe_type;
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    const auto& m(integration_proxy.m);
    const std::size_t dim(m.get_embedding_space_dimension());

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);


    // storage for the point-wise basis function evaluation
    const std::size_t n_dof(test_fe_type::n_dof_per_element);
    array<double> psi{n_q, n_dof, dim + 1};

    
    // loop over the elements
    for (unsigned int k(0); k < m.get_element_number(); ++k) {
      const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
      const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
									m.get_elements(),
									k, xq_hat));

      // prepare the basis function on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i)
	  psi.at(q, i, 0) = test_fe_type::phi(i, &xq_hat.at(q, 0));
      }
      
      // prepare the basis function derivatives on the quadrature points
      const array<double> jmt(m.get_jmt(k));
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i) {
	  for (std::size_t n(0); n < dim; ++n) {
	    psi.at(q, i, 1 + n) = 0.0;
	    for (std::size_t k(0); k < dim; ++k)
	      psi.at(q, i, 1 + n) += jmt.at(n, k) * test_fe_type::dphi(k, i, &xq_hat.at(q, 0));
	  }
	}
      }

      // evaluate the weak form
      const double volume(m.get_cell_volume(k));
      for (unsigned int j(0); j < n_dof; ++j) {
	double rhs_el(0.0);
	for (unsigned int q(0); q < n_q; ++q) {
	  rhs_el += volume * omega.at(q) * integration_proxy.f(k, &xq.at(q, 0), &xq_hat.at(q, 0),
							       &psi.at(q, j, 0));
	}
	f.at(test_fes.get_dof(integration_proxy.get_global_element_id(k), j)) += rhs_el;
      }
    }
  }

  expression<form<0,0,0> > get_test_function() const { return form<0,0,0>(); }

  const array<double>& get_components() const { return f; }
  
  void show(std::ostream& stream) {
    stream << "rhs = [" << f.at(0);
    for (std::size_t j(1); j < f.get_size(0); ++j)
      stream << "; " << f.at(j);
    stream << "];" << std::endl;
  }
private:
  const test_fes_type& test_fes;
  
  array<double> f;
};

template<typename test_fes_type, typename trial_fes_type>
class bilinear_form {
public:
  bilinear_form(const test_fes_type& te_fes,
		const trial_fes_type& tr_fes)
    : test_fes(te_fes), trial_fes(tr_fes),
      a(te_fes.get_dof_number(), tr_fes.get_dof_number()) {
    for (const auto& i: test_fes.get_dirichlet_dof())
      a.accumulate(i, i, 1.0);
  }

  template<typename T>
  void operator+=(const T& integration_proxy) {
    static_assert(T::form_type::rank == 2, "linear_form expects rank-2 expression.");
      
    typedef typename test_fes_type::fe_type test_fe_type;
    typedef typename trial_fes_type::fe_type trial_fe_type;
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;

    const auto& m(integration_proxy.m);
    const std::size_t dim(m.get_embedding_space_dimension());

    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);


    // storage for the point-wise basis function evaluation
    const std::size_t n_dof(test_fe_type::n_dof_per_element);
    array<double> phi{n_q, n_dof, dim + 1};


    // loop over the elements
    for (unsigned int k(0); k < m.get_element_number(); ++k) {
      const array<double> xq_hat(integration_proxy.get_quadrature_points(k));
      const array<double> xq(cell_type::map_points_to_space_coordinates(m.get_vertices(),
									m.get_elements(),
									k, xq_hat));
      // prepare the basis function on the quadrature points
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i)
	  phi.at(q, i, 0) = test_fe_type::phi(i, &xq_hat.at(q, 0));
      }
      
      // prepare the basis function derivatives on the quadrature points
      const array<double> jmt(m.get_jmt(k));
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i) {
	  for (std::size_t n(0); n < dim; ++n) {
	    phi.at(q, i, 1 + n) = 0.0;
	    for (std::size_t k(0); k < dim; ++k)
	      phi.at(q, i, 1 + n) += jmt.at(n, k) * test_fe_type::dphi(k, i, &xq_hat.at(q, 0));
	  }
	}
      }

      // evaluate the weak form
      const double volume(m.get_cell_volume(k));
      for (unsigned int i(0); i < n_dof; ++i) {
	for (unsigned int j(0); j < n_dof; ++j) {
	  double a_el(0.0);
	  for (unsigned int q(0); q < n_q; ++q) {
	    a_el += volume * omega.at(q) * integration_proxy.f(k,
							       &xq.at(q, 0), &xq_hat.at(q, 0),
							       &phi.at(q, i, 0),
							       &phi.at(q, j, 0));
	  }
	  accumulate(test_fes.get_dof(integration_proxy.get_global_element_id(k), j),
		     trial_fes.get_dof(integration_proxy.get_global_element_id(k), i), a_el);
	}
      }
    }
  }

  expression<form<0,0,0> > get_test_function() const { return form<0,0,0>(); }
  expression<form<1,0,0> > get_trial_function() const { return form<1,0,0>(); }

  typename trial_fes_type::element solve(const linear_form<test_fes_type>& form) const {
    array<double> f(form.get_components());
    for (const auto& i: test_fes.get_dirichlet_dof())
      f.at(i) = 0.0;
    
    linear_solver s;
    auto petsc_gmres_ilu(s.get_solver(solver::petsc,
				      method::gmres,
				      preconditioner::ilu,
				      test_fes.get_dof_number()));
    for (const auto& v: a.get_elements())
      petsc_gmres_ilu->add_value(v.first.first, v.first.second, v.second);
    petsc_gmres_ilu->assemble();

    typename trial_fes_type::element result(trial_fes, petsc_gmres_ilu->solve(f));
    delete petsc_gmres_ilu; petsc_gmres_ilu = nullptr;
    return result;
  }
  
  void show(std::ostream& stream) {
    a.show(stream);
  }
  
private:
  const test_fes_type& test_fes;
  const trial_fes_type& trial_fes;
  
  sparse_linear_system a;

  void accumulate(std::size_t i, std::size_t j, double value) {
    if (test_fes.get_dirichlet_dof().count(i) == 0)
      a.accumulate(i, j, value);
  }
};

#endif /* FORM_H */
