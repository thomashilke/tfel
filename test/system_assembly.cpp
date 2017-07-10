
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>

#include <spikes/array.hpp>

#include "../src/cell.hpp"
#include "../src/fe.hpp"
#include "../src/mesh.hpp"
#include "../src/fes.hpp"
#include "../src/quadrature.hpp"
#include "../src/timer.hpp"
#include "../src/sparse_linear_system.hpp"


int main(int argc, char *argv[]) {
  timer t;
  try {
    const std::size_t n(6);
    
    mesh<cell::edge> m(gen_segment_mesh(0.0, 1.0, n));
    submesh<cell::edge> dm(m.get_boundary_submesh());
    submesh<cell::edge> left_boundary(dm.query_elements([](const double* x) -> bool {
	  return *x < 1.0 / 2.0;
	}));
    submesh<cell::edge> right_boundary(dm.query_elements([](const double* x) -> bool {
	  return *x > 1.0 / 2.0;
	}));
    std::cout << std::setw(40) << "mesh: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
  
 
    typedef finite_element::edge_lagrange_p1 fe;
    finite_element_space<fe> fes(m, left_boundary);
    fes.show(std::cout);
    
    sparse_linear_system sys(fes.get_dof_number(),
			     fes.get_dof_number());

    std::cout << std::setw(40) << "finite element system: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;

    // integration over m:
    {
      const std::size_t dim(m.get_embedding_space_dimension());

      // prepare the quadrature points and weights
      typedef quad::edge::gauss3 m_quad;
      const std::size_t n_q(m_quad::n_point);
      array<double> xq{n_q, m_quad::n_point_space_dimension};
      xq.set_data(&m_quad::x[0][0]);
      array<double> omega{n_q};
      omega.set_data(&m_quad::w[0]);

      // prepare the basis function on the quadrature points
      const std::size_t n_dof(fe::n_dof_per_element);
      array<double> phi{n_dof, n_q}; // That can be cached before the loops
      for (unsigned int q(0); q < n_q; ++q) {
	for (std::size_t i(0); i < n_dof; ++i)
	  phi.at(i, q) = fe::phi(i, &xq.at(q, 0));
      }

      for (unsigned int k(0); k < m.get_element_number(); ++k) {

	// prepare the quadrature point space coordinates
	array<double> hat_xq(fe::cell_type::map_points_to_space_coordinates(m.get_vertices(),
									    m.get_elements(),
									    k, xq));

	// prepare the basis function derivatives on the quadrature points
	const array<double> jmt(m.get_jmt(k));
	array<double> dphi{dim, n_dof, n_q}; // And that too
	for (unsigned int q(0); q < n_q; ++q) {
	  for (std::size_t i(0); i < n_dof; ++i) {
	    for (std::size_t n(0); n < dim; ++n) {
	      dphi.at(n, i, q) = 0.0;
	      for (std::size_t k(0); k < dim; ++k)
		dphi.at(n, i, q) += jmt.at(n, k) * fe::dphi(k, i, &xq.at(q, 0));
	    }
	  }
	}

	// evaluate the weak form
	const double volume(m.get_cell_volume(k));
	for (unsigned int i(0); i < n_dof; ++i) {
	  for (unsigned int j(0); j < n_dof; ++j) {
	    double a_el(0.0);
	    for (unsigned int q(0); q < n_q; ++q) {
	      // here we have to make available k, x, x_hat, phi, dphi, psi, dpsi;
	      // k, x_hat, phi, dphi are already availiable. psi and dpsi is trivial,
	      // and is x.
	      for (unsigned int n(0); n < dim; ++n)
		a_el += volume * omega.at(q) * dphi.at(n, i, q) * dphi.at(n, j, q);
	    }
	    sys.accumulate(fes.get_dof(k, j),
			   fes.get_dof(k, i), a_el);
	  }
	}
      }
    }
    std::cout << std::setw(40) << "bilinear form volume integration: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;
  
    if (false) {
      // integration over dm:
      // well, there's nothing in this case
      for (std::size_t k(0); k < dm.get_element_number(); ++k) {
	for (std::size_t i(0); i < fe::n_dof_per_element; ++i) {
	  for (std::size_t j(0); j < fe::n_dof_per_element; ++j) {
	    typedef quad::point::eval dm_quad;

	    double a_el(0.0);
	    const std::size_t n_q(dm_quad::n_point);
	    for (std::size_t q(0); q < n_q; ++q) {
	      const std::size_t dim(m.get_embedding_space_dimension());
	      const std::size_t n_dof(fe::n_dof_per_element);

	      const double volume(dm.get_cell_volume(k));
	      const array<double> jmt(dm.get_jmt(k));

	      array<double> xq{n_q, dm_quad::n_point_space_dimension};
	      xq.set_data(&dm_quad::x[0][0]);
	      array<double> omega{n_q};
	      omega.set_data(&dm_quad::w[0]);

	      // Here we must transform the quadrature points on one
	      // of the possible boundary of the reference cell:
	  
	      array<double> hat_xq(fe::cell_type::map_points_to_subdomain(dm.get_subdomain_id(k), xq));


	      array<double> phi{n_dof};
	      for (std::size_t i(0); i < n_dof; ++i)
		phi.at(i) = fe::phi(i, &hat_xq.at(q, 0));

	      array<double> dphi{dim, n_dof};
	      for (std::size_t i(0); i < n_dof; ++i) {
		for (std::size_t n(0); n < dim; ++n) {
		  dphi.at(n, i) = 0.0;
		  for (std::size_t k(0); k < dim; ++k)
		    dphi.at(n, i) += jmt.at(n, k) * fe::dphi(k, i, &hat_xq.at(q, 0));
		}
	      }

	      a_el += 0.0;
	    }
	    sys.accumulate(fes.get_dof(dm.get_parent_element_id(k), j),
			   fes.get_dof(dm.get_parent_element_id(k), i), a_el);
	  }
	}
      }
    }



    array<double> rhs{fes.get_dof_number()};

    // assemble contribution of the volume to the rhs:
    {
      const std::size_t dim(m.get_embedding_space_dimension());
      for (unsigned int k(0); k < m.get_element_number(); ++k) {
	// prepare the quadrature points and weights
	typedef quad::edge::gauss3 m_quad;
	const std::size_t n_q(m_quad::n_point);
	array<double> xq{n_q, m_quad::n_point_space_dimension};
	xq.set_data(&m_quad::x[0][0]);
	array<double> omega{n_q};
	omega.set_data(&m_quad::w[0]);

	// prepare the basis functions on the quadrature points
	const std::size_t n_dof(fe::n_dof_per_element);
	array<double> phi{n_dof, n_q};
	for (unsigned int q(0); q < n_q; ++q) {
	  for (std::size_t i(0); i < n_dof; ++i)
	    phi.at(i, q) = fe::phi(i, &xq.at(q, 0));
	}

	// prepare the basis functions derivatives on the quadrature points
	const array<double> jmt(m.get_jmt(k));
	array<double> dphi{dim, n_dof, n_q};
	for (unsigned int q(0); q < n_q; ++q) {
	  for (std::size_t i(0); i < n_dof; ++i) {
	    for (std::size_t n(0); n < dim; ++n) {
	      dphi.at(n, i, q) = 0.0;
	      for (std::size_t k(0); k < dim; ++k)
		dphi.at(n, i, q) += jmt.at(n, k) * fe::dphi(k, i, &xq.at(q, 0));
	    }
	  }
	}

	// evaluate the weak form
	const double volume(m.get_cell_volume(k));
	for (unsigned int j(0); j < n_dof; ++j) {	
	  double rhs_el(0.0);
	  for (unsigned int q(0); q < n_q; ++q) {
	    //rhs_el += volume * omega.at(q) * force(x) * phi.at(j);
	    rhs_el += volume * omega.at(q) * 1.0 * phi.at(j, q);
	  }
	  rhs.at(fes.get_dof(k, j)) += rhs_el;
	}
      }
    }
    std::cout << std::setw(40) << "linear form volume integration: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;

    // assemble contribution of the right boundary to the rhs:
    {
      const std::size_t dim(m.get_embedding_space_dimension());
      for (std::size_t k(0); k < right_boundary.get_element_number(); ++k) {
	// prepare the quadrature points and weights
	typedef quad::point::eval dm_quad;
	const std::size_t n_q(dm_quad::n_point);
	array<double> xq{n_q, dm_quad::n_point_space_dimension};
	xq.set_data(&dm_quad::x[0][0]);
	array<double> omega{n_q};
	omega.set_data(&dm_quad::w[0]);

	// specific to the boundary
	array<double> hat_xq(fe::cell_type::map_points_to_subdomain(right_boundary.get_subdomain_id(k), xq));

	// prepare the basis functions on the quadrature points
	const std::size_t n_dof(fe::n_dof_per_element);
	array<double> phi{n_dof, n_q};
	for (std::size_t q(0); q < n_q; ++q) {
	  for (std::size_t i(0); i < n_dof; ++i)
	    phi.at(i, q) = fe::phi(i, &hat_xq.at(q, 0));
	}

	// prepare the basis functions derivative on the quadrature points
	const array<double> jmt(right_boundary.get_jmt(k));
	array<double> dphi{dim, n_dof, n_q};
	for (std::size_t q(0); q < n_q; ++q) {
	  for (std::size_t i(0); i < n_dof; ++i) {
	    for (std::size_t n(0); n < dim; ++n) {
	      dphi.at(n, i, q) = 0.0;
	      for (std::size_t k(0); k < dim; ++k)
		dphi.at(n, i, q) += jmt.at(n, k) * fe::dphi(k, i, &hat_xq.at(q, 0));
	    }
	  }
	}

	// evaluate the weak form
	const double volume(right_boundary.get_cell_volume(k));
	for (std::size_t j(0); j < fe::n_dof_per_element; ++j) {
	  double rhs_el(0.0);
	  for (std::size_t q(0); q < n_q; ++q) {
	    // rhs_el += volume * omega.at(q) * g_dot_n(x) * phi.at(j);
	    rhs_el += volume * omega.at(q) * 1.0 * phi.at(j, q);
	  }
	  // specific to the boundary
	  rhs.at(fes.get_dof(right_boundary.get_parent_element_id(k), j)) += rhs_el;
	}
      }
    }
    std::cout << std::setw(40) << "linear form boundary integration: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;

    return 0;
    std::cout << "rhs = [" << rhs.at(0);
    for (std::size_t j(1); j < rhs.get_size(0); ++j)
      std::cout << "; " << rhs.at(j);
    std::cout << "];" << std::endl;
  
    sys.show(std::cout);
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
    
  return 0;
}

