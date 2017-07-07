
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"
#include "mesh.hpp"
#include "fes.hpp"
#include "quadrature.hpp"
#include "expression.hpp"


class timer {
public:
  timer(): start(std::chrono::high_resolution_clock::now()) {}
  double tic() {
    std::chrono::time_point<std::chrono::high_resolution_clock>
      now(std::chrono::high_resolution_clock::now());

    const std::chrono::duration<double> elapsed(now - start);
    return std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  }

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
};


class failed_hypothesis {
public:
  failed_hypothesis(const std::string& msg) {}
};


class sparse_linear_system {
public:
  sparse_linear_system(std::size_t n_eq, std::size_t n_unk)
    :n_equation(n_eq), n_unknown(n_unk) {}
  void accumulate(std::size_t i, std::size_t j, double v) {
    elements[std::pair<std::size_t, std::size_t>(i, j)] += v;
  }

  void show(std::ostream& stream) const {
    array<double> profile{n_equation, n_unknown};
    for (const auto& elem: elements)
      profile.at(elem.first.first, elem.first.second) = elem.second;

    stream << "a = [";
    for (unsigned int i(0); i < n_equation - 1; ++i) {
      for (unsigned int j(0); j < n_unknown - 1; ++j)
	stream << std::setw(4) << std::setprecision(12) << profile.at(i, j) << ", ";
      stream << std::setw(4) << std::setprecision(12) << profile.at(i, n_unknown - 1) << "; " << std::endl;
      
    }
    for (unsigned int j(0); j < n_unknown - 1; ++j)
      stream << std::setw(4) << std::setprecision(12) << profile.at(n_equation - 1, j) << ", ";
    stream << std::setw(4) << std::setprecision(12) << profile.at(n_equation - 1, n_unknown - 1) << "];" << std::endl;
  }
  
private:
  std::size_t n_equation, n_unknown;
  std::map<std::pair<std::size_t, std::size_t>, double> elements;
};


mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n) {
  std::vector<double> vertices(n + 1);
  std::vector<unsigned int> elements(2 * n);

  for (unsigned int i(0); i < n + 1; ++i)
    vertices[i] = x_1 + (x_2 - x_1) / n * i;

  for (unsigned int k(0); k < n; ++k) {
    elements[2 * k] = k;
    elements[2 * k + 1] = k + 1;
  }

  return mesh<cell::edge>(&vertices[0], n + 1, 1,
			  &elements[0], n);
}

double force(double* x) {
  return 1.0;
}



//bilinear_form a(fes, fes);
//a += integrate<quad::triangle::>(d<1>(phi) * d<1>(psi), m);

//linear_form f(fes);
//f += integrate<quad::triangle::>(force * psi, m);
//f += integrate<quad::edge::gauss2>(g * psi, dm);


/*
template<typename cell_type, typename quadrature_type>
struct mesh_integration_proxy {
  const mesh<cell_type>& m;

  std::size_t get_global_element_id(std::size_t k) {
    return k;
  }

  array<double> get_quadrature_points(std::size_t k) {
    const std::size_t n_q(dm_quad::n_point);
    array<double> xq{n_q, dm_quad::n_point_space_dimension};
    xq.set_data(&dm_quad::x[0][0]);
    return hat_xq;
  }
};

template<typename cell_type, typename quadrature_type>
struct submesh_integration_proxy {
  const submesh<cell_type>& m;
  
  std::size_t get_global_element_id(std::size_t k) {
    return m.get_parent_element_id(k);
  }
  
  array<double> get_quadrature_points(std::size_t k) {
    const std::size_t n_q(dm_quad::n_point);
    array<double> xq{n_q, dm_quad::n_point_space_dimension};
    xq.set_data(&dm_quad::x[0][0]);
    array<double> hat_xq(cell_type::map_points_to_subdomain(m.get_subdomain_id(k), xq));
    return hat_xq;
  }
};

template<typename quadrature_type, typename cell_type>
T integrate(const expression<form_type>& , const mesh<cell_type>& m) {
  
}

template<typename quadrature_type, typename cell_type>
T integrate(const expression<form_type>&, const submesh<cell_type>& m) {
  
}

template<typename test_fes_type, typename trial_fes_type>
class bilinear_form {
public:
  bilinear_form(const test_fes_type& te_fes,
		const trial_fes_type& tr_fes)
    : a(te_fes.get_dof_number(), tr_fes.get_dof_number()) {}

  template<typename T>
  void operator+=(const T& integration_proxy) {
    // integration proxy encapsule la forme bilineaire, la quadrature et le maillage/sous-maillage
    // bilinear form encapsule le couple (test space, trial space)
    
  }
  
private:
  sparse_linear_system a;
};

template<typename test_fes_type>
class linear_form {
  linear_form(const test_fes_type& te_fes): f{te_fes.get_dof_number()} {}
  
private:
  array<double> f;
};
*/

int main(int argc, char *argv[]) {
  timer t;
  try {
    const std::size_t n(10000);
    
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
  
    //m.show(std::cout);
    //right_boundary.show(std::cout);
  
    typedef finite_element::edge_lagrange_p1 fe;
    finite_element_space<fe> fes(m, left_boundary);

    //fes.show(std::cout);

    sparse_linear_system sys(fes.get_dof_number(),
			     fes.get_dof_number());

    std::cout << std::setw(40) << "finite element system: "
	      << std::setw(6) << std::right << t.tic() << " [ms]" << std::endl;

    /*
    //trial_function_t phi((form<1,0,0>()));
    //test_function_t psi((form<0,0,0>()));
    auto phi(bilinear_form.trial_function());
    auto psi(bilinear_form.test_function());

    bilinear_form a(fes, fes);
    a += integrate<quad::triangle::>(d<1>(phi) * d<1>(psi), m);

    linear_form f(fes);
    f += integrate<quad::triangle::>(force * psi, m);
    f += integrate<quad::edge::gauss2>(g * psi, dm);

    fes::element u(a.solve(f));
    */

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

