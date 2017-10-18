#ifndef _FES_H_
#define _FES_H_

#include <unordered_set>
#include <ostream>

#include "cell.hpp"
#include "mesh.hpp"
#include "mesh_data.hpp"

template<typename fe>
class finite_element_space {
public:
  typedef fe fe_type;
  typedef typename fe_type::cell_type cell_type;
  struct element;

  finite_element_space(const fe_mesh<cell_type>& m)
    : m(m),
      dof_map{m.get_cell_number(),
      fe_type::n_dof_per_element},
      global_dof_to_local_dof{0},
      f_bc(default_f_bc) {
    const array<unsigned int>& elements(m.get_cells());
    
    using cell::subdomain_type;

    std::size_t global_dof_offset(0);
    std::size_t local_dof_offset(0);

    /*
     * Take care of the dofs on the different subdomains
     */
    for (unsigned int sd(0); sd < cell_type::n_subdomain_type; ++sd) {
      if(fe_type::n_dof_per_subdomain(sd)) {
	subdomain_list.push_back(cell_type::get_subdomain_list(elements, sd));

	// Copy the set into a vector for fast access
	std::vector<subdomain_type> subdomains(subdomain_list.back().size());
	std::copy(subdomain_list.back().begin(), subdomain_list.back().end(),
		  subdomains.begin());
	
	const std::size_t n(subdomains.size());
	const std::size_t hat_m(fe_type::n_dof_per_subdomain(sd));
	const std::size_t hat_n(cell_type::n_subdomain(sd));

	for (unsigned int k(0); k < m.get_cell_number(); ++k) {
	  for (unsigned int hat_j(0); hat_j < hat_n; ++hat_j) {
	    // for subdomain hat_j:
	    subdomain_type subdomain(cell_type::get_subdomain(elements, k, sd, hat_j));
	    // j is the global index of subdomain hat_j
	    const std::size_t j(std::distance(subdomains.begin(),
					      std::lower_bound(subdomains.begin(),
							       subdomains.end(),
							       subdomain)));

	    for (unsigned int hat_i(0); hat_i < hat_m; ++hat_i) {
	      // for each local dof hat_i of subdomain hat_j
	      dof_map.at(k, (hat_j * hat_m + hat_i) + local_dof_offset)
		= (j * hat_m + hat_i) + global_dof_offset;
	    }
	  }
	}
	global_dof_offset += hat_m * n;
	local_dof_offset += hat_m * hat_n;
      }
    }

    dof_number = global_dof_offset;

    global_dof_to_local_dof = array<unsigned int>{dof_number, 2};
    for (std::size_t k(0); k < dof_map.get_size(0); ++k)
      for (std::size_t n(0); n < dof_map.get_size(1); ++n) {
	global_dof_to_local_dof.at(dof_map.at(k, n), 0) = k;
	global_dof_to_local_dof.at(dof_map.at(k, n), 1) = n;
      }
  }

  template<typename c_cell_type>
  finite_element_space(const fe_mesh<cell_type>& m, const submesh<cell_type, c_cell_type>& dm)
    : finite_element_space(m) {
    f_bc = default_f_bc;
    set_dirichlet_boundary(dm);
  }

  finite_element_space(const fe_mesh<cell_type>& m,
		       const submesh<cell_type>& dm,
		       double (*f_bc)(const double*)): finite_element_space(m, dm) {
    this->f_bc = std::function<double(const double*)>(f_bc);
  }

  template<typename c_cell_type>
  void set_dirichlet_boundary(const submesh<cell_type, c_cell_type>& dm) {
    using cell::subdomain_type;

    dirichlet_dof.clear();
    
    std::size_t global_dof_offset(0);
    for (std::size_t sd(0); sd < submesh<cell_type>::cell_type::n_subdomain_type; ++sd) {
      const std::size_t hat_m(fe_type::n_dof_per_subdomain(sd));
      if (fe_type::n_dof_per_subdomain(sd)) {
	const array<unsigned int>& elements(dm.get_cells());
	std::set<subdomain_type> subdomains(submesh<cell_type>::cell_type::get_subdomain_list(elements, sd));
	for (const auto& subdomain: subdomains) {
	  const std::size_t j(std::distance(subdomain_list[sd].begin(),
					    subdomain_list[sd].find(subdomain)));
	  for (unsigned int hat_i(0); hat_i < hat_m; ++hat_i)
	    dirichlet_dof.insert((j * hat_m + hat_i) + global_dof_offset);
	}
	global_dof_offset += hat_m * subdomain_list[sd].size();
      }
    }    
  }

  void set_dirichlet_condition(const std::function<double(const double*)>& f_bc) {
    this->f_bc = f_bc;
  }

  template<typename c_cell_type>
  void set_dirichlet_boundary_condition(const submesh<cell_type, c_cell_type>& dm,
					const std::function<double(const double*)>& f_bc) {
    this->set_dirichlet_condition(f_bc);
    this->set_dirichlet_boundary(dm);
  }
  
  std::size_t get_dof_number() const { return dof_number; }

  unsigned int get_dof(std::size_t k, std::size_t i) const {
    return dof_map.at(k, i);
  }

  const std::unordered_set<unsigned int>& get_dirichlet_dof() const { return dirichlet_dof; }
  const std::vector<std::set<cell::subdomain_type> > get_subdomain_list() const { return subdomain_list; }
  
  void show(std::ostream& stream) {
    for (unsigned int k(0); k < dof_map.get_size(0); ++k) {
      stream << "element " << k << ": ";
      for (unsigned int i(0); i < dof_map.get_size(1); ++i) {
	stream << dof_map.at(k, i) << " ";
      }
      stream << std::endl;
    }

    stream << "dirichlet dofs: ";
    for (const auto& dof: dirichlet_dof)
      stream << dof << " ";
    stream << std::endl;

    stream << "dof map:" << std::endl;
    for (std::size_t k(0); k < m.get_cell_number(); ++k) {
      for (std::size_t n(0); n < fe_type::n_dof_per_element; ++n)
	stream << "local dof (k = " << k << ", n = " << n << ") -> " << get_dof(k, n) << std::endl;
    }
  }

  const fe_mesh<cell_type>& get_mesh() const { return m; }

  double boundary_value(const double* x) const { return f_bc(x); }

  array<double> get_dof_space_coordinate(unsigned int i) const {
    const std::size_t local_node_id(global_dof_to_local_dof.at(i, 1));
    
    array<double> x{1, cell_type::n_dimension};
    std::copy(&fe_type::x[local_node_id][0],
	      &fe_type::x[local_node_id][0] + cell_type::n_dimension,
	      &x.at(0, 0));

    return cell_type::map_points_to_space_coordinates(m.get_vertices(),
						      m.get_cells(),
						      global_dof_to_local_dof.at(i,0),
						      x);
  }

  std::size_t get_dof_element(std::size_t i) const { return global_dof_to_local_dof.at(i, 0); }
  std::size_t get_dof_local_id(std::size_t i) const { return global_dof_to_local_dof.at(i, 1); }
  
private:
  const fe_mesh<cell_type>& m;
  array<unsigned int> dof_map;
  array<unsigned int> global_dof_to_local_dof;
  std::size_t dof_number;
  
    
  std::unordered_set<unsigned int> dirichlet_dof;
  std::vector<std::set<cell::subdomain_type> > subdomain_list;

  static double default_f_bc(const double* x) { return 0.0; }
  std::function<double(const double*)> f_bc;
};


template<typename fe>
struct finite_element_space<fe>::element {
public:
  using fe_type = fe;
  using cell_type = typename fe_type::cell_type;
  
  element(const finite_element_space<fe>& fes)
    : coefficients{fes.get_dof_number()}, fes(fes) {}
  
  element(const finite_element_space<fe>& fes,
	  array<double>&& a)
    : coefficients(a), fes(fes) {}
  
  element(const finite_element_space<fe>& fes,
	  const array<double>& a)
    : coefficients(a), fes(fes) {}
  
  element(const element& e)
    : coefficients(e.coefficients), fes(e.fes) {}
  
  ~element() {}

  const finite_element_space<fe>& get_finite_element_space() const { return fes; }
  
  const fe_mesh<typename fe::cell_type>& get_mesh() const { return fes.get_mesh(); }

  const array<double>& get_coefficients() const { return coefficients; }

  double evaluate(std::size_t k, const double* x_hat) const {
    typedef fe fe_type;
    
    double value(0.0);
    for (std::size_t n(0); n < fe_type::n_dof_per_element; ++n) {
      value += coefficients.at(fes.get_dof(k, n)) * fe_type::phi(n, x_hat);
    }
    return value;
  }

  double evaluate(const double* x) const {
    std::size_t k(get_mesh().get_cell_at(x));
    
    const array<double>
      bc_coord(cell_type::get_barycentric_coordinates(
        get_mesh().get_vertices(),
	get_mesh().get_cells(),
	k, x));

    return evaluate(k, &bc_coord.at(1));
  }

  typename finite_element_space<fe>::element& operator=(const typename finite_element_space<fe>::element& e) {
    if (&fes != &(e.fes))
      throw std::string("Assigment of elements between different finite element spaces is not supported.");
    coefficients = e.coefficients;
    return *this;
  }

  typename finite_element_space<fe>::element
  restrict(const finite_element_space<fe>& sm_fes, const submesh<cell_type, cell_type>& sm) const {
    array<double> sm_coefficients{sm_fes.get_dof_number()};
    sm_coefficients.fill(0.0);

    for (std::size_t k(0); k < sm.get_cell_number(); ++k)
      for (std::size_t i(0); i < fe::n_dof_per_element; ++i)
	sm_coefficients.at(sm_fes.get_dof(k, i)) = this->coefficients.at(fes.get_dof(sm.get_parent_cell_id(k), i));
    
    return typename finite_element_space<fe>::element(sm_fes, sm_coefficients);
  }

  typename finite_element_space<fe>::element
  extend(const finite_element_space<fe>& m_fes, const submesh<cell_type, cell_type>& sm) const {
    array<double> m_coefficients{m_fes.get_dof_number()};
    m_coefficients.fill(0.0);

    for (std::size_t k(0); k < fes.get_mesh().get_cell_number(); ++k)
      for (std::size_t i(0); i < fe::n_dof_per_element; ++i)
	m_coefficients.at(m_fes.get_dof(sm.get_parent_cell_id(k), i)) = this->coefficients.at(fes.get_dof(k, i));

    return typename finite_element_space<fe>::element(m_fes, m_coefficients);
  }

  typename finite_element_space<fe>::element&
  transform(const std::function<double(double)>& f) {
    for (std::size_t i(0); i < coefficients.get_size(0); ++i)
      coefficients.at(i) = f(coefficients.at(i));

    return *this;
  }

  typename finite_element_space<fe>::element
  transformed(const std::function<double(double)>& f) {
    finite_element_space<fe>::element result(*this);
    result.transform(f);
    return result;
  }

  typename finite_element_space<fe>::element&
  operator*=(double s) {
    for (std::size_t i(0); i < coefficients.get_size(0); ++i)
      coefficients.at(i) *= s;

    return *this;
  }

  typename finite_element_space<fe>::element&
  operator+=(const finite_element_space<fe>::element& op) {
    for (std::size_t i(0); i < coefficients.get_size(0); ++i)
      coefficients.at(i) += op.coefficients.at(i);

    return *this;
  }

  typename finite_element_space<fe>::element&
  operator/=(double s) {
    for (std::size_t i(0); i < coefficients.get_size(0); ++i)
      coefficients.at(i) /= s;
    
    return *this;
  }

  typename finite_element_space<fe>::element&
  operator-=(const finite_element_space<fe>::element& op) {
    for (std::size_t i(0); i < coefficients.get_size(0); ++i)
      coefficients.at(i) -= op.coefficients.at(i);
    
    return *this;
  }
  
private:
  array<double> coefficients;
  const finite_element_space<fe>& fes;
};

template<typename fe>
typename finite_element_space<fe>::element
operator+(const typename finite_element_space<fe>::element& op1,
          const typename finite_element_space<fe>::element& op2) {
  typename finite_element_space<fe>::element result(op1);
  return op1 += op2;
}

template<typename fe>
typename finite_element_space<fe>::element
operator-(const typename finite_element_space<fe>::element& op1,
          const typename finite_element_space<fe>::element& op2) {
  typename finite_element_space<fe>::element result(op1);
  return op1 -= op2;
}

template<typename fe>
typename finite_element_space<fe>::element
operator*(double s,
          const typename finite_element_space<fe>::element& op) {
  typename finite_element_space<fe>::element result(op);
  return op *= s;
}

template<typename fe>
typename finite_element_space<fe>::element
operator*(const typename finite_element_space<fe>::element& op,
          double s) {
  return s * op;
}

template<typename fe>
typename finite_element_space<fe>::element
operator/(const typename finite_element_space<fe>::element& op,
          double s) {
  typename finite_element_space<fe>::element result(op);
  return op /= s;
}


template<typename cell_type>
typename finite_element_space<typename cell_type::fe::lagrange_p1>::element
to_p1_finite_element_function(const mesh_data<double, fe_mesh<cell_type> >& data,
                              const finite_element_space<typename cell_type::fe::lagrange_p1>& fes) {
  using return_type = typename finite_element_space<typename cell_type::fe::lagrange_p1>::element;

  if (data.get_kind() != mesh_data_kind::vertex)
    throw std::string("wrong data kind");

  if (data.get_component_number() != 1)
    throw std::string("wrong component number");

  const fe_mesh<cell_type>& m(fes.get_mesh());
  
  array<double> coefficients{fes.get_dof_number()};

  for (std::size_t k(0); k < fes.get_dof_number(); ++k) {
    const std::size_t element_id(fes.get_dof_element(k));
    const std::size_t dof_local_id(fes.get_dof_local_id(k));
    coefficients.at(k) = data.value(m.get_cells().at(element_id, dof_local_id), 0);
  }
  
  return return_type(fes, coefficients);
}

template<typename cell_type, typename expr_type>
typename finite_element_space<typename cell_type::fe::lagrange_p1>::element
to_p1_finite_element_function(const mesh_data_expression<expr_type>& expr,
                              const finite_element_space<typename cell_type::fe::lagrange_p1>& fes) {
  return to_p1_finite_element_function(mesh_data<double, fe_mesh<cell_type> >(fes.get_mesh(),
                                                                           mesh_data_kind::vertex,
                                                                           expr), fes);
}

template<typename cell_type>
typename finite_element_space<typename cell_type::fe::lagrange_p0>::element
to_p0_finite_element_function(const mesh_data<double, fe_mesh<cell_type> >& data,
                              const finite_element_space<typename cell_type::fe::lagrange_p0>& fes) {
  using return_type = typename finite_element_space<typename cell_type::fe::lagrange_p0>::element;

  if (data.get_kind() != mesh_data_kind::cell)
    throw std::string("wrong data kind");

  if (data.get_component_number() != 1)
    throw std::string("wrong component number");

  const fe_mesh<cell_type>& m(fes.get_mesh());
  
  array<double> coefficients{fes.get_dof_number()};

  for (std::size_t k(0); k < fes.get_dof_number(); ++k) {
    const std::size_t element_id(fes.get_dof_element(k));
    coefficients.at(k) = data.values(element_id, 0);
  }
  
  return return_type(fes, coefficients);
}

template<typename fe_type>
mesh_data<double, fe_mesh<typename fe_type::cell_type> >
to_mesh_cell_data(const typename finite_element_space<fe_type>::element& v) {
  using return_type = mesh_data<double, fe_mesh<typename fe_type::cell_type> >;
  using cell_type = typename fe_type::cell_type;

  const finite_element_space<fe_type>& fes(v.get_finite_element_space());
  const fe_mesh<cell_type>& m(fes.get_mesh());
  
  array<double> values{m.get_cell_number(), 1};

  const array<double> bc_hat(cell_type::barycenter());
  for (std::size_t k(0); k < fes.get_dof_number(); ++k) {
    values.at(k, 0) = v.evaluate(k, &bc_hat.at(0));
  }
  
  return return_type(m, mesh_data_kind::cell, values);
}

template<typename fe_type>
mesh_data<double, fe_mesh<typename fe_type::cell_type> >
to_mesh_vertex_data(const typename finite_element_space<fe_type>::element& v) {
  static_assert(is_continuous<fe_type>::value,
		"conversion to vertex data is only defined "
		"for continuous finite element spaces.");
  
  using return_type = mesh_data<double, fe_mesh<typename fe_type::cell_type> >;
  using cell_type = typename fe_type::cell_type;

  const finite_element_space<fe_type>& fes(v.get_finite_element_space());
  const fe_mesh<cell_type>& m(fes.get_mesh());
  
  array<double> values{m.get_vertex_number(), 1};

#warning "FIXME: This works only for finite element with nodes on the cell vertices"
  for (std::size_t k(0); k < m.get_cell_number(); ++k)
    for (std::size_t n(0); n < cell_type::n_vertex_per_cell; ++n)
      values.at(m.get_cells().at(k, n), 0) = v.get_coefficients().at(fes.get_dof(k, n));
  
  return return_type(m, mesh_data_kind::vertex, values);
}

#endif /* _FES_H_ */
