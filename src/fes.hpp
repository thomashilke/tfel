#ifndef _FES_H_
#define _FES_H_

#include <unordered_set>
#include <ostream>

#include "cell.hpp"
#include "mesh.hpp"

template<typename fe>
class finite_element_space {
public:
  typedef fe fe_type;
  typedef typename fe_type::cell_type cell_type;
  struct element;

  finite_element_space(const mesh<cell_type>& m)
    : m(m),
      dof_map{m.get_element_number(),
      fe_type::n_dof_per_element} {
    const array<unsigned int>& elements(m.get_elements());
    
    using cell::subdomain_type;

    std::size_t global_dof_offset(0);
    std::size_t local_dof_offset(0);

    /*
     * Take care of the dofs on the vertices
     */
    for (unsigned int sd(0); sd < cell_type::n_subdomain_type; ++sd) {
      if(fe_type::n_dof_per_subdomain(sd)) {
	subdomain_list.push_back(cell_type::get_subdomain_list(elements, sd));
	std::set<subdomain_type>& subdomains(subdomain_list.back());
	
	const std::size_t n(subdomains.size());
	const std::size_t hat_m(fe_type::n_dof_per_subdomain(sd));
	const std::size_t hat_n(cell_type::n_subdomain(sd));

	for (unsigned int k(0); k < m.get_element_number(); ++k) {
	  for (unsigned int hat_j(0); hat_j < hat_n; ++hat_j) {
	    // for subdomain hat_j:
	    subdomain_type subdomain(cell_type::get_subdomain(elements, k, sd, hat_j));
	    // j is the global index of subdomain hat_j
	    const std::size_t j(std::distance(subdomains.begin(),
					      subdomains.find(subdomain)));
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
  }

  finite_element_space(const mesh<cell_type>& m, const submesh<cell_type>& dm)
    : finite_element_space(m) {

    using cell::subdomain_type;
    
    std::size_t global_dof_offset(0);
    for (std::size_t sd(0); sd < submesh<cell_type>::cell_type::n_subdomain_type; ++sd) {
      const std::size_t hat_m(fe_type::n_dof_per_subdomain(sd));
      if (fe_type::n_dof_per_subdomain(sd)) {
	const array<unsigned int>& elements(dm.get_elements());
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
  }

  const mesh<cell_type>& get_mesh() const { return m; }

private:
  const mesh<cell_type>& m;
  array<unsigned int> dof_map;
  std::size_t dof_number;

  std::unordered_set<unsigned int> dirichlet_dof;
  std::vector<std::set<cell::subdomain_type> > subdomain_list;
};

template<typename fe>
struct finite_element_space<fe>::element {
public:
  element(const finite_element_space<fe>& fes,
	  array<double>&& a): coefficients(a), fes(fes) {}
  element(const finite_element_space<fe>& fes,
	  const array<double>& a): coefficients(a), fes(fes) {}
  element(const element& e): coefficients(e.coefficients), fes(e.fes) {}
  ~element() {}

  const finite_element_space<fe>& get_finite_element_space() const { return fes; }
  const mesh<typename fe::cell_type>& get_mesh() const { return fes.get_mesh(); }

  const array<double>& get_components() const { return coefficients; }

  double evaluate(std::size_t k, const double* x_hat) const {
    typedef fe fe_type;
    
    double value(0.0);
    for (std::size_t n(0); n < fe_type::n_dof_per_element; ++n) {
      value += coefficients.at(fes.get_dof(k, n)) * fe_type::phi(n, x_hat);
    }
    return value;
  }
  
private:
  array<double> coefficients;
  const finite_element_space<fe>& fes;
};

#endif /* _FES_H_ */
