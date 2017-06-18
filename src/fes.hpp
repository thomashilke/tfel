#ifndef _FES_H_
#define _FES_H_

#include "cell.hpp"
#include "mesh.hpp"

template<typename fe>
class finite_element_space {
public:
  typedef fe fe_type;
  typedef typename fe_type::cell_type cell_type;

  finite_element_space(const mesh<cell_type>& m)
    : dof_map{m.get_element_number(),
              fe_type::n_dof_per_element} {
    using cell::subdomain_type;

    std::size_t global_dof_offset(0);
    std::size_t local_dof_offset(0);

    /*
     * Take care of the dofs on the vertices
     */
    for (unsigned int sd(0); sd < cell_type::n_subdomain_type; ++sd) {
      if(fe_type::n_dof_per_subdomain(sd)) {
	std::set<subdomain_type> subdomains(cell_type::get_subdomain_list(m, sd));
	
	const std::size_t n(subdomains.size());
	const std::size_t hat_m(fe_type::n_dof_per_subdomain(sd));
	const std::size_t hat_n(cell_type::n_subdomain(sd));

	for (unsigned int k(0); k < m.get_element_number(); ++k) {
	  for (unsigned int hat_j(0); hat_j < hat_n; ++hat_j) {
	    // for subdomain hat_j:
	    subdomain_type subdomain(cell_type::get_subdomain(m, k, sd, hat_j));
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
  }


  void show(std::ostream& stream) {
    for (unsigned int k(0); k < dof_map.get_size(0); ++k) {
      for (unsigned int i(0); i < dof_map.get_size(1); ++i) {
	stream << dof_map.at(k, i) << " ";
      }
      stream << std::endl;
    }
  }

private:
  array<unsigned int> dof_map;
};

#endif /* _FES_H_ */
