
#include <type_traits>
#include <cstddef>
#include <iostream>
#include <typeinfo>


#include "../src/mesh.hpp"
#include "../src/fe.hpp"
#include "../src/fes.hpp"
#include "../src/meta.hpp"


/*
 *  Metafunction: return the cell_type of a fe_type
 */
struct get_cell_type {
  template<typename fe_type>
  struct apply: is_type<typename fe_type::cell_type> {};
};


/*
 *  A composite finite element is simply a wrapper around a typelist
 *  with an interface to access each of the finite element types in the list.
 */
template<typename ... fe_pack>
struct composite_finite_element {
  using fe_list = type_list<fe_pack...>;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;
};


/*
 *  Declaration of a composite finite element space.
 *  A composite finite element space is a wrapper around an std::tuple
 *  where each component is a simple finite element space.
 */
template<typename cfe_type>
class composite_finite_element_space;

template<typename ... fe_pack>
class composite_finite_element_space<composite_finite_element<fe_pack...> > {
public:
  using cfe_type = composite_finite_element<fe_pack...>;
  using cell_list = unique_t<transform<get_cell_type, typename cfe_type::fe_list> >;
  using cell_type = get_element_at_t<0, cell_list>;

  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-value"

  composite_finite_element_space(const mesh<cell_type>& m)
    : fe_instances((sizeof(fe_pack), m)...) {}
  
  #pragma clang diagnostic pop
  
  template<std::size_t n>
  std::size_t get_dof_number() const {
    return std::get<n>(fe_instances).get_dof_number();
  }

  template<std::size_t n>
  unsigned int get_dof(std::size_t k, std::size_t i) const {
    return std::get<n>(fe_instances).get_dof();
  }

  template<std::size_t n>
  const std::unordered_set<unsigned int>& get_dirichlet_dof() const {
    return std::get<n>(fe_instances).get_dirichlet_dof();
  }
	
  template<std::size_t n>
  const std::vector<std::set<cell::subdomain_type> > get_subdomain_list() const {
    return std::get<n>(fe_instances).get_subdomain_list();
  }

  template<std::size_t n>
  void show(std::ostream& stream) {
    std::get<n>(fe_instances).show(stream);
  }
    
  template<std::size_t n>
  const mesh<cell_type>& get_mesh() const {
    return std::get<n>(fe_instances).get_mesh();
  }

  template<std::size_t n>
  double boundary_value(const double* x) const {
    return std::get<n>(fe_instances).boundary_value(x);
  }

  template<std::size_t n>
  array<double> get_dof_space_coordinate(unsigned int i) const {
    return std::get<n>(fe_instances).get_dof_space_coordinate(i);
  }

private:
  std::tuple<finite_element_space<fe_pack>...> fe_instances;
};


void test_meta() {
  {
    using tl = type_list<double, char, int>;
    using tlt = transform<identity, tl>;

    using tl_element_1 = get_element_at_t<0, tl>;
    using tl_element_2 = get_element_at_t<1, tl>;
    using tl_element_3 = get_element_at_t<2, tl>;
    //using tl_element_3 = get_t<3, tl>; // will not compile

    print_type<tl_element_1>("element 0: ");
    print_type<tl_element_2>("element 1: ");
    print_type<tl_element_3>("element 2: ");

    print_type<tlt>("transformed by the identity: ");
    
    std::cout << is_member<int, tl>::value << std::endl;
    std::cout << is_member<double, tl>::value << std::endl;
    std::cout << is_member<char, tl>::value << std::endl;
    std::cout << is_member<short, tl>::value << std::endl;
  }

  {
    using tl = type_list<double, char, double>;

    print_type<unique_t<tl> >();
  }
}


void test_cfe() {
  using cell_type = cell::triangle;
  using fe_0_type = finite_element::triangle_lagrange_p0;
  using fe_1_type = finite_element::triangle_lagrange_p1;
  using fe_type = composite_finite_element<fe_0_type, fe_1_type>;
  using fes_type = composite_finite_element_space<fe_type>;
  
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 3, 3));
  fes_type fes(m);

  fes.show<0>(std::cout);
  fes.show<1>(std::cout);
}


void test_cfes() {
  
}


int main(int argc, char *argv[]) {

  test_meta();
  test_cfe();
  test_cfes();
  
  return 0;
}
