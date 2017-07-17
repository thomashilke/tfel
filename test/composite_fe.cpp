
#include <type_traits>
#include <cstddef>
#include <iostream>
#include <typeinfo>

#include <cxxabi.h>

#include "../src/mesh.hpp"
#include "../src/fe.hpp"
#include "../src/fes.hpp"

template<typename T>
void print_type(std::string msg = "") {
  int status(0);
  char *name(abi::__cxa_demangle(typeid(T).name(), 0, 0, &status));
  
  std::cout << msg << name << std::endl;

  delete name;
  name = nullptr;
}

/*
template<typename T, T v>
struct integral_constant { static constexpr T value = v; };

using true_type = bool_constant<true>;
using false_type = bool_constant<false>;
*/


/*
 *  TMP utilities
 */
using std::integral_constant;
using std::true_type;
using std::false_type;

template<bool b>
using bool_constant = integral_constant<bool, b>;

template<typename T>
struct is_type { typedef T type; };


/*
 *  Type list definition
 */
template<typename ... Ts> struct type_list {};


/*
 * Type list utility size: return the number of elements in the list
 */
template<typename TL>
struct list_size;

template<typename ... Ts>
struct list_size<type_list<Ts...> >: integral_constant<std::size_t, sizeof...(Ts)> {};

/*
 *  Type list utility append: append an element in front of a type list
 */
template<typename T, typename TL>
struct append;

template<typename T, typename ... Ts>
struct append<T, type_list<Ts...> >: is_type<type_list<T, Ts...> > {};

template<typename T, typename TL>
using append_t = typename append<T, TL>::type;
  
/*
 *  Type list utility get_element_at: access element of a type list
 */
template<std::size_t n, typename ... T>
struct get_element_at;

template<typename T, typename ... Ts>
struct get_element_at<0, type_list<T, Ts...> >: is_type<T> {};

template<std::size_t n, typename T, typename ... Ts>
struct get_element_at<n, type_list<T, Ts...> >: public get_element_at<n - 1, type_list<Ts...> > {};

template<std::size_t n, typename ... Ts>
using get_element_at_t = typename get_element_at<n, Ts...>::type;


/*
 *  Type list utility is_member: check type membership of a type list
 */
template<typename T, typename TL>
struct is_member;

template<typename T, typename ... Ts>
struct is_member<T, type_list<T, Ts...> >: true_type {};

template<typename T, typename U, typename ... Ts>
struct is_member<T, type_list<U, Ts...> >: is_member<T, type_list<Ts...> > {};

template<typename T>
struct is_member<T, type_list<> >: false_type {};


/*
 * Foldr meta algorithm
 */

template<typename F, typename S, typename L>
struct foldr;

template<typename F, typename S, typename T, typename ... Ts>
struct foldr<F, S, type_list<T, Ts...> >: is_type<typename F::template apply<T, typename foldr<F, S, type_list<Ts...> >::type >::type > {};

template<typename F, typename S>
struct foldr<F, S, type_list<> >: is_type<S> {};

template<typename F, typename S, typename L>
using foldr_t = typename foldr<F, S, L>::type;


/*
 * Foldl meta algorithm
 */

template<typename F, typename S, typename L>
struct foldl;

template<typename F, typename S, typename T, typename ... Ts>
struct foldl<F, S, type_list<T, Ts...> >: is_type<typename foldl<F, typename F::template apply<S, T>::type, type_list<Ts...> >::type> {};

template<typename F, typename S, typename T>
struct foldl<F, S, type_list<T> >: is_type<typename F::template apply<S, T>::type> {};

template<typename F, typename S, typename L>
using foldl_t = typename foldl<F, S, L>::type;


/*
 * Type list utility unique: generate a list of unique type from a type list
 */
namespace unique_impl {
  struct append_if_absent {
    template<typename T, typename TL>
    struct apply: std::conditional<is_member<T, TL>::value,
				   is_type<TL>,
				   is_type<append_t<T, TL> > >::type {};
  };
}

template<typename TL>
using unique_t = typename foldr<unique_impl::append_if_absent, type_list<>, TL>::type;


template<typename F>
struct transform_impl {
  template<typename A, typename TL>
  struct apply: is_type<append_t<typename F::template apply<A>::type, TL> > {};
};

template<typename F, typename TL>
using transform = typename foldr<transform_impl<F>, type_list<>, TL>::type;

struct identity {
  template<typename T>
  struct apply: is_type<T> {};
};



template<typename ... fe_pack>
struct composite_finite_element {
  using fe_list = type_list<fe_pack...>;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;
};

struct get_cell_type {
  template<typename fe_type>
  struct apply: is_type<typename fe_type::cell_type> {};
};

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
  
  return 0;
}
