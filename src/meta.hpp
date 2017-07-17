#ifndef _META_H_
#define _META_H_

#include <cxxabi.h>

template<typename T>
void print_type(std::string msg = "") {
  int status(0);
  char *name(abi::__cxa_demangle(typeid(T).name(), 0, 0, &status));
  
  std::cout << msg << name << std::endl;

  delete name;
  name = nullptr;
}

template<typename T>
void print_type(const T&, std::string msg = "") {
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
 *  Type list utility cat_list: concatenate two lists
 */
template<typename TL1, typename TL2>
struct cat_list;

template<typename ... Ts1, typename ... Ts2>
struct cat_list<type_list<Ts1...>, type_list<Ts2...> >: is_type<type_list<Ts1..., Ts2...> > {};

template<typename TL1, typename TL2>
using cat_list_t = typename cat_list<TL1, TL2>::type;


/*
 *  Type list utility get_element_at: access element of a type list
 */
template<std::size_t n, typename ... T>
struct get_element_at;

template<typename T, typename ... Ts>
struct get_element_at<0, type_list<T, Ts...> >: is_type<T> {};

template<std::size_t n, typename T, typename ... Ts>
struct get_element_at<n, type_list<T, Ts...> >: get_element_at<n - 1, type_list<Ts...> > {};

template<std::size_t n, typename ... Ts>
using get_element_at_t = typename get_element_at<n, Ts...>::type;


/*
 *  Type list utility get_index_of_element: return the index of the element in the type list TL
 */
template<std::size_t n, typename T, typename TL>
struct get_index_of_element_impl;

template<std::size_t n, typename T, typename U, typename ... Ts>
struct get_index_of_element_impl<n, T, type_list<U, Ts...> >: get_index_of_element_impl<n + 1, T, type_list<Ts...> > {};

template<std::size_t n, typename T, typename ... Ts>
struct get_index_of_element_impl<n, T, type_list<T, Ts...> >: integral_constant<std::size_t, n> {};

template<typename T, typename TL>
using get_index_of_element = get_index_of_element_impl<0, T, TL>;


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
 * Type list utility reverse_list: return the list with the elements in reverse order
 */
struct append_to_list {
  template<typename TL, typename A>
  struct apply: is_type<append_t<A, TL> > {};
};

template<typename TL>
using reverse_list_t = typename foldl<append_to_list, type_list<>, TL>::type;


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


/*
 * Type list utility transform: return a list which items are 
 * the result of F applied to each element of TL
 */
template<typename F>
struct transform_impl {
  template<typename A, typename TL>
  struct apply: is_type<append_t<typename F::template apply<A>::type, TL> > {};
};

template<typename F, typename TL>
using transform = typename foldr<transform_impl<F>, type_list<>, TL>::type;


/*
 * Metafunction: the identity metafunction
 */
struct identity {
  template<typename T>
  struct apply: is_type<T> {};
};


/*
 * Build a list of increasing indices of the same length
 * as the type list TL
 */
struct enumerate_list_elements {
  template<typename A, typename IL>
  struct apply {
    using type =
      append_t<integral_constant<std::size_t, list_size<IL>::value>,
	       IL>;
  };
};

template<typename TL>
using make_index_list_t = reverse_list_t<typename foldr<enumerate_list_elements,
							type_list<>,
							TL>::type>;

#endif /* _META_H_ */
