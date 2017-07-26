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
 *  Type list utility tail: return the tail of the list
 */
template<typename TL>
struct tail;

template<typename T, typename ... Ts>
struct tail<type_list<T, Ts...> >: is_type<type_list<Ts...> > {};

template<typename TL>
using tail_t = typename tail<TL>::type;


/*
 *  Type list utility head: return the head of the list
 */
template<typename TL>
struct head;

template<typename T, typename ... Ts>
struct head<type_list<T, Ts...> >: is_type<T> {};

template<typename TL>
using head_t = typename head<TL>::type;


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


/*
 *  Represent a list of compile-time integer values
 */
template<typename T, T ... ts>
struct integral_sequence {};


/*
 *  Concatenate two integral_sequences
 */
template<typename IS1, typename IS2>
struct concat_integral_sequence;

template<typename T, T ... ts1, T ... ts2>
struct concat_integral_sequence<integral_sequence<T, ts1...>,
				integral_sequence<T, ts2...> >
  : is_type<integral_sequence<T, ts1..., ts2...> > {};


/*
 *  Generate an integral sequence of length n, starting from 0
 */
template<typename IC>
struct make_integral_sequence;

template<typename T, T n>
struct make_integral_sequence<integral_constant<T, n> > {
  using type = typename concat_integral_sequence<typename make_integral_sequence<integral_constant<T, n-1> >::type,
						 integral_sequence<T, n - 1> >::type;
};

template<typename T>
struct make_integral_sequence<integral_constant<T, 0> >: is_type<integral_sequence<T> > {};

template<typename T, T n>
using make_integral_sequence_t = typename make_integral_sequence<integral_constant<T, n> >::type;



/*
 *  Convert an integral_sequence to a list of integral_constant's
 */
template<typename IS>
struct sequence_to_list;

template<typename T, T t, T ... ts>
struct sequence_to_list<integral_sequence<T, t, ts...> >
  : is_type< append_t<integral_constant<T, t>,
		      typename sequence_to_list<integral_sequence<T, ts...> >::type> > {};

template<typename T, T t>
struct sequence_to_list<integral_sequence<T, t> >
  : is_type<type_list<integral_constant<T, t> > > {};

template<typename IS>
using sequence_to_list_t = typename sequence_to_list<IS>::type;


template<typename T, T n>
using make_integral_list_t = sequence_to_list_t<make_integral_sequence_t<T, n> >;

/*
 *  Wrap all types in a type_list with the template T
 */
template<template<typename...> class T>
struct wrap_impl {
  template<typename U>
  struct apply {
    using type = T<U>;
  };
};

template<template<typename...> class T, typename TL>
using wrap_t = transform<wrap_impl<T>, TL>;



/*
 *  Append an integral constant IC to every element (which are type_lists) of the list ICLL
 */
template<typename IC>
struct append_to_each_element_impl {
  template<typename ICL, typename ICLL>
  struct apply {
    using type = append_t<append_t<IC, ICL>, ICLL>;
  };
};

template<typename IC, typename ICLL>
using append_to_each_element_t = typename foldr<append_to_each_element_impl<IC>,type_list<>, ICLL>::type;


/*
 * Flatten a list of lists
 */
struct flatten_impl {
  template<typename L1, typename L2>
  struct apply {
    using type = cat_list_t<L1, L2>;
  };
};

template<typename L>
using flatten_list_t = typename foldr<flatten_impl, type_list<>, L>::type;

/*
 *  Generate the tensor product of two list of integral_constants
 */
template<typename IL>
struct tensor_product_impl {
  template<typename IC, typename ICLL>
  struct apply {
    using type = append_t<append_to_each_element_t<IC, wrap_t<type_list, IL> >, ICLL>;
  };
};

template<typename IL1, typename IL2>
using tensor_product_of_lists_t = flatten_list_t<typename foldr<tensor_product_impl<IL1>, type_list<>, IL2>::type>;


/*
 *  Helper function to fill an array with the return value of a set of calls to template functions
 */
template<typename R, typename F, std::size_t n, std::size_t n_max>
struct fill_array_with_return_values {
  template<typename ... As>
  static void fill(R* ptr, As... as) {
    *ptr = F::template call<n>(as...);
    fill_array_with_return_values<R, F, n + 1, n_max>::template fill<As...>(ptr + 1, as...);
  }
};

template<typename R, typename F, std::size_t n_max>
struct fill_array_with_return_values<R, F, n_max, n_max> {
  template<typename ... As>
  static void fill(R* ptr, As... as) {}
};

/*
 *  Helper function to access the n-th value of a function parameter pack
 */
template<std::size_t n, typename T, typename ... Ts>
struct get_arg {
  using return_type = get_element_at_t<n - 1, type_list<Ts...> >;
  static return_type call(T t, Ts ... ts) { return get_arg<n - 1, Ts...>::call(ts...); }
};

template<typename T, typename ... Ts>
struct get_arg<0, T, Ts...> {
  static T call(T t, Ts ... ts) { return t; }
};

template<std::size_t n, typename ... Ts>
get_element_at_t<n, type_list<Ts...> >
argument(Ts ... ts) {
  return get_arg<n, Ts...>::call(ts...);
}

#endif /* _META_H_ */
