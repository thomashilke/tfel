#ifndef FE_VALUE_MANAGER_H
#define FE_VALUE_MANAGER_H


template<template<typename> F, typename TL, typename ... As>
struct call_for_each;

template<template<typename> F, typename T, typename ... Ts, typename ... As>
struct call_for_each<F, type_list<T, Ts...>, As...> {
  static void call_for_each(F f, As&& ... as) {
    f<T>(std::forward(as)...);
  }
};

template<template<typename> F, typename ... Ts, typename ... As>
struct call_for_each<F, type_list<>, As...> {
  static void call_for_each(F f, As&& ... as) {}
};
  


template<typename FE_LIST>
struct fe_value_manager;

template<typename ... fe_pack>
struct fe_value_manager<type_list<fe_pack...> > {
  using fe_list = type_list<fe_pack...>;
  using cell_list = unique_t<transform<get_cell_type, typename cfe_type::fe_list> >;
  using cell_type = get_element_at_t<0, cell_list>;

  template<std::size_t n>
  using fe_type = get_element_at_t<n, fe_list>;
  
  fe_value_manager(std::size_t n_quadrature_point):
    n_q(n_q),
    values({n_quadrature_point,
	  fe_pack::n_dof_per_element,
	  fe_pack::cell_type::n_dimension + 1}...) {
    static_assert(list_size<cell_list>::value == 1,
		  "the cells types must be homogeneous");
  }

  void prepare(const array<double>& jmt,
	       const array<double>& xq) {
    using index_sequence = make_index_list_t<fe_list>;
    call_for_each<prepare_impl, index_sequence>::call(values, jmt, xq);
  }
  
  template<std::size_t n>
  const array<double>& get_values() const {
    return std::get<n>(values);
  }
  
private:
  template<typename A, typename B> struct return_2nd: is_type<B> {};
  template<typename A, typename B> using return_2nd_t = typename return_2nd<A, B>::type;
  
  const std::size_t n_q;
  std::tuple<return_2nd_t<array<double>, fe_pack>... > values;

  template<typename integer_value>
  static void prepare_impl(std::tuple<array<double> >& values,
			   const array<double>& jmt,
			   const array<double>& xq) {
    const std::size_t n(integer_value::value);
    
    array<double>& phi(std::get<n>(values));
    
    const std::size_t n_dim(cell_type::n_dimension);
    const std::size_t n_dof(fe_type<n>::n_dof_per_element);
    
    // prepare the basis function on the quadrature points
    for (unsigned int q(0); q < n_q; ++q) {
      for (std::size_t i(0); i < n_dof; ++i)
	phi.at(q, i, 0) = fe_type<n>::phi(i, &xq_hat.at(q, 0));
    }
      
    // prepare the basis function derivatives on the quadrature points
    const array<double> jmt(m.get_jmt(k));
    for (unsigned int q(0); q < n_q; ++q) {
      for (std::size_t i(0); i < n_dof; ++i) {
	for (std::size_t s(0); s < dim; ++s) {
	  phi.at(q, i, 1 + s) = 0.0;
	  for (std::size_t t(0); t < dim; ++t)
	    phi.at(q, i, 1 + s) += jmt.at(s, t) * fe_type<n>::dphi(t, i, &xq_hat.at(q, 0));
	}
      }
    }
  }
};



#endif /* FE_VALUE_MANAGER_H */
