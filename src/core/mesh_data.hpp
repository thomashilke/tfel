#ifndef _MESH_DATA_H_
#define _MESH_DATA_H_

#include <spikes/array.hpp>

#include "operator.hpp"
#include "meta.hpp"
#include "mesh.hpp"

template<typename value_t, typename mesh_t>
class mesh_data;


template<typename ast>
struct mesh_data_expression {
  using ast_type = ast;
  using value_type = typename ast_type::value_type;
  ast_type expression;

  std::size_t component_number() const { return expression.component_number(); }
  
  mesh_data_expression(ast_type expr): expression(expr) {}
  value_type operator()(std::size_t k, std::size_t n) const {
    return expression(k, n);
  }
};

template<typename value_t, typename mesh_t>
struct mesh_data_value {
  using value_type = value_t;
  using mesh_type = mesh_t;
  using cell_type = typename mesh_type::cell_type;

  std::size_t component_number() const { return data.component_number(); }
  
  mesh_data_value(const mesh_data<value_type, mesh_type>& d) : data(d) {}
  
  const mesh_data<value_type, cell_type>& data;
  value_type operator()(std::size_t k, std::size_t n) const { return data.value(k, n); }
};

template<typename value_t, typename mesh_t>
struct mesh_data_component_value {
  using value_type = value_t;
  using mesh_type = mesh_t;
  using cell_type = typename mesh_type::cell_type;

  std::size_t component_number() const { return 1ul; }

  mesh_data_component_value(const mesh_data<value_type, mesh_type>& d,
		      std::size_t n): data(d), component(n) {}
  
  const mesh_data<value_type, mesh_type>& data;
  const std::size_t component;
  value_type operator()(std::size_t k, std::size_t n) const { return data.value(k, component); }
};

template<typename value_t>
struct mesh_data_constant {
  using value_type = value_t;

  std::size_t component_number() const { return 1ul; }
  
  mesh_data_constant(double c): constant(c) {}

  const double constant;
  value_type operator()(std::size_t k, std::size_t n) const { return constant; }
};

template<typename left, typename right, typename op>
struct mesh_data_binary_op {
  using left_type = left;
  using right_type = right;
  using value_type = typename op::value_type;

  std::size_t component_number() const { return std::max(left_expr.component_number(),
							 right_expr.component_number()); }
  
  left_type left_expr;
  right_type right_expr;
  mesh_data_binary_op(const left_type& l, const right_type& r): left_expr(l), right_expr(r) {}

  value_type operator()(std::size_t k, std::size_t n) const { return op::apply(left_expr(k, n), right_expr(k, n)); }
};


#define DECL_NUMERIC_BINARY_OP(OP_SYMBOL, OP_IMPL) \
template<typename left_ast_type, typename right_ast_type> \
mesh_data_expression< \
  mesh_data_binary_op<left_ast_type, \
		      right_ast_type, \
		      OP_IMPL<typename left_ast_type::value_type> > > \
operator OP_SYMBOL(const mesh_data_expression<left_ast_type>& l, \
	   const mesh_data_expression<right_ast_type>& r) { \
  return mesh_data_expression< \
           mesh_data_binary_op<left_ast_type, \
			       right_ast_type, \
			       OP_IMPL<typename left_ast_type::value_type> > >( \
                                 mesh_data_binary_op<left_ast_type, \
			                             right_ast_type, \
						     OP_IMPL<typename left_ast_type::value_type> >( \
						       l.expression, r.expression)); \
}


#define DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(OP_SYMBOL, OP_IMPL) \
template<typename left_type, typename right_ast_type> \
mesh_data_expression< \
  mesh_data_binary_op<mesh_data_constant<left_type>, \
		      right_ast_type, \
		      OP_IMPL<typename right_ast_type::value_type> > > \
operator OP_SYMBOL( \
  const left_type& l,					    \
  const mesh_data_expression<right_ast_type>& r) { \
\
  return mesh_data_expression< \
           mesh_data_binary_op<mesh_data_constant<left_type>, \
			       right_ast_type, \
			       OP_IMPL<typename right_ast_type::value_type> > >( \
                                 mesh_data_binary_op<mesh_data_constant<left_type>, \
			                             right_ast_type, \
						     OP_IMPL<typename right_ast_type::value_type> >( \
				                       mesh_data_constant<left_type>(l), r.expression)); \
}

#define DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(OP_SYMBOL, OP_IMPL) \
template<typename left_ast_type, typename right_type> \
mesh_data_expression< \
  mesh_data_binary_op<left_ast_type, \
		      mesh_data_constant<right_type>,		       \
		      OP_IMPL<typename left_ast_type::value_type> > > \
operator OP_SYMBOL( \
  const mesh_data_expression<left_ast_type>& l,	\
  const right_type& r) { \
\
  return mesh_data_expression< \
           mesh_data_binary_op<left_ast_type, \
			       mesh_data_constant<right_type>, \
			       OP_IMPL<typename left_ast_type::value_type> > >( \
                                 mesh_data_binary_op<left_ast_type, \
			                             mesh_data_constant<right_type>,	\
				                     OP_IMPL<typename left_ast_type::value_type> >( \
                                                       l.expression, mesh_data_constant<right_type>(r)));	\
}

#define DECL_NUMERIC_BINARY_LOGICAL_OP(OP_SYMBOL, OP_IMPL) \
template<typename left_ast_type, typename right_ast_type> \
mesh_data_expression< \
  mesh_data_binary_op<left_ast_type, \
		      right_ast_type, \
		      OP_IMPL> > \
operator OP_SYMBOL(const mesh_data_expression<left_ast_type>& l, \
	   const mesh_data_expression<right_ast_type>& r) { \
  return mesh_data_expression< \
           mesh_data_binary_op<left_ast_type, \
			       right_ast_type, \
			       OP_IMPL> >(l.expression, r.expression); \
}

// homogeneous version
DECL_NUMERIC_BINARY_OP(*, multiply)
DECL_NUMERIC_BINARY_OP(/, divide)
DECL_NUMERIC_BINARY_OP(+, add)
DECL_NUMERIC_BINARY_OP(-, substract)

DECL_NUMERIC_BINARY_OP(>, greater_than)
DECL_NUMERIC_BINARY_OP(<, less_than)

DECL_NUMERIC_BINARY_OP(>=, greater_equal_than)
DECL_NUMERIC_BINARY_OP(<=, less_equal_than)

DECL_NUMERIC_BINARY_OP(==, equal)
DECL_NUMERIC_BINARY_OP(!=, not_equal)

// left constant version
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(*, multiply)
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(/, divide)
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(+, add)
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(-, substract)

DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(>, greater_than)
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(<, less_than)

DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(>=, greater_equal_than)
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(<=, less_equal_than)

DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(==, equal)
DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT(!=, not_equal)

// right constant version
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(*, multiply)
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(/, divide)
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(+, add)
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(-, substract)

DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(>, greater_than)
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(<, less_than)

DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(>=, greater_equal_than)
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(<=, less_equal_than)

DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(==, equal)
DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT(!=, not_equal)

// logical operators
DECL_NUMERIC_BINARY_LOGICAL_OP(&&, logical_and)
DECL_NUMERIC_BINARY_LOGICAL_OP(||, logical_or)

#undef DECL_NUMERIC_BINARY_OP
#undef DECL_NUMERIC_BINARY_OP_LEFT_CONSTANT
#undef DECL_NUMERIC_BINARY_OP_RIGHT_CONSTANT
#undef DECL_NUMERIC_BINARY_LOGICAL_OP


enum class mesh_data_kind {cell, vertex};

template<typename value_type, typename ... array_types>
struct concat_arrays;

template<typename value_type, typename array_type, typename ... array_types>
struct concat_arrays<value_type, array_type, array_types...> {
  static void call(array<value_type>& dst, const array_type& a, const array_types& ... as) {
    assert(a.get_rank() == 2 and a.get_size(1) == 1);
    assert(a.get_size(0) == dst.get_size(0));
    assert(dst.get_rank() == 2);
    assert(dst.get_size(1) >= sizeof...(array_types));
    
    const std::size_t m(dst.get_size(1) - (sizeof...(array_types) + 1));
    for (std::size_t n(0); n < a.get_size(0); ++n)
      dst.at(n, m) = a.at(n, 0);

    concat_arrays<value_type, array_types...>::call(dst, as...);
  }
};

template<typename value_type>
struct concat_arrays<value_type> {
  static void call(array<value_type>& dst) {}
};

template<typename value_t, typename mesh_t>
class mesh_data {
public:
  using mesh_type = mesh_t;
  using cell_type = typename mesh_type::cell_type;
  using value_type = value_t;

  
  mesh_data(const mesh_type& m, mesh_data_kind t, std::size_t n_component = 1ul)
    : m(m), type(t), values{value_number(t), n_component} {}

  mesh_data(const mesh_type& m, mesh_data_kind t, const array<double>& values)
    : m(m), type(t), values(values) {}

  mesh_data(const mesh_type& m, mesh_data_kind t, array<double>&& values)
    : m(m), type(t), values(std::move(values)) {}

  
  
  template<typename ast_type>
  mesh_data(const mesh_type& m, mesh_data_kind t, const mesh_data_expression<ast_type>& expr)
    : m(m), type(t), values{value_number(t), expr.component_number()} {
      
      for (std::size_t k(0); k < values.get_size(0); ++k)
	for (std::size_t n(0); n < values.get_size(1); ++n)
	  values.at(k, n) = expr(k, n);
  }

  template<typename ...  mesh_data_types>
  mesh_data(const mesh_type& m, mesh_data_kind t, const mesh_data_types& ... datas)
    : m(m), type(t), values{value_number(t), sizeof...(mesh_data_types)}{
    concat_arrays<double, array<typename mesh_data_types::value_type>...>::call(values, (datas.get_values())...);
  }
  
  const value_type& value(std::size_t k, std::size_t n) const { return values.at(k, n); }
  value_type& value(std::size_t k, std::size_t n) { return values.at(k, n); }

  
  mesh_data_expression<mesh_data_component_value<value_type, mesh_type > >
  operator[](std::size_t n) const {
    return mesh_data_expression<mesh_data_component_value<value_type, mesh_type > >(
      mesh_data_component_value<value_type, mesh_type>(*this, n));
  }
  
  template<typename ast_type>
  mesh_data& operator=(const mesh_data_expression<ast_type>& expr) {
    if (expr.component_number() != values.get_size(1))
      values = array<double>{value_number(type), expr.component_number()};
	
    for (std::size_t k(0); k < values.get_size(0); ++k)
      for (std::size_t n(0); n < values.get_size(1); ++n)
	values.at(k, n) = expr(k, n);
    
    return *this;
  }

  const array<value_type>& get_values() const { return values; }

  std::size_t get_component_number() const { return values.get_size(1); }

  mesh_data_kind get_kind() const { return type; }

  const mesh_type& get_mesh() const { return m; }
  
private:
  const mesh_type& m;
  mesh_data_kind type;
  array<value_type> values;

  std::size_t value_number(mesh_data_kind type) const {
    switch (type) {
    case mesh_data_kind::cell: return m.get_element_number();
    case mesh_data_kind::vertex: return m.get_vertex_number();
    }
  }
};

template<typename cell_type>
mesh_data<double, submesh<cell_type, typename cell_type::boundary_cell_type> >
compute_boundary_normals(const submesh<cell_type>& dm) {
  mesh_data<double, submesh<cell_type> >
    result(dm, mesh_data_kind::cell,
	   dm.get_embedding_space_dimension());
  
  for (std::size_t k(0); k < dm.get_element_number(); ++k) {
    const auto normal(cell_type::subdomain_normal(dm.get_vertices(),
						  dm.get_mesh().get_elements(),
						  dm.get_parent_element_id(k),
						  dm.get_subdomain_id(k)));

    for (std::size_t n(0); n < dm.get_embedding_space_dimension(); ++n)
      result.value(k, n) = normal.at(n);
  }

  return result;
}


template<typename mesh_type, typename ... Fs>
mesh_data<double, mesh_type>
evaluate_on_cells(const mesh_type& m, Fs... fs) {
  using cell_type = typename mesh_type::cell_type;

  mesh_data<double, mesh_type> result(m, mesh_data_kind::cell,
				      sizeof...(Fs));

  for (std::size_t k(0); k < m.get_element_number(); ++k) {
    const auto x(cell_type::barycenter(m.get_vertices(), m.get_elements(), k));
    double values[sizeof...(Fs)];
    array_from<double, return_2nd_t<Fs, double>...>::parameters(values, fs(&x.at(0))...);
    for (std::size_t n(0); n < sizeof...(Fs); ++n)
      result.value(k, n) = values[n];
  }

  return result;
}

template<typename mesh_type, typename ... Fs>
mesh_data<double, mesh_type>
evaluate_on_vertices(const mesh_type& m, Fs... fs) {
  mesh_data<double, mesh_type> result(m, mesh_data_kind::vertex,
				      sizeof...(Fs));

  for (std::size_t k(0); k < m.get_vertex_number(); ++k) {
    double values[sizeof...(Fs)];
    array_from<double, return_2nd_t<Fs, double>...>::parameters(values, fs(&m.get_vertices().at(k, 0))...);
    for (std::size_t n(0); n < sizeof...(Fs); ++n)
      result.value(k, n) = values[n];
  }

  return result;
}
    

#endif /* _MESH_DATA_H_ */
