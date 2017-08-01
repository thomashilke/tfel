
#include <type_traits>
#include <cstddef>
#include <iostream>
#include <typeinfo>


#include "../src/core/mesh.hpp"
#include "../src/core/fe.hpp"
#include "../src/core/fes.hpp"
#include "../src/core/meta.hpp"
#include "../src/core/fe_value_manager.hpp"
#include "../src/core/form.hpp"
#include "../src/core/quadrature.hpp"
#include "../src/core/export.hpp"
#include "../src/core/fes.hpp"


#include "../src/core/composite_fe.hpp"
#include "../src/core/composite_fes.hpp"
#include "../src/core/composite_form.hpp"


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

    using i_tl = make_index_list_t<tl>;
    print_type<i_tl>("index list generated from tl: ");
    
    std::cout << is_member<int, tl>::value << std::endl;
    std::cout << is_member<double, tl>::value << std::endl;
    std::cout << is_member<char, tl>::value << std::endl;
    std::cout << is_member<short, tl>::value << std::endl;

    std::cout << "index of double is: " << get_index_of_element<double, tl>::value << std::endl;
    std::cout << "index of char is: " << get_index_of_element<char, tl>::value << std::endl;
    std::cout << "index of int is: " << get_index_of_element<int, tl>::value << std::endl;
  }

  {
    using tl = type_list<double, char, double>;

    print_type<unique_t<tl> >();
  }
}


void test_meta_integral_sequence() {
  using is = make_integral_sequence_t<std::size_t, 2>;
  using il = sequence_to_list_t<is>;
  using ill = wrap_t<type_list, il>;
  using ic = integral_constant<int, 7>;
  using ic_ill = append_to_each_element_t<ic, ill>;
  using tp = tensor_product_of_lists_t<il, il>;
  
  print_type<il>("list of integral_constant: ");

  print_type<wrap_t<type_list, il> >("wrap il elements in type_list's: ");
  print_type<ic_ill>("ic appended: ");

  print_type<tp>("tensor_product of the list with itself: ");
  print_type<append_to_each_element_t<int, tp> >("tensor product with int appended: ");
}


double f_0(const double* x) { return 1.0; }
double f_1(const double* x) { return -1.0; }

void test_cfe() {
  using cell_type = cell::triangle;
  using fe_0_type = finite_element::triangle_lagrange_p1;
  using fe_1_type = finite_element::triangle_lagrange_p1;
  using fe_type = composite_finite_element<fe_0_type, fe_1_type>;
  using fes_type = composite_finite_element_space<fe_type>;
  
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 100, 100));
  submesh<cell_type> dm(m.get_boundary_submesh());
  fes_type fes(m);

  fes.set_dirichlet_boundary<0>(dm);
  fes.set_dirichlet_boundary<1>(dm);
  
  //fes.show<0>(std::cout);
  //fes.show<1>(std::cout);


  bilinear_form<fes_type, fes_type> a(fes, fes); {
    const auto u_0(a.get_trial_function<0>());
    const auto u_1(a.get_trial_function<1>());
    const auto v_0(a.get_test_function<0>());
    const auto v_1(a.get_test_function<1>());

    a += integrate<quad::triangle::qf5pT>(d<1>(u_0) * d<1>(v_0) + d<2>(u_0) * d<2>(v_0)
					  + d<1>(u_1) * d<1>(v_1) + d<2>(u_1) * d<2>(v_1), m);
  }

  linear_form<fes_type> f(fes); {
    const auto v_0(f.get_test_function<0>());
    const auto v_1(f.get_test_function<1>());

    f += integrate<quad::triangle::qf5pT>(f_0 * v_0 + f_1 * v_1, m);
  }
  //f.show(std::cout);

  const fes_type::element x(a.solve(f));

  exporter::ensight6("laplacien_0", x.get_component<0>(), "solution");
  exporter::ensight6("laplacien_1", x.get_component<1>(), "solution");
}


void test_fe_value_manager() {
  using fe_list = type_list<finite_element::triangle_lagrange_p0,
			    finite_element::triangle_lagrange_p1>;

  const std::size_t n_q(3);
  array<double> xq{3, 2};
  xq.at(0, 0) = 0.0;
  xq.at(0, 0) = 0.0;
    
  xq.at(0, 0) = 1.0;
  xq.at(0, 0) = 0.0;
	
  xq.at(0, 0) = 0.0;
  xq.at(0, 0) = 1.0;

  array<double> jmt{2, 2};
  jmt.at(0, 0) = 1.0;
  jmt.at(0, 1) = 0.0;
  jmt.at(1, 0) = 0.0;
  jmt.at(1, 1) = 1.0;
  
  fe_value_manager<fe_list> fe_values(n_q);
  fe_values.set_points(xq);
  fe_values.prepare(jmt);

  std::cout << fe_values.get_values<0ul>().at(0,0,0) << std::endl;
  const array<double>& phi(fe_values.get_values<0>());
  const array<double>& psi(fe_values.get_values<1>());
  
  std::cout << phi.at(0,0,0) << std::endl;
  std::cout << psi.at(0,0,0) << std::endl;
}


void test_cfes() {
  
}


int main(int argc, char *argv[]) {
  try {
  //test_meta();
  test_cfe();
  //test_fe_value_manager();
  //test_cfes();
  //test_meta_integral_sequence();
  } catch (std::string& e) {
    std::cout << e << std::endl;
  }

  return 0;
}
