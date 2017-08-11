#include "../src/core/mesh_data.hpp"
#include "../src/core/export.hpp"
#include "../src/core/fe.hpp"
#include "../src/core/fes.hpp"
#include "../src/core/composite_fe.hpp"
#include "../src/core/composite_fes.hpp"
#include "../src/core/projector.hpp"

double b0(const double* x) {
  return -x[1];
}

double b1(const double* x) {
  return x[0];
}

double b2(const double* x) {
  return 0.0;
}


void test_1() {
  using cell_type = cell::triangle;
  using mesh_type = fe_mesh<cell_type>;
  
  fe_mesh<cell_type> m(gen_square_mesh(1.0, 1.0, 10, 10));

  mesh_data<double, mesh_type> data(evaluate_on_vertices(m, b0, b1, b2));

  exporter::ensight6("multicomponent_export", data, "data");
}

void test_2() {
  using cell_type = cell::triangle;
  using mesh_type = fe_mesh<cell_type>;
  
  mesh_type m(gen_square_mesh(1.0, 1.0, 10, 10));
  using fe_type = cell_type::fe::lagrange_p1;
  using cfe_type = composite_finite_element<fe_type, fe_type, fe_type>;
  using cfes_type = composite_finite_element_space<cfe_type>;
  using element_type = cfes_type::element;
  
  cfes_type cfes(m);
  element_type v(cfes,
		 projector::lagrange<fe_type>(b0, cfes.get_finite_element_space<0>()),
		 projector::lagrange<fe_type>(b1, cfes.get_finite_element_space<1>()),
		 projector::lagrange<fe_type>(b2, cfes.get_finite_element_space<2>()));

  
  mesh_data<double, mesh_type> data(to_mesh_vertex_data<cfe_type>(v));

  exporter::ensight6("export_fes_element", data, "data");
}

void test_3() {
  using cell_type = cell::triangle;
  using mesh_type = fe_mesh<cell_type>;
  
  mesh_type m(gen_square_mesh(1.0, 1.0, 10, 10));
  using fe_type = cell_type::fe::lagrange_p1;
  using cfe_type = composite_finite_element<fe_type, fe_type, fe_type>;
  using cfes_type = composite_finite_element_space<cfe_type>;
  using element_type = cfes_type::element;

  mesh_data<double, mesh_type> data_a(evaluate_on_vertices(m, b0, b1, b2));
  
  cfes_type cfes(m);
  element_type v(to_composite_p1_finite_element_function<3>(data_a, cfes));

  
  mesh_data<double, mesh_type> data(to_mesh_vertex_data<cfe_type>(v));

  exporter::ensight6("export_mesh_data_to_fes_element", data, "data");
}

int main(int argc, char *argv[]) {

  test_1();
  test_2();
  test_3();
  
  return 0;
}

