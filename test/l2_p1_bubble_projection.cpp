
#include "../src/core/projector.hpp"
#include "../src/core/export.hpp"

double f(const double* x) {
  return std::sin(x[0] * M_PI) + std::sin(x[1] * M_PI);
  //return - (x[0] - x[1] + 1.0) * (x[0] - x[1] - 1.0);
}


int main(int argc, char *argv[]) {
  //using fe_type = cell::triangle::fe::lagrange_p1;
  //using fe_type = cell::triangle::fe::lagrange_p1_bubble;
  using fe_type = cell::triangle::fe::lagrange_p2;

  fe_mesh<cell::triangle> m(gen_square_mesh(1.0, 1.0, 10, 10));
  finite_element_space<fe_type> fes(m);

  const auto f_h(projector::l2<fe_type, quad::triangle::qf5pT>(f, fes));
  
  exporter::ensight6("l2_proj_p1_bubble",
                     to_mesh_vertex_data<fe_type>(f_h), "f_h");

  std::cout << integrate<quad::triangle::qf5pT>((make_expr<fe_type>(f_h) - make_expr(f))*
                                                (make_expr<fe_type>(f_h) - make_expr(f)), m)
            << std::endl;
  
  return 0;
}
