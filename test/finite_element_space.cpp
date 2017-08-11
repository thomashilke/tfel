
#include <iostream>

#include <spikes/array.hpp>

#include "../src/core/cell.hpp"
#include "../src/core/fe.hpp"
#include "../src/core/mesh.hpp"
#include "../src/core/fes.hpp"

void test_triangle_lagrange_p1() {
  double vertices[] = { 0.0, 0.0,
		        1.0, 0.0,
		       -0.5, 1.0,
		        0.5, 1.0};
  unsigned int elements[] = { 0, 1, 3,
			      0, 2, 3};

  using cell_type = cell::triangle;
  fe_mesh<cell_type> m(vertices, 4, 2,
			 elements, 2);

  using fe_type = cell_type::fe::lagrange_p1;
  finite_element_space<fe_type> fes(m);

  fes.show(std::cout);
}

void test_triangle_lagrange_p1_bubble() {
  double vertices[] = { 0.0, 0.0,
		        1.0, 0.0,
		       -0.5, 1.0,
		        0.5, 1.0};
  unsigned int elements[] = { 0, 1, 3,
			      0, 2, 3};

  using cell_type = cell::triangle;
  fe_mesh<cell_type> m(vertices, 4, 2,
			 elements, 2);

  using fe_type = cell_type::fe::lagrange_p1_bubble;
  finite_element_space<fe_type> fes(m);

  fes.show(std::cout);
}

void test_edge_lagrange_p1() {
  double vertices[] = { 0.0, 1.0, 2.0};
  unsigned int elements[] = { 0, 1,
			      1, 2};

  using cell_type = cell::edge;
  fe_mesh<cell_type> m(vertices, 3, 1,
		     elements, 2);

  using fe_type = cell_type::fe::lagrange_p1;
  finite_element_space<fe_type> fes(m);

  fes.show(std::cout);
}

void test_edge_lagrange_p1_bubble() {
  double vertices[] = { 0.0, 1.0, 2.0};
  unsigned int elements[] = { 0, 1,
			      1, 2};

  using cell_type = cell::edge;
  fe_mesh<cell_type> m(vertices, 3, 1,
		     elements, 2);

  using fe_type = cell_type::fe::lagrange_p1_bubble;
  finite_element_space<fe_type> fes(m);

  fes.show(std::cout);
}

int main(int argc, char *argv[]) {
  test_triangle_lagrange_p1();
  test_triangle_lagrange_p1_bubble();
  test_edge_lagrange_p1();
  test_edge_lagrange_p1_bubble();

  return 0;
}
