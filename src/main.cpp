
#include <iostream>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"
#include "mesh.hpp"
#include "fes.hpp"

void test_triangle_lagrange_p1() {
  double vertices[] = { 0.0, 0.0,
		        1.0, 0.0,
		       -0.5, 1.0,
		        0.5, 1.0};
  unsigned int elements[] = { 0, 1, 3,
			      0, 2, 3};

  mesh<cell::triangle> m(vertices, 4, 2,
			 elements, 2);

  typedef finite_element::triangle_lagrange_p1 fe;
  finite_element_space<fe> fes(m);

  fes.show(std::cout);
}

void test_triangle_lagrange_p1_bubble() {
  double vertices[] = { 0.0, 0.0,
		        1.0, 0.0,
		       -0.5, 1.0,
		        0.5, 1.0};
  unsigned int elements[] = { 0, 1, 3,
			      0, 2, 3};

  mesh<cell::triangle> m(vertices, 4, 2,
			 elements, 2);

  typedef finite_element::triangle_lagrange_p1_bubble fe;
  finite_element_space<fe> fes(m);

  fes.show(std::cout);
}

void test_edge_lagrange_p1() {
  double vertices[] = { 0.0, 1.0, 2.0};
  unsigned int elements[] = { 0, 1,
			      1, 2};

  mesh<cell::edge> m(vertices, 3, 1,
		     elements, 2);

  typedef finite_element::edge_lagrange_p1 fe;
  finite_element_space<fe> fes(m);

  fes.show(std::cout);
}

void test_edge_lagrange_p1_bubble() {
  double vertices[] = { 0.0, 1.0, 2.0};
  unsigned int elements[] = { 0, 1,
			      1, 2};

  mesh<cell::edge> m(vertices, 3, 1,
		     elements, 2);

  typedef finite_element::edge_lagrange_p1_bubble fe;
  finite_element_space<fe> fes(m);

  fes.show(std::cout);
}

int main(int argc, char *argv[]) {
  test_triangle_lagrange_p1();
  test_triangle_lagrange_p1_bubble();
  test_edge_lagrange_p1();
  test_edge_lagrange_p1_bubble();

  return 0;
}
