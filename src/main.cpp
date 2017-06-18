
#include <iostream>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"
#include "mesh.hpp"
#include "fes.hpp"

int main(int argc, char *argv[]) {

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
  
  return 0;
}
