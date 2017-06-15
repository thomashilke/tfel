
#include <iostream>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"
#include "mesh.hpp"

int main(int argc, char *argv[]) {

  double vertices[] = { 0.0, 0.0,
		        1.0, 0.0,
		       -0.5, 1.0,
		        0.5, 1.0};
  unsigned int elements[] = { 0, 1, 3,
			      0, 2, 3};

  mesh<cell::triangle> m(vertices, 4, 2,
			 elements, 2);
  
  return 0;
}
