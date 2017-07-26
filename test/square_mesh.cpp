
#include <iostream>

#include "../src/core/mesh.hpp"

int main(int argc, char *argv[]) {
  mesh<cell::triangle> m(gen_square_mesh(1.0, 1.0, 2, 3));
  m.show(std::cout);
  
  return 0;
}
