
#include <iostream>

#include <spikes/array.hpp>

#include "cell.hpp"
#include "fe.hpp"


int main(int argc, char *argv[]) {
  
  std::cout << "hello, world!" << std::endl;

  double x(0.5);
  unsigned int d(1);

  std::cout << finite_element::edge_lagrange_p1::basis_function(1, &d, &x) << std::endl;
  
  return 0;
}
