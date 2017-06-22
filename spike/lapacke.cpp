
#include <string>
#include <algorithm>
#include <iostream>

#include "../src/linear_algebra.hpp"

int main(int argc, char *argv[]) {
  double a[2][2] = {{0.90947, 0.55757},
		    {0.22113,  0.17065}};
  double b[2][2] = {};

  matrix_inverse(&b[0][0], &a[0][0], 2);

  std::cout.precision(16);
  
  std::cout << b[0][0] << " " << b[0][1] << std::endl;
  std::cout << b[1][0] << " " << b[1][1] << std::endl;

  std::cout << std::endl
	    << "Product of a and b should be the identity" << std::endl;
  std::cout << a[0][0] * b[0][0] + a[0][1] * b[1][0] << " ";
  std::cout << a[0][0] * b[0][1] + a[0][1] * b[1][1] << std::endl;
  std::cout << a[1][0] * b[0][0] + a[1][1] * b[1][0] << " ";
  std::cout << a[1][0] * b[0][1] + a[1][1] * b[1][1] << std::endl;
}
