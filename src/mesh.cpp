#include "mesh.hpp"

mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n) {
  std::vector<double> vertices(n + 1);
  std::vector<unsigned int> elements(2 * n);

  for (unsigned int i(0); i < n + 1; ++i)
    vertices[i] = x_1 + (x_2 - x_1) / n * i;

  for (unsigned int k(0); k < n; ++k) {
    elements[2 * k] = k;
    elements[2 * k + 1] = k + 1;
  }

  return mesh<cell::edge>(&vertices[0], n + 1, 1,
			  &elements[0], n);
}
