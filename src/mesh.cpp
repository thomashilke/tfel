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

mesh<cell::triangle> gen_square_mesh(double x_1, double x_2,
				     unsigned int n_1, unsigned int n_2) {
  std::vector<double> vertices((n_1 + 1) * (n_2 + 1) * 2);
  std::vector<unsigned int> elements(2 * n_1 * n_2 * 3);

  for (unsigned int j(0); j <= n_2; ++j)
    for (unsigned int i(0); i <= n_1; ++i) {
      vertices[2 * (j * (n_1 + 1) + i)] = x_1 / n_1 * i;
      vertices[2 * (j * (n_1 + 1) + i) + 1] = x_2 / n_2 * j;
    }

  for (unsigned int j(0); j < n_2; ++j)
    for (unsigned int i(0); i < n_1; ++i) {
      elements[6 * (n_1 * j + i) + 0] = (j    ) * (n_1 + 1) + (i    );
      elements[6 * (n_1 * j + i) + 1] = (j    ) * (n_1 + 1) + (i + 1);
      elements[6 * (n_1 * j + i) + 2] = (j + 1) * (n_1 + 1) + (i    );
      elements[6 * (n_1 * j + i) + 3] = (j    ) * (n_1 + 1) + (i + 1);
      elements[6 * (n_1 * j + i) + 4] = (j + 1) * (n_1 + 1) + (i    );
      elements[6 * (n_1 * j + i) + 5] = (j + 1) * (n_1 + 1) + (i + 1);
    }

  return mesh<cell::triangle>(&vertices[0], (n_1 + 1) * (n_2 + 1), 2,
			      &elements[0], n_1 * n_2 * 2);
}
