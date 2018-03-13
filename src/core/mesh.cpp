#include "mesh.hpp"

fe_mesh<cell::edge> gen_segment_mesh(double x_1, double x_2, unsigned int n) {
  std::vector<double> vertices(n + 1);
  std::vector<unsigned int> elements(2 * n);

  for (unsigned int i(0); i < n + 1; ++i)
    vertices[i] = x_1 + (x_2 - x_1) / n * i;

  for (unsigned int k(0); k < n; ++k) {
    elements[2 * k] = k;
    elements[2 * k + 1] = k + 1;
  }

  return fe_mesh<cell::edge>(&vertices[0], n + 1, 1,
			  &elements[0], n);
}

fe_mesh<cell::triangle> gen_square_mesh(double x_1, double x_2,
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
      // first diagonal (south-west north-east)
      //elements[6 * (n_1 * j + i) + 0] = (j    ) * (n_1 + 1) + (i    );
      //elements[6 * (n_1 * j + i) + 1] = (j    ) * (n_1 + 1) + (i + 1);
      //elements[6 * (n_1 * j + i) + 2] = (j + 1) * (n_1 + 1) + (i    );
      //elements[6 * (n_1 * j + i) + 3] = (j    ) * (n_1 + 1) + (i + 1);
      //elements[6 * (n_1 * j + i) + 4] = (j + 1) * (n_1 + 1) + (i    );
      //elements[6 * (n_1 * j + i) + 5] = (j + 1) * (n_1 + 1) + (i + 1);

      // second diagonal (south-east north-west)
      elements[6 * (n_1 * j + i) + 0] = (j    ) * (n_1 + 1) + (i    );
      elements[6 * (n_1 * j + i) + 1] = (j    ) * (n_1 + 1) + (i + 1);
      elements[6 * (n_1 * j + i) + 2] = (j + 1) * (n_1 + 1) + (i + 1);
      elements[6 * (n_1 * j + i) + 3] = (j    ) * (n_1 + 1) + (i    );
      elements[6 * (n_1 * j + i) + 4] = (j + 1) * (n_1 + 1) + (i    );
      elements[6 * (n_1 * j + i) + 5] = (j + 1) * (n_1 + 1) + (i + 1);
    }

  return fe_mesh<cell::triangle>(&vertices[0], (n_1 + 1) * (n_2 + 1), 2,
                                 &elements[0], n_1 * n_2 * 2);
}


fe_mesh<cell::tetrahedron> gen_cube_mesh(double x_1, double x_2, double x_3,
                                         unsigned int n_1, unsigned int n_2, unsigned int n_3) {
  array<double> vertices{n_3 + 1, n_2 + 1, n_1 + 1, 3};
  array<unsigned int> v_ids{n_3 + 1, n_2 + 1, n_1 + 1};
  array<unsigned int> elements{n_3, n_2, n_1, 5, 4};

  const double
    h_x(x_1 / static_cast<double>(n_1)),
    h_y(x_2 / static_cast<double>(n_2)),
    h_z(x_3 / static_cast<double>(n_3));

  unsigned int v_id(0);
  for (unsigned int k(0); k <= n_3; ++k)
    for (unsigned int j(0); j <= n_2; ++j)
      for (unsigned int i(0); i <= n_1; ++i) {
	v_ids.at(k, j, i) = v_id;
	v_id += 1;

	vertices.at(k, j, i, 0) = i * h_x;
	vertices.at(k, j, i, 1) = j * h_y;
	vertices.at(k, j, i, 2) = k * h_z;
      }

  for (unsigned int k(0); k < n_3; ++k)
    for (unsigned int j(0); j < n_2; ++j)
      for (unsigned int i(0); i < n_1; ++i) {
	if ((i + j + k) % 2 == 0) {
	  elements.at(k, j, i, 0, 0) = v_ids.at(k,     j,     i);
	  elements.at(k, j, i, 0, 1) = v_ids.at(k,     j,     i + 1);
	  elements.at(k, j, i, 0, 2) = v_ids.at(k,     j + 1, i);
	  elements.at(k, j, i, 0, 3) = v_ids.at(k + 1, j,     i);

	  elements.at(k, j, i, 1, 0) = v_ids.at(k,     j,     i + 1);
	  elements.at(k, j, i, 1, 1) = v_ids.at(k + 1, j,     i);
	  elements.at(k, j, i, 1, 2) = v_ids.at(k + 1, j,     i + 1);
	  elements.at(k, j, i, 1, 3) = v_ids.at(k + 1, j + 1, i + 1);

	  elements.at(k, j, i, 2, 0) = v_ids.at(k, j, i + 1);
	  elements.at(k, j, i, 2, 1) = v_ids.at(k, j + 1, i);
	  elements.at(k, j, i, 2, 2) = v_ids.at(k, j + 1, i + 1);
	  elements.at(k, j, i, 2, 3) = v_ids.at(k + 1, j + 1, i + 1);

	  elements.at(k, j, i, 3, 0) = v_ids.at(k,     j + 1, i);
	  elements.at(k, j, i, 3, 1) = v_ids.at(k + 1, j,     i);
	  elements.at(k, j, i, 3, 2) = v_ids.at(k + 1, j + 1, i);
	  elements.at(k, j, i, 3, 3) = v_ids.at(k + 1, j + 1, i + 1);

	  elements.at(k, j, i, 4, 0) = v_ids.at(k,     j,     i + 1);
	  elements.at(k, j, i, 4, 1) = v_ids.at(k,     j + 1, i);
	  elements.at(k, j, i, 4, 2) = v_ids.at(k + 1, j,     i);
	  elements.at(k, j, i, 4, 3) = v_ids.at(k + 1, j + 1, i + 1);
	} else {
	  elements.at(k, j, i, 0, 0) = v_ids.at(k,     j,     i);
	  elements.at(k, j, i, 0, 1) = v_ids.at(k,     j,     i + 1);
	  elements.at(k, j, i, 0, 2) = v_ids.at(k,     j + 1, i + 1);
	  elements.at(k, j, i, 0, 3) = v_ids.at(k + 1, j,     i + 1);

	  elements.at(k, j, i, 1, 0) = v_ids.at(k,     j,     i);
	  elements.at(k, j, i, 1, 1) = v_ids.at(k    , j + 1, i);
	  elements.at(k, j, i, 1, 2) = v_ids.at(k    , j + 1, i + 1);
	  elements.at(k, j, i, 1, 3) = v_ids.at(k + 1, j + 1, i);

	  elements.at(k, j, i, 2, 0) = v_ids.at(k,     j,     i);
	  elements.at(k, j, i, 2, 1) = v_ids.at(k + 1, j,     i);
	  elements.at(k, j, i, 2, 2) = v_ids.at(k + 1, j,     i + 1);
	  elements.at(k, j, i, 2, 3) = v_ids.at(k + 1, j + 1, i);

	  elements.at(k, j, i, 3, 0) = v_ids.at(k,     j + 1, i + 1);
	  elements.at(k, j, i, 3, 1) = v_ids.at(k + 1, j,     i + 1);
	  elements.at(k, j, i, 3, 2) = v_ids.at(k + 1, j + 1, i);
	  elements.at(k, j, i, 3, 3) = v_ids.at(k + 1, j + 1, i + 1);

	  elements.at(k, j, i, 4, 0) = v_ids.at(k,     j,     i);
	  elements.at(k, j, i, 4, 1) = v_ids.at(k,     j + 1, i + 1);
	  elements.at(k, j, i, 4, 2) = v_ids.at(k + 1, j,     i + 1);
	  elements.at(k, j, i, 4, 3) = v_ids.at(k + 1, j + 1, i);
	}
      }
  
  return fe_mesh<cell::tetrahedron>(&vertices.at(0,0,0,0), (n_1 + 1) * (n_2 + 1) * (n_3 + 1), 3,
                                    &elements.at(0,0,0,0,0), n_1 * n_2 * n_3 * 5);
}
