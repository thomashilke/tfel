#ifndef BASIC_SIMPLEX_H
#define BASIC_SIMPLEX_H

namespace cell {

  template<std::size_t n>
  struct basic_simplex {
    static const std::size_t n_dimension = n;
    static const std::size_t n_vertex_per_element = n + 1;
    static const bool is_simplicial = true;

    static array<double> map_points_to_space_coordinates(const array<double>& vertices,
							 const array<unsigned int>& elements,
							 std::size_t k,
							 const array<double>& xs) {
      array<double> hat_xs{xs.get_size(0), vertices.get_size(1)};

      for (std::size_t i(0); i < xs.get_size(0); ++i) {
	for (std::size_t n(0); n < xs.get_size(1); ++n) {
	  hat_xs.at(i, n) = vertices.at(elements.at(k, 0), n);
	  for (std::size_t m(0); m < n_vertex_per_element - 1; ++m)
	    hat_xs.at(i, n) += xs.at(i, m) * (vertices.at(elements.at(k, m + 1), n)
					      - vertices.at(elements.at(k, 0), n));
	}
      }
      
      return hat_xs;
    }

    static array<double> get_jmt(const array<double>& vertices,
				 const array<unsigned int>& elements,
				 unsigned int k) {
      assert(vertices.get_size(1) == n_dimension);
      
      array<double> jmt{n_dimension, n_dimension};

      for (std::size_t i(0); i < n_dimension; ++i)
	for (std::size_t j(0); j < n_dimension; ++j) {
	  jmt.at(j, i) = vertices.at(elements.at(k, 1 + i), j) - vertices.at(elements.at(k, 0), j);
	}

      matrix_inverse(&jmt.at(0,0), n_dimension);

      for (std::size_t i(0); i < n_dimension - 1; ++i)
	for (std::size_t j(i + 1); j < n_dimension; ++j)
	  std::swap(jmt.at(i, j), jmt.at(j, i));

      return jmt;
    }

    static double element_diameter(const array<double>& vertices,
				   const array<unsigned int>& elements,
				   unsigned int k) {
      double longest_side(0.0);
      for (unsigned int i(0); i < n_vertex_per_element - 1; ++i)
	for (unsigned int j(i + 1); j < n_vertex_per_element; ++j) {
	  double side_length(0.0);
	  for (unsigned int n(0); n < vertices.get_size(1); ++n)
	    side_length += std::pow(vertices.at(elements.at(k, i), n) -
				    vertices.at(elements.at(k, j), n),
				    2);
	  side_length = std::sqrt(side_length);

	  longest_side = std::max(longest_side, side_length);
	}
      return longest_side;
    }

    static array<double> get_barycentric_coordinate_map(const array<double>& vertices,
							const array<unsigned int>& elements,
							std::size_t k) {
      array<double> map{n_vertex_per_element, n_vertex_per_element};

      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	map.at(0, j) = 1.0;
      
      for (std::size_t j(0); j < n_vertex_per_element; ++j)
	for (std::size_t i(0); i < n_dimension; ++i)
	  map.at(i + 1, j) = vertices.at(elements.at(k, j), i);

      matrix_inverse(&map.at(0, 0), n_vertex_per_element);
      
      return map;
    }

    static array<double> get_barycentric_coordinates(const array<double>& vertices,
						     const array<unsigned int>& elements,
						     std::size_t k,
						     const double* x) {
      const array<double> bc_map(get_barycentric_coordinate_map(vertices, elements, k));
      array<double>
	x_coord{n_vertex_per_element},
	bc_coord{n_vertex_per_element};

	x_coord.at(0) = 1.0;
	for (std::size_t i(0); i < n_dimension; ++i)
	  x_coord.at(i + 1) = x[i];

	bc_coord.fill(0.0);
	for (std::size_t i(0); i < n_vertex_per_element; ++i)
	  for (std::size_t j(0); j < n_vertex_per_element; ++j)
	    bc_coord.at(i) += bc_map.at(i, j) * x_coord.at(j);

	return bc_coord;
    }
  }

}

#endif /* BASIC_SIMPLEX_H */
