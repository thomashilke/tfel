#include "cell.hpp"

const std::size_t cell::point::n_subdomain_of_type[1] = {1};
const std::size_t cell::edge::n_subdomain_of_type[2] = {2, 1};
const std::size_t cell::triangle::n_subdomain_of_type[3] = {3, 3, 1};
