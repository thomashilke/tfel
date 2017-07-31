#ifndef SPARSE_LINEAR_SYSTEM_H
#define SPARSE_LINEAR_SYSTEM_H

#include <cstddef>
#include <utility>
#include <map>
#include <ostream>
#include <iomanip>

class sparse_linear_system {
public:
  sparse_linear_system(std::size_t n_eq, std::size_t n_unk)
    :n_equation(n_eq), n_unknown(n_unk) {}

  void set_sizes(std::size_t n_eq, std::size_t n_unk) {
    n_equation = n_eq;
    n_unknown = n_unk;
  }
  
  void accumulate(std::size_t i, std::size_t j, double v) {
    if (i >= n_equation or j >= n_unknown)
      throw std::string("sparse_linear_system::out of bound access");
      
    elements[std::pair<std::size_t, std::size_t>(i, j)] += v;
  }

  void clear() {
    elements.clear();
  }
  
  void show(std::ostream& stream) const {
    auto p(stream.precision(2));
    array<double> profile{n_equation, n_unknown};
    profile.fill(0.0);
    for (const auto& elem: elements)
      profile.at(elem.first.first, elem.first.second) = elem.second;

    stream << "a = [";
    for (unsigned int i(0); i < n_equation - 1; ++i) {
      for (unsigned int j(0); j < n_unknown - 1; ++j)
	stream << std::setw(4) << profile.at(i, j) << ", ";
      stream << std::setw(4) << profile.at(i, n_unknown - 1) << "; " << std::endl;
      
    }
    for (unsigned int j(0); j < n_unknown - 1; ++j)
      stream << std::setw(4) << profile.at(n_equation - 1, j) << ", ";
    stream << std::setw(4) << profile.at(n_equation - 1, n_unknown - 1) << "];" << std::endl;

    stream.precision(p);
  }

  std::size_t get_equation_number() const { return n_equation; }
  std::size_t get_equation_unknown() const { return n_unknown; }

  const std::map<std::pair<std::size_t, std::size_t>, double> get_elements() const { return elements; }

  std::size_t get_element_number() const { return elements.size(); }
  
private:
  std::size_t n_equation, n_unknown;
  std::map<std::pair<std::size_t, std::size_t>, double> elements;
};

#endif /* SPARSE_LINEAR_SYSTEM_H */
