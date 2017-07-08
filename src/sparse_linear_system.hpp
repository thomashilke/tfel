#ifndef SPARSE_LINEAR_SYSTEM_H
#define SPARSE_LINEAR_SYSTEM_H

#include <cstddef>
#include <utility>
#include <map>
#include <ostream>

class sparse_linear_system {
public:
  sparse_linear_system(std::size_t n_eq, std::size_t n_unk)
    :n_equation(n_eq), n_unknown(n_unk) {}
  void accumulate(std::size_t i, std::size_t j, double v) {
    elements[std::pair<std::size_t, std::size_t>(i, j)] += v;
  }

  void show(std::ostream& stream) const {
    array<double> profile{n_equation, n_unknown};
    for (const auto& elem: elements)
      profile.at(elem.first.first, elem.first.second) = elem.second;

    stream << "a = [";
    for (unsigned int i(0); i < n_equation - 1; ++i) {
      for (unsigned int j(0); j < n_unknown - 1; ++j)
	stream << std::setw(4) << std::setprecision(12) << profile.at(i, j) << ", ";
      stream << std::setw(4) << std::setprecision(12) << profile.at(i, n_unknown - 1) << "; " << std::endl;
      
    }
    for (unsigned int j(0); j < n_unknown - 1; ++j)
      stream << std::setw(4) << std::setprecision(12) << profile.at(n_equation - 1, j) << ", ";
    stream << std::setw(4) << std::setprecision(12) << profile.at(n_equation - 1, n_unknown - 1) << "];" << std::endl;
  }
  
private:
  std::size_t n_equation, n_unknown;
  std::map<std::pair<std::size_t, std::size_t>, double> elements;
};

#endif /* SPARSE_LINEAR_SYSTEM_H */
