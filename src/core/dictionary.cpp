#include "dictionary.hpp"

template<>
void dictionary::value<std::string>::print(std::ostream& stream) const {
  stream << "\"" << v << "\"";
}

template<>
void dictionary::value<bool>::print(std::ostream& stream) const {
  stream << std::boolalpha << v;
}
