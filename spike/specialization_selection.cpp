
#include <iostream>

template<typename T, typename V = void>
struct a {
  static constexpr int value = 1;
};

template<typename S>
struct a<S, void> {
  static constexpr int value = 2;
};

int main(int argc, char *argv[])
{
  std::cout << a<int>::value << std::endl;
  return 0;
}

