
#include <spikes/timer.hpp>

#include "../src/core/fes.hpp"
#include "../src/core/fe.hpp"
#include "../src/core/mesh.hpp"
#include "../src/core/cell.hpp"

void build_fes(std::size_t n) {
  const fe_mesh<cell::triangle> m(gen_square_mesh(1.0, 1.0, n, n));
  volatile finite_element_space<cell::triangle::fe::lagrange_p1> fes(m);
}

int main(int argc, char *argv[]) {
  std::vector<std::size_t> ns{10, 20, 30, 40, 50,
                              60, 70, 80, 90, 100,
      110, 120, 130, 140, 150,
      160, 170, 180, 190, 200,
      210, 220, 230, 240, 250,
      300, 350, 400, 450, 500,
      600, 700, 800, 1000};
  for (auto n: ns) {
    timer t;
    build_fes(n);
    std::cout << n << " " << t.tic() << std::endl;
  }
  
  return 0;
}
