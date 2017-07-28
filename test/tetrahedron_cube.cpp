
#include "../src/core/mesh.hpp"
#include "../src/core/export.hpp"


int main(int argc, char *argv[]) {
  mesh<cell::tetrahedron> m(gen_cube_mesh(1.0, 1.0, 1.0, 10, 10, 10));

  exporter::ensight6_geometry("cube", m);
  return 0;
}
