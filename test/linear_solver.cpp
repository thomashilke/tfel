#include "../src/core/linear_solver.hpp"

int main(int argc, char *argv[]) {
  linear_solver s;
  bool export_to_stdout(false);

  unsigned int n(10000);
  auto petsc_gmres_ilu(s.get_solver(solver::petsc, method::gmres, preconditioner::ilu));
  petsc_gmres_ilu->set_size(n);

  std::vector<int> nz(n, 3);
  petsc_gmres_ilu->preallocate(&nz[0]);
  
  const double h(1.0 / n);
  
  array<double> rhs{n};
  for(unsigned int i(1); i < n - 1; ++i) {
    petsc_gmres_ilu->add_value(i, i-1, -1.0);
    petsc_gmres_ilu->add_value(i, i, 2.0);
    petsc_gmres_ilu->add_value(i, i+1, -1.0);
    rhs.at(i) = 1.0 * h * h;
  }
  rhs.at(0) = 1.0 * h * h;
  petsc_gmres_ilu->add_value(0, 0, 2.0);
  petsc_gmres_ilu->add_value(0, 1, -1.0);

  rhs.at(n - 1) = 1.0 * h * h;
  petsc_gmres_ilu->add_value(n - 1, n - 2, -1.0);
  petsc_gmres_ilu->add_value(n - 1, n - 1, 2.0);
  
  petsc_gmres_ilu->assemble();

  array<double> x(petsc_gmres_ilu->solve(rhs));

  if (export_to_stdout) {
    std::cout.precision(12);
    for(std::size_t i(0); i < n; ++i)
      std::cout << i * h << ' ' << x.at(i) << '\n';
  }

  delete petsc_gmres_ilu;
  solver::petsc::initialize::release();
  
  return 0;
}
