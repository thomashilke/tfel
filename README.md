# Template Finite Element Library
By Thomas Hilke. More informations can be found at
[https://thomashilke.github.io/tfel/](https://thomashilke.github.io/tfel/).

## Requirements
To compile an use the TFEL, you must have the following packages which
are probably available through your OS package system:

 - OpenMPI v2.1.1
 - liblapacke-dev v3.5.0
 - PETSc 3.7.0
 - clang++ 3.4
 
and the following packages must be compiled and installed:

 - https://github.com/thomashilke/alucell-db-tools.git
 - https://github.com/thomashilke/spikes.git

Location of the PETSc installation and `liblapacke-dev` header files
can be configured in the file `site-config.mk`. Run `make` in the root
directory of the project to build the static library and a set of
binaries, which are stored in the subdirectory bin/. Each binary
uses and tests some feature of the library. Run `make install` to
install the binaries, static library and header files in the
directories ~/.local/{bin/,lib/,include/tfel/}.

Some knowledge of the finite element theory and the intricacies of the
implementation of these methods helps.

## Brief Intro

The TFEL is a set of tools to help in rapid prototyping of numerical
algorithms built around finite element methods. The TFEL provide a
domain specific language which matches closely the mathematical
notation of finite element problem formulations. The TFEL models the
finite element method concepts using the C++ type system and template
metaprogramming techniques to offer a user-friendly interface an
prevent compilation of expressions and constructs that are
mathematically non-sensical.

## Features

The TFEL is a set of building blocks that are organized to be easily
extensible.  Nevertheless the toolbox is already rich enough to
implement and solve some not-so-trivial-problems, such as the
Navier-Stokes equations.

The following concepts and features, among others, are already available:

 - Conformal, simplicial meshes built from edges, triangle, tetrahedron cells,
 - Polynomial Lagrange finite element of degree 0, 1 and 2 on edges, triangle, tetrahedron cells,
 - Scalar and vector finite element spaces,
 - Linear and bilinear forms expressed as integrals over meshes and submeshes,
 - Various numerical quadrature formulas,
 - Interface with PETSc's sparse solvers, and LAPACK dense solver,
 - Export in Ensight6 file format,

## Hello World: The Poisson Equation in 2D
One of the simplest elliptical partial differential equation is the
 [Poisson's equation](https://en.wikipedia.org/wiki/Poisson%27s_equation). The
following example solve the Poisson's equation with homogeneous
boundary conditions and a constant right hand side on a square domain
in 2 space dimensions.

```c++
#include <tfel/tfel.hpp>

double f(const double* x);

const std::size_t n = 10;

int main() {
  using cell_type = cell::triangle;
  using mesh_type = fe_mesh<cell_type>;
  using fe_type = cell_type::fe::lagrange_p1;
  using fes_type = finite_element_space<fe_type>;
  using quad_type = quad::triangle::qf5pT;
  
  mesh_type m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());

  fes_type fes(m);
  fes.add_dirichlet_boundary(dm);

  bilinear_form<fes_type, fes_type> a(fes, fes); {
    auto u(a.get_trial_function());
    auto v(a.get_test_function());
  
    a += integrate<quad_type>(d<1>(u) * d<1>(v) +
                              d<2>(u) * d<2>(v)
                              , m);
  }
  
  linear_form<fes_type> b(fes); {
    auto v(b.get_test_function());
    
    b += integrate<quad_type>(f * v, m);
  }
  
  dictionary param(dictionary()
                     .set("maxits",  2000u)
                     .set("restart", 1000u)
                     .set("rtol",    1.e-8)
                     .set("atol",    1.e-50)
                     .set("dtol",    1.e20)
                     .set("ilufill", 2u));
  solver::petsc::gmres_ilu s(param);
    
  const fes_type::element u(a.solve(b, s));

  exporter::ensight6("poisson",
                     to_mesh_vertex_data<fe_type>(u), "u");

  return 0;
}

double f(const double* x) {
  return 1.0;
}
```

Assuming the code is saved in a file poisson.cpp, it can be compiled with the command

```shell
clang++ -std=c++14 poisson.cpp -I$HOME/.local/include -I/usr//local/Cellar/lapack/3.7.1/include/ -I$PETSC_DIR/include -L$HOME/.local/lib -L$PETSC_DIR/lib -L/usr/local/Cellar/lapack/3.7.1/lib -llapacke -ltfel  -lpetsc -lmpi -o poisson
```

(Your setup may vary. Here LAPACK is installed through homebrew on macos, and PETSc is compiled from source)

Run the computation with

```shell
./poisson
```

and the results are stored in the output file poisson.case.
