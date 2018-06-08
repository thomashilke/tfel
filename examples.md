---
layout: page
title: Examples
permalink: /examples/
toc: true
---

The strength of the TFEL is to provide a DSL close to the mathematical
language typically used to express and describe linear partial
differential equations and their solutions.

This page go through several examples, going form the mathematical
description of continuous problem to it's discretization and
implementation with the help of the TFEL.

We start from the most basic problem such as the Poisson equation with
Dirichlet boundary conditions to more advanced problems related to
transport equations and fluid flows.

# Poisson problem with Dirichlet boundary conditions
In this example we solve the Poisson equation in a square in the 2D
plane. The Poisson equation can be thought as a model for the elastic
deformation a thin membrane under a surfacic load. The problem is stated as
follow. Let \\( \Omega \\) be the area occupied by the membrane. For
simplicity we well consider a square membrane:

$$
  \Omega = (0, 1)^2 \subset \mathbb R^2, 
$$

Let \\( f:\Omega\to\mathbb R \\) a function in \\( \mathrm L^2(\Omega)
\\). The values of the function \\( f \\) model the force per unit
surface applied on the membrane. Then, the vertical displacment of the
membrane \\( u:\Omega\to\mathbb R \\) is a solution of the Poisson
equation

$$
  -\Delta\, u(x) = f(x) \quad \forall x\in \Omega.
$$

For the problem to be well-posed, it is well known that we need
additional boundary conditions. Here we consider the homogeneous Dirichlet
boundary conditions 

$$
  u(x) = 0 \quad \forall x \in \partial\Omega,
$$

where \\( \partial\Omega \\) is the boundary of the open set \\(
\Omega \\).

We now give the weak formulation of the Poisson problem. We are
looking for the vertical displacement \\( u \in H\_0^1(\Omega) \\) such
that 

$$
  \int_\Omega \nabla u \cdot \nabla v \,\mathrm dx = \int_\Omega f v
  \,\mathrm dx \tag{1}\label{eq:weak-form}
$$

for all test function \\( v \in H_0^1(\Omega) \\). The functional
space \\( H_0^1(\Omega) \\) is the subspace of the Sobolev space \\(
W^{1, 2} \\) whose elements vanish on the boundary of \\( \Omega \\):

$$
  H^1_0(\Omega) = \left\{ v \in W^{1,2}\ \mid\ v = 0 \text{ on } \partial \Omega  \right\}.
$$

Expanding the intregand, the expression \\(\ref{eq:weak-form}\\) is
equivalent to

$$
  \int_\Omega \displaystyle\frac{\partial u}{\partial x_1}\displaystyle\frac{\partial v}{\partial x_1} + \displaystyle\frac{\partial u}{\partial x_2}\displaystyle\frac{\partial v}{\partial x_2} \,\mathrm dx = \int_\Omega f v \,\mathrm dx \tag{2}\label{eq:weak-form-exp}
$$

We now proceed to the discretization of the weak formulation of the
Poisson problem \\(\ref{eq:weak-form}\\). Let \\( \mathcal M_h \\) be a conformal triangular mesh of the domain
\\( \Omega \\), where \\( h \\) is the largest element diameter:

$$
  h = \max_{\mathcal K \in\ M_h} \mathrm{diam}(\mathcal K),
$$

and let \\( V_h \\) be the function space

$$
  V_h = \left\{ v \in H_0^1(\Omega) \mid v|_\mathcal K \in \mathbb P_1(\mathcal K)\ \forall \mathcal K \in \mathcal M_h\right\}.
$$

We obtain a discrete version of the weak problem with the following
Galerkin approximation. The discrete problem is a follow. We look for
a function \\( u_h \in V_h \\) such that

$$
  \int_\Omega \nabla u_h \cdot \nabla v_h \,\mathrm dx =
  \int_\Omega f v_h \,\mathrm dx
$$

for all test function \\( v_h \in V_h \\). It is well known that this
last problem is well-posed and that the function \\( u_h \\) is an
approximation of \\( u \\) in the sense that

$$
    \left\lVert u - u_h \right\rVert = O(h^2).
$$

We know proceed to the implementation of discrete weak formulation
with the help of the TFEL. All the following code should go in a single file,
called `poisson.cpp` for example.

We start by including the library header. For clarity, the function
`f` is only declared here, and will be defined later in the
source. The constant `n` hold the number of subdivision for the
discretization of \\(\Omega\\). We then open the `main` code block.
{% highlight c++ linenos %}
#include <tfel/tfel.hpp>

double f(const double* x);

const std::size_t n = 10;

int main() {
{% endhighlight %}
The value of the function `f` is the force that we apply on
the membrane, as a function of space coordinates given in the argument
`x`. Spaces coordinates are always `const double*` in the parameter
list of functions.

We know need to declare the type of mesh that we wish to use, the type
of finite element and finite element space that will be used for the
discretization of the weak form. In the TFEL, all these model
parameters are represented by types and combinaison of types through
templates parameters. For clarity, it is a good idea to declare
aliases to these new types as follow:
{% highlight c++ linenos %}
  using cell_type = cell::triangle;
  using mesh_type = fe_mesh<cell_type>;
  using fe_type = cell_type::fe::lagrange_p1;
  using fes_type = finite_element_space<fe_type>;
  using quad_type = quad::triangle::qf5pT;
{% endhighlight %}
Here we declare `cell_type`, the type of cells we use to discretize
the domain \\(\Omega\\) (triangles), `mesh_type` the type of mesh
datatype which is parametrized by the cell type, `fe_type` the finite
element we picked (continuous, linear and Lagrange with nodes on the
vertices of the triangles), `fes_type` the finite element space type
which describe the space where the solution and test functions live,
and `quad_type` the quadrature formula that will be used to integrate
the linear and bilinear forms.

We can then declare a mesh and its boundary submesh:
{% highlight c++ linenos %}
  mesh_type m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());
{% endhighlight %}
We use the utility function `gen_square_mesh` to generate the
structured mesh of a rectangle with side \\((1,1)\\) into \\(2n^2 =
20\\) triangles. The variable `dm` hold the boundary submesh of `m`.

We then allocate a finite element space based on the mesh `m` with the
finite element `fe_type`, and specify the Dirichlet boundary:
{% highlight c++ linenos %}
  fes_type fes(m);
  fes.add_dirichlet_boundary(dm);
{% endhighlight %}
The `add_dirichlet_boundary` method of the class `fes_type` specify a
homogeneous Dirichlet condition by default, but arbitraty function
defined on \\(\partial\Omega\\) can be specified.

We can now proceed to the allocation and integration of the bilinear
and linear forms:
{% highlight c++ linenos %}
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
{% endhighlight %}
On the lines 1 and 10 we allocate bilinear and linear forms
parametrized by the finite element space `fes_type`. For convenience
and to avoid name conflicts for the test and trial functions in
particular we wrap the definition of each forms in a code block
`{...}`.

Lines 2, 3 and 11 declare instances of the test and trial functions of
the problem, and which are used in the expressions in the integrands of the
integrals that define the linear and bilinear forms. We then use a
natural syntax to express the expression to be integrated on the lines
5-7 and 13. The first argument of the `integrate<quad_type>` function
call is an expression to be integrated over the domain discretized by
the mesh `m`. This expression is a function of the variables `u` (for
the bilinear form) and `v`, and must be linear in these variables
(although no test is done at compile time or runtime). Factors and
terms can be composed with the usual `C++` operators `+`, `-` and
`*`. Derivatives of the variables `u` and `v` are expressed with the
function `d<position>(u)` where `position` is a `std::size_t` and
represent the argument of `u` to derive with respect to. For exemple
`d<1>(u)` is correspond to \\(\frac{\partial u}{\partial x\_1}\\),
where \\(x_1\\) is the first argument of the function \\(u\\). You can
compare the expression on the lines 5-7 and 13 with the integrands in
equation \\(\ref{eq:weak-form-exp}\\).
This notation allow for a high readability of the code and help to
prevent most mistakes in the assembly routines.

Once the two forms are built, we can solve the resulting linear
system. So far two linear solver are available: a dense matrix LU
solver which call LAPACK routines, and a sparse matrix GMRES solver
with an ILU preconditioner implemented with PETSc.

{% highlight c++ linenos %}
  dictionary param(dictionary()
                     .set("maxits",  2000u)
                     .set("restart", 1000u)
                     .set("rtol",    1.e-8)
                     .set("atol",    1.e-50)
                     .set("dtol",    1.e20)
                     .set("ilufill", 2u));
  solver::petsc::gmres_ilu s(param);
    
  const fes_type::element u(a.solve(b, s));
{% endhighlight %}
Here we first declare a `dictionary` which keys are the parameter of
the solver declare on line 7. Then we call the `solve` method of the
bilinear form with the linear form `f` and the solver `s` as
parameters.

The solution `u` is a finite element function, that is an element of
the finite element space `fes`. The solution can eventually be
exported for visualization in the Ensight6 file format and return:
{% highlight c++ linenos %}
  exporter::ensight6("poisson",
                     to_mesh_vertex_data<fe_type>(u), "u");

  return 0;
}
{% endhighlight %}
Since most file format only support data defined on the cell or on the
vertices of a mesh, we first need to convert the finite element
function `u` to a `mesh_data<double, mesh_type>`. It is done
temporarily on line 2.

The last piece missing is the definition of the function `f`. For
simplicity, we choose a constant function \\( f(x) = 1\ \forall
x\in\Omega\\):
{% highlight c++ linenos %}
double f(const double* x) {
  return 1.0;
}
{% endhighlight %}

The full source code corresponding to this example fits in just under
55 lines, and is ready to be copy-pasted here:
{% highlight c++ linenos %}
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
{% endhighlight %}

The file `poisson.cpp` can be compiled with the command
```bash
$ clang++ -std=c++14 poisson.cpp -I$HOME/.local/include -I/usr//local/Cellar/lapack/3.7.1/include/ -I$PETSC_DIR/include -L$HOME/.local/lib -L$PETSC_DIR/lib -L/usr/local/Cellar/lapack/3.7.1/lib -llapacke -ltfel  -lpetsc -lmpi -o poisson
```
(can vary depending on your setup) and run:
```bash
$ ./poisson
```
The result is stored in an Ensight6 case file, which can be visualized with Paraview. Below is a screeshot of the solution where the color/x3-coordinate is the value of \\(u\\), the displacment of the membrane under the load \\(f = 1\\) on \\(\Omega\\).
![Solution to the Poisson problem]({{ "/assets/poisson.png" | absolute_url }})


# Transport equation with an SUPG numerical stabilisation method
ToDo

# Stokes problem with Dirichlet boundary conditions
ToDo

# Navier-Stokes problem with Dirichlet boundary conditions
ToDo
