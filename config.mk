include site-config.mk

CXX = mpicxx
DEPS_BIN = g++
DEPSFLAGS = -I$(SITE_INCLUDE_DIR) -I$(SITE_PETSC_INCLUDE_DIR) -I$(SITE_LAPACK_INCLUDE_DIR) $(shell mpicxx --showme:compile)
CXXFLAGS += -Wall -Wextra -Wno-unused-parameter -std=c++11 -I$(SITE_INCLUDE_DIR) -I$(SITE_PETSC_INCLUDE_DIR) -I$(SITE_LAPACK_INCLUDE_DIR)
LDFLAGS += -Wall -Wextra -L$(SITE_PETSC_LIB_DIR) -L$(SITE_LAPACK_LIB_DIR) -L./lib/ -L$(SITE_LIB_DIR)
LDLIBS = -lpetsc -llapacke -ltfel -lalucelldb
AR = ar
ARFLAGS = rc
MKDIR = mkdir
MKDIRFLAGS = -p

PREFIX = ~/.local/
BIN_DIR = bin/
INCLUDE_DIR = include/
LIB_DIR = lib/


SOURCES = \
	test/finite_element_space.cpp \
	test/quadrature_1d.cpp \
	test/expression.cpp \
	test/system_assembly.cpp \
	test/linear_solver.cpp \
	test/square_mesh.cpp \
	test/l2_projection.cpp \
	test/triangle_cell.cpp \
	test/fes_cost.cpp \
	test/laplacian_transient.cpp \
	test/basic_finite_element_formulation.cpp \
	test/quadrature_2d.cpp \
	src/core/quadrature.cpp \
	src/core/linear_solver.cpp \
	test/composite_fe.cpp \
	test/stokes_2d.cpp \
	test/l2_p1_bubble_projection.cpp \
	test/alucell_import.cpp \
	test/tetrahedron_cube.cpp \
	test/linear_constraint.cpp \
	test/mesh_data.cpp \
	src/protocols/stokes_2d/driven_cavity.cpp \
	src/protocols/steady_advection_diffusion_2d/step.cpp \
	src/protocols/unsteady_advection_diffusion_2d/rotating_hill.cpp \
	src/protocols/unsteady_advection_diffusion_2d/front.cpp \
	src/protocols/steady_advection_diffusion_1d/stabilisation.cpp \
	src/main.cpp \
	src/core/cell.cpp \
	src/core/mesh.cpp \
	src/core/fe.cpp

HEADERS = \
	include/tfel/tfel.hpp \
	include/tfel/core/basic_fe_formulation.hpp \
	include/tfel/core/bilinear_form.hpp \
	include/tfel/core/cell.hpp \
	include/tfel/core/composite_bilinear_form.hpp \
	include/tfel/core/composite_fe.hpp \
	include/tfel/core/composite_fes.hpp \
	include/tfel/core/composite_form.hpp \
	include/tfel/core/composite_linear_form.hpp \
	include/tfel/core/element_diameter.hpp \
	include/tfel/core/errors.hpp \
	include/tfel/core/export.hpp \
	include/tfel/core/expression.hpp \
	include/tfel/core/fe.hpp \
	include/tfel/core/fes.hpp \
	include/tfel/core/fe_value_manager.hpp \
	include/tfel/core/form.hpp \
	include/tfel/core/linear_algebra.hpp \
	include/tfel/core/linear_form.hpp \
	include/tfel/core/linear_solver.hpp \
	include/tfel/core/mesh.hpp \
	include/tfel/core/meta.hpp \
	include/tfel/core/projector.hpp \
	include/tfel/core/quadrature.hpp \
	include/tfel/core/sparse_linear_system.hpp \
	include/tfel/core/timer.hpp \
	include/tfel/formulations/steady_advection_diffusion_2d.hpp \
	include/tfel/formulations/steady_advection_diffusion_1d.hpp \
	include/tfel/formulations/stokes_2d.hpp \
	include/tfel/formulations/unsteady_advection_diffusion_2d.hpp \
	include/tfel/formulations/unsteady_diffusion_2d.hpp \
	include/tfel/utility/importer.hpp \
	include/tfel/core/vector_operation.hpp \
	include/tfel/core/subdomain.hpp \
	include/tfel/core/mesh_data.hpp \
	include/tfel/core/operator.hpp


BIN = \
	bin/test_finite_element_space \
	bin/test_expression \
	bin/test_quadrature_1d \
	bin/test_system_assembly \
	bin/test_linear_solver \
	bin/test_square_mesh \
	bin/test_l2_projection \
	bin/test_triangle_cell \
	bin/test_fes_cost \
	bin/test_laplacian_transient \
	bin/test_basic_finite_element_formulation \
	bin/test_composite_fe \
	bin/test_quadrature_2d \
	bin/test_stokes_2d \
	bin/test_tetrahedron_cube \
	bin/test_l2_p1_bubble_projection \
	bin/test_linear_constraint \
	bin/test_mesh_data \
	bin/prot_stoke_2d_driven_cavity \
	bin/prot_steady_advection_diffusion_2d_step \
	bin/prot_unsteady_advection_diffusion_2d_rotating_hill \
	bin/prot_unsteady_advection_diffusion_2d_front \
	bin/prot_steady_advection_diffusion_1d_stabilisation \
	bin/test_alucell_import \
	bin/main

bin/test_finite_element_space: build/test/finite_element_space.o 
bin/main: build/src/main.o 
bin/test_quadrature_1d: build/test/quadrature_1d.o
bin/test_expression: build/test/expression.o
bin/test_system_assembly: build/test/system_assembly.o
bin/test_linear_solver: build/test/linear_solver.o
bin/test_square_mesh: build/test/square_mesh.o
bin/test_l2_projection: build/test/l2_projection.o
bin/test_triangle_cell: build/test/triangle_cell.o
bin/test_fes_cost: build/test/fes_cost.o
bin/test_laplacian_transient: build/test/laplacian_transient.o
bin/test_basic_finite_element_formulation: build/test/basic_finite_element_formulation.o
bin/test_composite_fe: build/test/composite_fe.o
bin/test_quadrature_2d: build/test/quadrature_2d.o
bin/test_stokes_2d: build/test/stokes_2d.o
bin/test_l2_p1_bubble_projection: build/test/l2_p1_bubble_projection.o
bin/prot_stoke_2d_driven_cavity: build/src/protocols/stokes_2d/driven_cavity.o
bin/prot_steady_advection_diffusion_2d_step: build/src/protocols/steady_advection_diffusion_2d/step.o
bin/prot_unsteady_advection_diffusion_2d_rotating_hill: build/src/protocols/unsteady_advection_diffusion_2d/rotating_hill.o
bin/prot_unsteady_advection_diffusion_2d_front: build/src/protocols/unsteady_advection_diffusion_2d/front.o
bin/prot_steady_advection_diffusion_1d_stabilisation: build/src/protocols/steady_advection_diffusion_1d/stabilisation.o
bin/test_alucell_import: build/test/alucell_import.o
bin/test_tetrahedron_cube: build/test/tetrahedron_cube.o
bin/test_linear_constraint: build/test/linear_constraint.o
bin/test_mesh_data: build/test/mesh_data.o

LIB = lib/libtfel.a

lib/libtfel.a: \
	build/src/core/fe.o \
	build/src/core/mesh.o \
	build/src/core/quadrature.o \
	build/src/core/cell.o \
	build/src/core/linear_solver.o
