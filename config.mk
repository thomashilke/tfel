include site-config.mk

export OMPI_CXX = clang++
CXX = mpicxx
DEPS_BIN = g++
DEPSFLAGS = -I$(SITE_INCLUDE_DIR) -I$(SITE_PETSC_INCLUDE_DIR) -I$(SITE_LAPACK_INCLUDE_DIR) $(shell mpicxx --showme:compile)
CXXFLAGS = -g -Wall -Wextra -std=c++11 -I$(SITE_INCLUDE_DIR) -I$(SITE_PETSC_INCLUDE_DIR) -I$(SITE_LAPACK_INCLUDE_DIR)
LDFLAGS = -g -Wall -Wextra -L$(SITE_PETSC_LIB_DIR) -L$(SITE_LAPACK_LIB_DIR)
LDLIBS = -lpetsc -llapacke
AR = ar
ARFLAGS = rc

PREFIX = ~/.local/
BIN_DIR = bin/
INCLUDE_DIR = include/
LIB_DIR = lib/


SOURCES = test/finite_element_space.cpp \
	test/quadrature_1d.cpp \
	test/expression.cpp \
	test/system_assembly.cpp \
	test/linear_solver.cpp \
	test/square_mesh.cpp \
	test/l2_projection.cpp \
	test/triangle_cell.cpp \
	src/quadrature.cpp \
	src/main.cpp \
	src/cell.cpp \
	src/mesh.cpp \
	src/fe.cpp

HEADERS = 

BIN = bin/test_finite_element_space \
	bin/test_expression \
	bin/test_quadrature_1d \
	bin/test_system_assembly \
	bin/test_linear_solver \
	bin/test_square_mesh \
	bin/test_l2_projection \
	bin/test_triangle_cell \
	bin/main
bin/test_finite_element_space: build/test/finite_element_space.o build/src/fe.o
bin/main: build/src/main.o build/src/quadrature.o build/src/cell.o build/src/mesh.o build/src/fe.o
bin/test_quadrature_1d: build/test/quadrature_1d.o build/src/quadrature.o
bin/test_expression: build/test/expression.o
bin/test_system_assembly: build/test/system_assembly.o build/src/quadrature.o build/src/mesh.o build/src/fe.o
bin/test_linear_solver: build/test/linear_solver.o build/src/fe.o
bin/test_square_mesh: build/test/square_mesh.o build/src/mesh.o
bin/test_l2_projection: build/test/l2_projection.o build/src/quadrature.o build/src/cell.o build/src/mesh.o build/src/fe.o
bin/test_triangle_cell: build/test/triangle_cell.o build/src/mesh.o build/src/quadrature.o

LIB = 

#lib/...: ...
