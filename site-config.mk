ifndef PETSC_DIR
  $(error PETSC_DIR environment variable is not set)
endif

SITE_INCLUDE_DIR=$(HOME)/.local/include/
SITE_LIB_DIR=$(HOME)/.local/lib/

SITE_LAPACK_INCLUDE_DIR=/usr/include/
SITE_LAPACK_LIB_DIR=/usr/lib/

SITE_PETSC_INCLUDE_DIR=$(PETSC_DIR)/include
SITE_PETSC_LIB_DIR=$(PETSC_DIR)/lib

CXXFLAGS = -O2
LDFLAGS = -O2

WITH_ALUCELL = on
WITH_FREEFEM = off
WITH_GMSH = off
