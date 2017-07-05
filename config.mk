include site-config.mk

CXX = g++
DEPS_BIN = g++
DEPSFLAGS = -I$(SITE_INCLUDE_DIR)
CXXFLAGS = -g -std=c++11 -I$(SITE_INCLUDE_DIR)
LDFLAGS = -g
LDLIBS = 
AR = ar
ARFLAGS = rc

PREFIX = ~/.local/
BIN_DIR = bin/
INCLUDE_DIR = include/
LIB_DIR = lib/


SOURCES = test/finite_element_space.cpp \
	test/quadrature_1d.cpp \
	src/quadrature.cpp \
	src/main.cpp \
	src/cell.cpp

HEADERS = 

BIN = bin/test_finite_element_space \
	bin/test_quadrature_1d \
	bin/main
bin/test_finite_element_space: build/test/finite_element_space.o
bin/main: build/src/main.o build/src/quadrature.o build/src/cell.o
bin/test_quadrature_1d: build/test/quadrature_1d.o build/src/quadrature.o

LIB = 

#lib/...: ...
