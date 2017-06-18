CXX = g++
DEPS_BIN = g++
DEPS_FLAGS = -I/home/thomas/.local/include/
CXXFLAGS = -g -std=c++11 -I/home/thomas/.local/include/
LDFLAGS = -g
LDLIBS = 
AR = ar
ARFLAGS = rc

PREFIX = ~/.local/
BIN_DIR = bin/
INCLUDE_DIR = include/
LIB_DIR = lib/


SOURCES = test/finite_element_space.cpp \
	src/main.cpp

HEADERS = 

BIN = bin/test_finite_element_space \
	bin/main
bin/test_finite_element_space: build/test/finite_element_space.o
bin/main: build/src/main.o

LIB = 

#lib/...: ...
