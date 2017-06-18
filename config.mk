CXX = g++
DEPS_BIN = g++
DEPS_FLAGS = -I/Users/thomashilke/.local/include/
CXXFLAGS = -g -std=c++11 -I/Users/thomashilke/.local/include/
LDFLAGS = -g
LDLIBS = 
AR = ar
ARFLAGS = rc

PREFIX = ~/.local/
BIN_DIR = bin/
INCLUDE_DIR = include/
LIB_DIR = lib/


SOURCES = src/main.cpp

HEADERS = 

BIN = bin/main
bin/main: build/src/main.o


LIB = 

#lib/...: ...
