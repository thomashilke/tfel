CXX = g++
DEPS_BIN = g++
CXXFLAGS = -g -std=c++11 -I/home/foetisch/.local/include/
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
