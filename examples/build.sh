#!/bin/bash

CXX=clang++
CXXFLAGS="-std=c++11 -O3 -flto -I/Users/thomashilke/.local/include/ -I/usr/local/opt/lapack/include/"
LDLIB=" -lalucelldb -llapacke -lpetsc -ltfel -lmpi "
LDFLAGS=" -L/Users/thomashilke/.local/lib/ -L/usr/local/opt/lapack/lib/"

$CXX $CXXFLAGS $LDFLAGS -o poisson poisson.cpp $LDLIB
