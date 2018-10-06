---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults
layout: home
---

# Introduction
TFEL is a C++ library which provide a framework and basic tools to
implement and experiment with finite element methods and solve partial
differential equations. TFEL leverage C++ template mechanisms and the
C++ type system to assemble algorithms and help avoid common mistakes
while experimenting with numerical algorithms.

This project started as an experiment and its development is a
hobby. It's not meant to be used in a production environment, and I
provide absolutly no guarantees of correctness, or that it won't kill
your cat(s).

The TFEL does not try to compete with peta-scale-massively-parallel
codes which run on hundreds of thousands of cluster cores. It was
meant for small-scale experiments that could have been written in
MatLAB, but runs with the performances of a compiled language.

# Development

I have a few stuff on my current TODO list:

  * Move FES element definition out of the FES namespace.
  * Wrap all large array allocations in C++11 smart pointers.
  * Avoid evaluation of null block in finite element system matrices.

However, there is no guarantee any of these items will ever be
checked, as it depends on the time I can allocate for this project.
