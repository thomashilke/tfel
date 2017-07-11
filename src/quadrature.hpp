#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include <cstddef>

#include "cell.hpp"

namespace quad {
  namespace point {
    struct eval {
      typedef cell::point cell_type;
      
      static const std::size_t n_point_space_dimension = 1;
      static const std::size_t n_point = 1;
      static const double x[1][1];
      static const double w[1];
    };
  }

  namespace edge {
    // see wikipedia on gauss quadrature rules
    struct gauss1 {
      typedef cell::edge cell_type;
      
      static const std::size_t n_point_space_dimension = 1;
      static const std::size_t n_point = 1;
      static const double x[1][1];
      static const double w[1];
    };
  
    struct gauss2 {
      typedef cell::edge cell_type;
      
      static const std::size_t n_point_space_dimension = 1;
      static const std::size_t n_point = 2;
      static const double x[2][1];
      static const double w[2];
    };

    struct gauss3 {
      typedef cell::edge cell_type;
      
      static const std::size_t n_point_space_dimension = 1;
      static const std::size_t n_point = 3;
      static const double x[3][1];
      static const double w[3];
    };

    struct gauss4 {
      typedef cell::edge cell_type;
      
      static const std::size_t n_point_space_dimension = 1;
      static const std::size_t n_point = 4;
      static const double x[4][1];
      static const double w[4];
    };

    struct gauss5 {
      typedef cell::edge cell_type;
      
      static const std::size_t n_point_space_dimension = 1;
      static const std::size_t n_point = 5;
      static const double x[5][1];
      static const double w[5];
    };
  }

  namespace triangle {
    // see freefem++ manual, pp190
    struct qf1pT {
      typedef cell::triangle cell_type;
      
      static const std::size_t n_point_space_dimension = 2;
      static const std::size_t n_point = 1;
      static const double x[1][2];
      static const double w[1];
    };

    struct qf2pT {
      typedef cell::triangle cell_type;
      
      static const std::size_t n_point_space_dimension = 2;
      static const std::size_t n_point = 3;
      static const double x[3][2];
      static const double w[3];
    };

    struct qf5pT {
      typedef cell::triangle cell_type;
      
      static const std::size_t n_point_space_dimension = 2;
      static const std::size_t n_point = 7;
      static const double x[7][2];
      static const double w[7];
    };

    struct qf1pTlump {
      typedef cell::triangle cell_type;
      
      static const std::size_t n_point_space_dimension = 2;
      static const std::size_t n_point = 3;
      static const double x[3][2];
      static const double w[3];
    };
  }
};

template<typename Q, typename F>
double integrate(const F& f, double a, double b, double n) {
  double result(0.0);

  const double h((b - a) / n);
  for (unsigned int k(0); k < n; ++k) {
    const double x1(k * h), x2((k + 1) * h);
    for (unsigned int q(0); q < Q::n_point; ++q)
      result += (x2 - x1) * Q::w[q] * f(x1 + (x2 - x1) * Q::x[q][0]);
  }

  return result;
}

#endif /* _QUADRATURE_H_ */
