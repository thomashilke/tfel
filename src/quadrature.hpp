#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include <cstddef>

namespace quad {
  namespace point {
    struct eval {
      static const std::size_t number_of_points = 1;
      static const double x[1];
      static const double w[1];
    };
  }

  namespace edge {
    // see wikipedia on gauss quadrature rules
    struct gauss1 {
      static const std::size_t number_of_points = 1;
      static const double x[1];
      static const double w[1];
    };
  
    struct gauss2 {
      static const std::size_t number_of_points = 2;
      static const double x[2];
      static const double w[2];
    };

    struct gauss3 {
      static const std::size_t number_of_points = 3;
      static const double x[3];
      static const double w[3];
    };

    struct gauss4 {
      static const std::size_t number_of_points = 4;
      static const double x[4];
      static const double w[4];
    };

    struct gauss5 {
      static const std::size_t number_of_points = 5;
      static const double x[5];
      static const double w[5];
    };
  }

  namespace triangle {
    // see freefem++ manual, pp190
    struct qf1pT {
      static const std::size_t number_of_points = 1;
      static const double x[1][2];
      static const double w[1];
    };

    struct qf2pT {
      static const std::size_t number_of_points = 3;
      static const double x[3][2];
      static const double w[3];
    };

    struct qf5pT {
      static const std::size_t number_of_points = 7;
      static const double x[7][2];
      static const double w[7];
    };

    struct qf1pTlump {
      static const std::size_t number_of_points = 3;
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
    for (unsigned int q(0); q < Q::number_of_points; ++q)
      result += (x2 - x1) / 2.0 * Q::w[q] * f((x1 + x2)/2.0 + (x2 - x1)/2.0 * Q::x[q]);
  }

  return result;
}

#endif /* _QUADRATURE_H_ */
