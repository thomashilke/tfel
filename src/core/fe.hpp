#ifndef _FE_H_
#define _FE_H_

#include "cell.hpp"

template<typename fe_type>
struct is_continuous {
  static const bool value = fe_type::is_continuous;
};

template<typename fe_type>
struct is_lagrangian {
  static const bool value = fe_type::is_lagrangian;
};


namespace cell {
  struct edge::fe::lagrange_p0 {
    typedef cell::edge cell_type;

    static const std::size_t n_dof_per_element = 1;
    static constexpr std::size_t n_dof[2] = {0, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = false;
    static const bool is_lagrangian = true;

    // node coordinates
    static constexpr double x[1][1] = {{0.5}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivatives,
				 const double* x) {
      if (derivatives[0] == 0) {
	return 1.0;
      } else {
	return 0.0;
      }
    }

    static double phi(unsigned int i, const double* x) {
      return 1.0;
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 1)
	throw std::string("finite_element::edge_lagrange_p0::dphi: space dimension is 1,"
			  " therefore k must be in {0}.");
      return 0.0;
    }
  };

  struct edge::fe::lagrange_p1 {
    typedef cell::edge cell_type;

    static const std::size_t n_dof_per_element = 2;
    static constexpr std::size_t n_dof[2] = {1, 0};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = true;

    static constexpr double x[2][1] = {{0.0}, {1.0}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivatives,
				 const double* x) {
      typedef double (*bf_type)(const double*);
      static const bf_type bf[2][2] = {{bf_0_1, bf_0_2},
				       {bf_1_1, bf_1_2}};
      if (i >= 2)
	throw std::string("finite_element::edge_lagrange_p1::basis_function: "
			  "i out of range");

      if (derivatives[0] > 1)
	return 0.0;

      return bf[derivatives[0]][i](x);
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d(0);
      return basis_function(i, &d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k != 0)
	throw std::string("finite_element::edge_lagrange_p1::dphi: space dimension is 1,"
			  " therefore k must be in {0}.");
      unsigned int d(1);
      return basis_function(i, &d, x);
    }

  protected:
    static double bf_0_1(const double* x) { return 1.0 - x[0]; }
    static double bf_0_2(const double* x) { return       x[0]; }
    static double bf_1_1(const double* x) { return - 1.0; }
    static double bf_1_2(const double* x) { return   1.0; }
  };


  struct edge::fe::lagrange_p1_bubble: public edge::fe::lagrange_p1 {
    using edge::fe::lagrange_p1::cell_type;

    static const std::size_t n_dof_per_element = 3;
    static constexpr std::size_t n_dof[] = {1, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = false;

    static constexpr double x[3][1] = {{0.0}, {1.0}, {0.5}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivatives,
				 const double* x) {
      typedef double (*bf_type)(const double*);
      static const bf_type bf[2][3] = {{bf_0_1, bf_0_2, bf_0_3},
				       {bf_1_1, bf_1_2, bf_1_3}};
      if (i >= 3)
	throw std::string("finite_element::edge_lagrange_p1_bubble::basis_function: "
			  "i out of range");

      if (derivatives[0] > 1) {
	if (i < 2)
	  return 0.0;
	else
	  return bf_2_3(x);
      }

      return bf[derivatives[0]][i](x);
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d(0);
      return basis_function(i, &d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k != 0)
	throw std::string("finite_element::edge_lagrange_p1::dphi: space dimension is 1,"
			  " therefore k must be in {0}.");
      unsigned int d(1);
      return basis_function(i, &d, x);
    }

  protected:
    using edge::fe::lagrange_p1::bf_0_1;
    using edge::fe::lagrange_p1::bf_0_2;

    using edge::fe::lagrange_p1::bf_1_1;
    using edge::fe::lagrange_p1::bf_1_2;

    static double bf_0_3(const double* x) {
      return bf_0_1(x) * bf_0_2(x);
    }

    static double bf_1_3(const double* x) {
      return (bf_1_1(x) * bf_0_2(x) +
	      bf_0_1(x) * bf_1_2(x));
    }

    static double bf_2_3(const double* x) {
      return (bf_1_1(x) * bf_1_2(x) +
	      bf_1_1(x) * bf_1_2(x));
    }
  };

  struct edge::fe::lagrange_p2 {
    typedef cell::edge cell_type;

    static const std::size_t n_dof_per_element = 3;
    static constexpr std::size_t n_dof[2] = {1, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = true;

    static constexpr double x[3][1] = {{0.0}, {1.0}, {0.5}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivatives,
				 const double* x) {
      typedef double (*bf_type)(const double*);
      static const bf_type bf[2][3] = {{bf_0_1, bf_0_2, bf_0_3},
				       {bf_1_1, bf_1_2, bf_1_3}};
      if (i >= 3)
	throw std::string("finite_element::edge_lagrange_p1::basis_function: "
			  "i out of range");

      if (derivatives[0] > 1)
	throw std::string("edge::fe::lagrange_p2: higher order derivatives are not implemented.");

      return bf[derivatives[0]][i](x);
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d(0);
      return basis_function(i, &d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k != 0)
	throw std::string("finite_element::edge_lagrange_p1::dphi: space dimension is 1,"
			  " therefore k must be in {0}.");
      unsigned int d(1);
      return basis_function(i, &d, x);
    }

  protected:
    static double bf_0_1(const double* x) { return 2.0 * (x[0] - 0.5) * (x[0] - 1.0); }
    static double bf_0_2(const double* x) { return 2.0 * x[0] * (x[0] - 0.5); }
    static double bf_0_3(const double* x) { return 4.0 * x[0] * (x[0] - 1.0); }
    static double bf_1_1(const double* x) { return 4.0 * x[0] - 3.0; }
    static double bf_1_2(const double* x) { return 4.0 * x[0] - 1.0; }
    static double bf_1_3(const double* x) { return 8.0 * x[0] - 4.0; }
  };


  struct triangle::fe::lagrange_p0 {
    typedef cell::triangle cell_type;

    static const std::size_t n_dof_per_element = 1;
    static constexpr std::size_t n_dof[3] = {0, 0, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = false;
    static const bool is_lagrangian = true;

    static constexpr double x[1][2] = {{1.0 / 3.0, 1.0 / 3.0}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivatives,
				 const double* x) {
      if (derivatives[0] == 0 and derivatives[1] == 0) {
	return 1.0;
      } else {
	return 0.0;
      }
    }

    static double phi(unsigned int i, const double* x) {
      return 1.0;
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 2)
	throw std::string("finite_element::edge_lagrange_p0::dphi: space dimension is 1,"
			  " therefore k must be in {0, 1}.");
      return 0.0;
    }
  };


  struct triangle::fe::lagrange_p1 {
    typedef cell::triangle cell_type;

    static const std::size_t n_dof_per_element = 3;
    static constexpr std::size_t n_dof[] = {1, 0, 0};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = true;

    static constexpr double x[3][2] = {{0.0, 0.0},
				       {1.0, 0.0},
				       {0.0, 1.0}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivative,
				 const double* x) {
      typedef double (*bf_type)(const double*);
      static const bf_type bf[3][3] = {{bf_00_1, bf_00_2, bf_00_3},
				       {bf_10_1, bf_10_2, bf_10_3},
				       {bf_01_1, bf_01_2, bf_01_3}};
      if (derivative[0] == 0 and derivative[1] == 0) {
	return bf[0][i](x);
      } else if(derivative[0] == 1 and  derivative[1] == 0) {
	return bf[1][i](x);
      } else if (derivative[0] == 0 and  derivative[1] == 1) {
	return bf[2][i](x);
      } else {
	return 0.0;
      }
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d[2] = {0, 0};
      return basis_function(i, d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 2)
	throw std::string("finite_element::triangle_lagrange_p1::dphi: space dimension is 2,"
			  " therefore k must be in {0, 1}.");
      unsigned int d[2] = {0, 0};
      d[k] = 1;
      return basis_function(i, d, x);
    }

  protected:
    static double bf_00_1(const double* x) { return 1.0 - x[0] - x[1]; }
    static double bf_00_2(const double* x) { return x[0]; }
    static double bf_00_3(const double* x) { return x[1]; }
    static double bf_10_1(const double* x) { return -1.0; }
    static double bf_10_2(const double* x) { return  1.0; }
    static double bf_10_3(const double* x) { return  0.0; }
    static double bf_01_1(const double* x) { return -1.0; }
    static double bf_01_2(const double* x) { return  0.0; }
    static double bf_01_3(const double* x) { return  1.0; }
  };

  struct triangle::fe::lagrange_p1_bubble: public triangle::fe::lagrange_p1 {
    using triangle::fe::lagrange_p1::bf_00_1;
    using triangle::fe::lagrange_p1::bf_00_2;
    using triangle::fe::lagrange_p1::bf_00_3;
    using triangle::fe::lagrange_p1::bf_10_1;
    using triangle::fe::lagrange_p1::bf_10_2;
    using triangle::fe::lagrange_p1::bf_10_3;
    using triangle::fe::lagrange_p1::bf_01_1;
    using triangle::fe::lagrange_p1::bf_01_2;
    using triangle::fe::lagrange_p1::bf_01_3;


    using triangle::fe::lagrange_p1::cell_type;

    static const std::size_t n_dof_per_element = 4;
    static constexpr std::size_t n_dof[] = {1, 0, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = false;

    static constexpr double x[4][2] = {{0.0, 0.0},
				       {1.0, 0.0},
				       {0.0, 1.0},
				       {1.0 / 3.0, 1.0 / 3.0}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivative,
				 const double* x) {
      typedef double(*bf_type)(const double*);
      static const bf_type bf[3][4] = {{bf_00_1, bf_00_2, bf_00_3, bf_00_4},
				       {bf_10_1, bf_10_2, bf_10_3, bf_10_4},
				       {bf_01_1, bf_01_2, bf_01_3, bf_01_4}};

      if (derivative[0] == 0 and derivative[1] == 0) {
	return bf[0][i](x);
      } else if (derivative[0] == 1 and derivative[1] == 0) {
	return bf[1][i](x);
      } else if (derivative[0] == 0 and derivative[1] == 1) {
	return bf[2][i](x);
      } else {
	throw std::string("triangle_lagrange_p1_bubble: unsupported derivative order.");
      }
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d[2] = {0, 0};
      return basis_function(i, d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 2)
	throw std::string("finite_element::triangle_lagrange_p1::dphi: space dimension is 2,"
			  " therefore k must be in {0, 1}.");
      unsigned int d[2] = {0, 0};
      d[k] = 1;
      return basis_function(i, d, x);
    }

  protected:
    static double bf_00_4(const double* x) { return    27.0 * bf_00_1(x) * bf_00_2(x) * bf_00_3(x); }
    static double bf_10_4(const double* x) { return 27.0 * (  bf_10_1(x) * bf_00_2(x) * bf_00_3(x)
							    + bf_00_1(x) * bf_10_2(x) * bf_00_3(x)
							    + bf_00_1(x) * bf_00_2(x) * bf_10_3(x)); }
    static double bf_01_4(const double* x) { return 27.0 * (  bf_01_1(x) * bf_00_2(x) * bf_00_3(x)
							    + bf_00_1(x) * bf_01_2(x) * bf_00_3(x)
							    + bf_00_1(x) * bf_00_2(x) * bf_01_3(x)); }
  };

  struct triangle::fe::lagrange_p2 {
    typedef cell::triangle cell_type;

    static const std::size_t n_dof_per_element = 6;
    static constexpr std::size_t n_dof[] = {1, 1, 0};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }

    static const bool is_continuous = true;
    static const bool is_lagrangian = true;

    static constexpr double x[6][2] = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
                                       {0.5, 0.0}, {0.0, 0.5}, {0.5, 0.5}};
    static double basis_function(unsigned int i,
                                 const unsigned int* derivative,
                                 const double* x) {
      typedef double (*bf_type)(const double*);
      static const bf_type bf[3][6] = {{bf_00_1, bf_00_2, bf_00_3, bf_00_4, bf_00_5, bf_00_6},
                                       {bf_10_1, bf_10_2, bf_10_3, bf_10_4, bf_10_5, bf_10_6},
                                       {bf_01_1, bf_01_2, bf_01_3, bf_01_4, bf_01_5, bf_01_6}};

      if (derivative[0] == 0 and derivative[1] == 0) {
	return bf[0][i](x);
      } else if(derivative[0] == 1 and  derivative[1] == 0) {
	return bf[1][i](x);
      } else if (derivative[0] == 0 and  derivative[1] == 1) {
	return bf[2][i](x);
      } else {
        throw std::string("triangle::fe::lagrange_p2::basis_function: higher order derivatives not implemented.");
      }
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d[2] = {0, 0};
      return basis_function(i, d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 2)
	throw std::string("finite_element::triangle_lagrange_p2::dphi: space dimension is 2,"
			  " therefore k must be in {0, 1}.");
      unsigned int d[2] = {0, 0};
      d[k] = 1;
      return basis_function(i, d, x);
    }

  protected:
    static double bf_00_1(const double* x) { return   2.0 * (x[0] + x[1] - 1.0) * (x[0] + x[1] - 0.5); }
    static double bf_00_2(const double* x) { return   2.0 * x[0] * (x[0] - 0.5); }
    static double bf_00_3(const double* x) { return   2.0 * x[1] * (x[1] - 0.5); }
    static double bf_00_4(const double* x) { return - 4.0 * x[0] * (x[0] + x[1] - 1.0); }
    static double bf_00_5(const double* x) { return - 4.0 * x[1] * (x[0] + x[1] - 1.0); }
    static double bf_00_6(const double* x) { return   4.0 * x[0] * x[1]; }

    static double bf_10_1(const double* x) { return   4.0 * (x[0] + x[1]) - 3.0; }
    static double bf_10_2(const double* x) { return   4.0 * x[0] - 1.0; }
    static double bf_10_3(const double* x) { return   0.0; }
    static double bf_10_4(const double* x) { return - 4.0 * (2.0 * x[0] + x[1] - 1.0); }
    static double bf_10_5(const double* x) { return - 4.0 * x[1]; }
    static double bf_10_6(const double* x) { return   4.0 * x[1]; }

    static double bf_01_1(const double* x) { return   4.0 * (x[0] + x[1]) - 3.0; }
    static double bf_01_2(const double* x) { return   0.0; }
    static double bf_01_3(const double* x) { return   4.0 * x[1] - 1.0; }
    static double bf_01_4(const double* x) { return - 4.0 * x[0]; }
    static double bf_01_5(const double* x) { return - 4.0 * (2.0 * x[1] + x[0] - 1.0); }
    static double bf_01_6(const double* x) { return   4.0 * x[0]; }
  };

  struct tetrahedron::fe::lagrange_p0 {
    typedef cell::tetrahedron cell_type;

    static const std::size_t n_dof_per_element = 1;
    static constexpr std::size_t n_dof[] = {0, 0, 0, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = false;
    static const bool is_lagrangian = true;

    static constexpr double x[1][3] = {{1.0/4.0, 1.0/4.0, 1.0/4.0}};

    static double basis_function(unsigned int i,
				 unsigned int* derivatives,
				 const double* x) {
      if (derivatives[0] == 0 and derivatives[1] == 0 and derivatives[2] == 0) {
	return 1.0;
      } else {
	return 0.0;
      }
    }

    static double phi(unsigned int i, const double* x) {
      return 1.0;
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 3)
	throw std::string("finite_element::edge_lagrange_p0::dphi: space dimension is 1,"
			  " therefore k must be in {0, 1}.");
      return 0.0;
    }
  };

  struct tetrahedron::fe::lagrange_p1 {
    typedef cell::tetrahedron cell_type;

    static const std::size_t n_dof_per_element = 4;
    static constexpr std::size_t n_dof[] = {1, 0, 0, 0};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = true;

    static constexpr double x[4][3] = {{0.0, 0.0, 0.0},
				       {1.0, 0.0, 0.0},
				       {0.0, 1.0, 0.0},
				       {0.0, 0.0, 1.0}};

    static double basis_function(unsigned int i,
				 const unsigned int* derivative,
				 const double* x) {
      typedef double (*bf_type)(const double*);
      static const bf_type bf[4][4] = {{bf_000_1, bf_000_2, bf_000_3, bf_000_4},
				       {bf_100_1, bf_100_2, bf_100_3, bf_100_4},
				       {bf_010_1, bf_010_2, bf_010_3, bf_010_4},
				       {bf_001_1, bf_001_2, bf_001_3, bf_001_4}};
      if (derivative[0] == 0 and derivative[1] == 0 and derivative[2] == 0) {
	return bf[0][i](x);
      } else if(derivative[0] == 1 and  derivative[1] == 0 and derivative[2] == 0) {
	return bf[1][i](x);
      } else if (derivative[0] == 0 and  derivative[1] == 1 and derivative[2] == 0) {
	return bf[2][i](x);
      } else if (derivative[0] == 0 and  derivative[1] == 0 and derivative[2] == 1) {
	return bf[3][i](x);
      } else {
	return 0.0;
      }
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d[3] = {0, 0, 0};
      return basis_function(i, d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 3)
	throw std::string("finite_element::tetrahedron_lagrange_p1::dphi: space dimension is 3,"
			  " therefore k must be in {0, 2}.");
      unsigned int d[3] = {0, 0, 0};
      d[k] = 1;
      return basis_function(i, d, x);
    }

  protected:
    static double bf_000_1(const double* x) { return 1.0 - x[0] - x[1] - x[2]; }
    static double bf_000_2(const double* x) { return x[0]; }
    static double bf_000_3(const double* x) { return x[1]; }
    static double bf_000_4(const double* x) { return x[2]; }
    static double bf_100_1(const double* x) { return -1.0; }
    static double bf_100_2(const double* x) { return  1.0; }
    static double bf_100_3(const double* x) { return  0.0; }
    static double bf_100_4(const double* x) { return  0.0; }
    static double bf_010_1(const double* x) { return -1.0; }
    static double bf_010_2(const double* x) { return  0.0; }
    static double bf_010_3(const double* x) { return  1.0; }
    static double bf_010_4(const double* x) { return  0.0; }
    static double bf_001_1(const double* x) { return -1.0; }
    static double bf_001_2(const double* x) { return  0.0; }
    static double bf_001_3(const double* x) { return  0.0; }
    static double bf_001_4(const double* x) { return  1.0; }
  };

  struct tetrahedron::fe::lagrange_p1_bubble: public tetrahedron::fe::lagrange_p1 {
    using tetrahedron::fe::lagrange_p1::bf_000_1;
    using tetrahedron::fe::lagrange_p1::bf_000_2;
    using tetrahedron::fe::lagrange_p1::bf_000_3;
    using tetrahedron::fe::lagrange_p1::bf_000_4;

    using tetrahedron::fe::lagrange_p1::bf_100_1;
    using tetrahedron::fe::lagrange_p1::bf_100_2;
    using tetrahedron::fe::lagrange_p1::bf_100_3;
    using tetrahedron::fe::lagrange_p1::bf_100_4;

    using tetrahedron::fe::lagrange_p1::bf_010_1;
    using tetrahedron::fe::lagrange_p1::bf_010_2;
    using tetrahedron::fe::lagrange_p1::bf_010_3;
    using tetrahedron::fe::lagrange_p1::bf_010_4;

    using tetrahedron::fe::lagrange_p1::bf_001_1;
    using tetrahedron::fe::lagrange_p1::bf_001_2;
    using tetrahedron::fe::lagrange_p1::bf_001_3;
    using tetrahedron::fe::lagrange_p1::bf_001_4;

    using tetrahedron::fe::lagrange_p1::cell_type;

    static const std::size_t n_dof_per_element = 5;
    static constexpr std::size_t n_dof[] = {1, 0, 0, 1};
    static constexpr std::size_t n_dof_per_subdomain(unsigned int i) {
      return n_dof[i];
    }
    static const bool is_continuous = true;
    static const bool is_lagrangian = false;

    static constexpr double x[5][3] = {{0.0, 0.0, 0.0},
                                       {1.0, 0.0, 0.0},
                                       {0.0, 1.0, 0.0},
                                       {0.0, 0.0, 1.0},
                                       {1.0/4.0, 1.0/4.0, 1.0/4.0}};

    static double basis_function(unsigned int i,
                                 const unsigned int* derivative,
                                 const double* x) {
      typedef double(*bf_type)(const double* x);
      static const bf_type bf[4][5] = {{bf_000_1, bf_000_2, bf_000_3, bf_000_4, bf_000_5},
                                       {bf_100_1, bf_100_2, bf_100_3, bf_100_4, bf_100_5},
                                       {bf_010_1, bf_010_2, bf_010_3, bf_010_4, bf_010_5},
                                       {bf_001_1, bf_001_2, bf_001_3, bf_001_4, bf_001_5}};

      if        (derivative[0] == 0 and  derivative[1] == 0 and derivative[2] == 0) {
	return bf[0][i](x);
      } else if (derivative[0] == 1 and  derivative[1] == 0 and derivative[2] == 0) {
	return bf[1][i](x);
      } else if (derivative[0] == 0 and  derivative[1] == 1 and derivative[2] == 0) {
	return bf[2][i](x);
      } else if (derivative[0] == 0 and  derivative[1] == 0 and derivative[2] == 1) {
	return bf[3][i](x);
      } else {
	throw std::string("tetrahedron::fe::lagrange_p1_bubble: unsupported derivative order.");
      }
    }

    static double phi(unsigned int i, const double* x) {
      unsigned int d[3] = {0, 0, 0};
      return basis_function(i, d, x);
    }

    static double dphi(unsigned int k, unsigned int i, const double* x) {
      if (k >= 3)
	throw std::string("tetrahedron::fe::lagrange_p1_bubble::dphi: space dimension is 3,"
			  " therefore k must be in {0, 1, 2}.");
      unsigned int d[3] = {0, 0, 0};
      d[k] = 1;
      return basis_function(i, d, x);
    }

  protected:
    static double bf_000_5(const double* x) { return 256.0 * bf_000_1(x) * bf_000_2(x) * bf_000_3(x) * bf_000_4(x); }
    static double bf_100_5(const double* x) { return 256.0 * (  bf_100_1(x) * bf_000_2(x) * bf_000_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_100_2(x) * bf_000_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_000_2(x) * bf_100_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_000_2(x) * bf_000_3(x) * bf_100_4(x)); }
    
    static double bf_010_5(const double* x) { return 256.0 * (  bf_010_1(x) * bf_000_2(x) * bf_000_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_010_2(x) * bf_000_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_000_2(x) * bf_010_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_000_2(x) * bf_000_3(x) * bf_010_4(x)); }
    
    
    static double bf_001_5(const double* x) { return 256.0 * (  bf_001_1(x) * bf_000_2(x) * bf_000_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_001_2(x) * bf_000_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_000_2(x) * bf_001_3(x) * bf_000_4(x)
                                                              + bf_000_1(x) * bf_000_2(x) * bf_000_3(x) * bf_001_4(x)); }
  };
}

#endif /* _FE_H_ */
