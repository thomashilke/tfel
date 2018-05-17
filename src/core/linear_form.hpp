#ifndef _LINEAR_FORM_H_
#define _LINEAR_FORM_H_

#include <spikes/thread_pool.hpp>

template<typename test_fes_type>
class linear_form {
public:
  linear_form(const test_fes_type& te_fes,
              std::size_t algebraic_dof_number = 0,
	      std::size_t n_thread = std::thread::hardware_concurrency())
    : tp(n_thread),
      test_fes(te_fes),
      f{te_fes.get_dof_number()},
      constraint_values(algebraic_dof_number, 0.0) {
    f.fill(0.0);
  }

  template<typename T>
  void operator+=(T integration_proxy) {
    typedef typename test_fes_type::fe_type test_fe_type;
    const std::size_t n_test_dof(test_fe_type::n_dof_per_element);

    std::size_t n_element(integration_proxy.m.get_cell_number());
    std::size_t n_thread(tp.size());

    std::vector<array<double> > rhs_els(n_thread, array<double>{});
    for (std::size_t n(0); n < n_thread; ++n) {
      std::size_t
	k_begin(n_element * n / n_thread ),
	k_end(std::min(n_element * (n + 1) / n_thread, n_element));

      rhs_els[n] = array<double>{k_end - k_begin, n_test_dof};
    }

    std::vector<std::future<void> > futures;
    for (std::size_t n(0); n < n_thread; ++n) {
      std::size_t
	k_begin(n_element * n / n_thread ),
	k_end(std::min(n_element * (n + 1) / n_thread, n_element));
      futures.push_back(
        tp.enqueue(
	  [this, &rhs_els, integration_proxy, k_begin, k_end, n]
	  () {
	    this->assemble_element_range<T>(rhs_els[n],	k_begin, k_end, integration_proxy);
	  }
        )
      );

      //assemble_element_range<T>(rhs_els[n], k_begin, k_end, integration_proxy);
    }

    for (auto& f: futures) f.wait();

    for (std::size_t n(0); n < n_thread; ++n) {
      std::size_t
	k_begin(n_element * n / n_thread ),
	k_end(std::min(n_element * (n + 1) / n_thread, n_element));

      for (std::size_t k(k_begin); k < k_end; ++k)
	for (std::size_t j(0); j < n_test_dof; ++j)
	  f.at(test_fes.get_dof(integration_proxy.get_global_cell_id(k), j)) += rhs_els[n].at(k - k_begin, j);
    }
  }

  template<typename T>
  void assemble_element_range(array<double>& rhs_el,
			      std::size_t k_begin, std::size_t k_end,
			      T integration_proxy) {
    rhs_el.fill(0.0);

    static_assert(T::form_type::rank == 1, "linear_form expects rank-1 expression.");

    typedef typename test_fes_type::fe_type test_fe_type;
    using test_fe_type = typename test_fes_type::fe_type;
    typedef typename T::quadrature_type quadrature_type;
    typedef typename T::cell_type cell_type;
    using form_type = typename T::form_type;

    const auto& m(integration_proxy.m);


    // prepare the quadrature weights
    const std::size_t n_q(quadrature_type::n_point);
    array<double> omega{n_q};
    omega.set_data(&quadrature_type::w[0]);

    // storage for the quadrature points
    array<double> xq_hat{n_q, test_fe_type::cell_type::n_dimension};
    array<double> xq{n_q, test_fe_type::cell_type::n_dimension};

    // storage for the point-wise basis function evaluation
    using fe_list = type_list<test_fe_type>;
    using unique_fe_list = unique_t<fe_list>;

    const std::size_t test_fe_index(get_index_of_element<test_fe_type, unique_fe_list>::value);
    fe_value_manager<unique_fe_list> fe_values(n_q);

    if (T::point_set_number == 1) {
      xq_hat = integration_proxy.get_quadrature_points(0);
      fe_values.set_points(xq_hat);
    }

    const std::size_t n_test_dof(test_fe_type::n_dof_per_element);

    // loop over the elements
    for (std::size_t k(k_begin); k < k_end; ++k) {

      // prepare the quadrature points if necessary
      if (T::point_set_number > 1) {
	xq_hat = integration_proxy.get_quadrature_points(k);
	fe_values.set_points(xq_hat);
      }

      if (form_type::require_space_coordinates)
	cell_type::map_points_to_space_coordinates(xq, m.get_vertices(),
						   m.get_cells(),
						   k, xq_hat);

      // prepare the basis function values if necessary
      if (form_type::differential_order == 1ul) {
	const array<double>& jmt(m.get_jmt(k));
	fe_values.prepare(jmt);
      }

      const array<double>& psi(fe_values.template get_values<test_fe_index>());

      // evaluate the weak form
      const double volume(m.get_cell_volume(k));
      for (std::size_t q(0); q < n_q; ++q) {
	integration_proxy.f.prepare(k, &xq.at(q, 0ul), &xq_hat.at(q, 0ul));

	for (std::size_t j(0); j < n_test_dof; ++j) {
	  rhs_el.at(k - k_begin, j) += omega.at(q) * volume *
	    integration_proxy.f(k, &xq.at(q, 0ul),
				&xq_hat.at(q, 0ul),
				&psi.at(q, j, 0ul));
	}
      }

    }
  }

  double& algebraic_equation_value(std::size_t a_dof) { return constraint_values[a_dof]; }

  expression<form<0,1,0> > get_test_function() const { return form<0,1,0>(); }

  const array<double>& get_coefficients() const { return f; }
  const std::vector<double>& get_constraint_values() const { return constraint_values; }

  void clear() {
    std::fill(f.get_data(),
	      f.get_data() + f.get_size(0),
	      0.0);

    constraint_values.clear();
  }

  void show(std::ostream& stream) const {
    stream << "rhs = [" << f.at(0);
    for (std::size_t j(1); j < f.get_size(0); ++j)
      stream << "; " << f.at(j);
    stream << "];" << std::endl;
  }

  void export_data(std::ostream& stream) const {
    auto p(stream.precision(16));

    for (std::size_t j(0); j < f.get_size(0); ++j)
      stream << f.at(j) << std::endl;

    stream.precision(p);
  }

private:
  thread_pool tp;

  const test_fes_type& test_fes;

  array<double> f;
  std::vector<double> constraint_values;
};

#endif /* _LINEAR_FORM_H_ */
