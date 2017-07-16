

// So far this syntax is completely hypothetical. Will not compile.

double f_1(const double* x) { return x[1] > 0.5 ? 1.0 : 0.0; }
double f_2(const double* x) { return 0.0; }

int main(int argc, char *argv[]) {
  typedef cell::triangle cell_type;
  typedef finite_element::triangle_lagrange_p1_bubble velocity_fe_type;
  typedef finite_element::triangle_lagrange_p1 pressure_fe_type;
  typedef finite_element::composite_finite_element<velocity_fe_type,
						   velocity_fe_type,
						   pressure_fe_type>
    fe_type;
  typedef finite_element_space<fe_type> fes_type;
  typedef typename finite_element_space<fe_type>::element element_type;
  typedef quad::triangle::qf5pT q_type;

  const std::size_t n(100);
  mesh<cell_type> m(gen_square_mesh(1.0, 1.0, n, n));
  submesh<cell_type> dm(m.get_boundary_submesh());

  fes_type fes(m);
  fes.set_dirichlet_boundary<0>(dm);
  fes.set_dirichlet_boundary<1>(dm);

  double mu(1.0);
  
  bilinear_form<fes_type, fes_type> a(fes, fes); {
    const auto
      u_0(a.get_trial_function<0>()),
      u_1(a.get_trial_function<1>()),
      p(a.get_trial_function<2>());

    const auto
      v_0(a.get_trial_function<0>()),
      v_1(a.get_trial_function<1>()),
      q(a.get_trial_function<2>());
    
    a += integrate<q_type>(mu * (d<1>(u_0) * d<1>(v_0) +
				 d<2>(u_0) * d<2>(v_0) +
				 d<1>(u_1) * d<1>(v_1) +
				 d<2>(u_1) * d<2>(v_1))
			   - p * (d<1>(v_0) + d<2>(v_1))
			   - q * (d<1>(u_0) + d<2>(u_1)), m);

    linear_form<fes_type> f(fes); {
      const auto
	v_0(f.get_trial_function<0>()),
	v_1(f.get_trial_function<1>()),
	q(f.get_trial_function<2>());

      f += integrate<q_type>(f_1 * v_0 + f_2 * v_1, m);
    }

    const auto solution(a.solve(f));

    exporter::ensight6("stokes", solution, "up");
  }
    
  return 0;
}
