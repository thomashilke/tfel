
#include "../src/core/projector.hpp"
#include "../src/utility/importer.hpp"
#include "../src/core/export.hpp"
#include "../src/core/composite_fe.hpp"
#include "../src/core/composite_fes.hpp"
#include "../src/core/fes.hpp"


double f(const double* x) {
  return std::sin(x[0] * M_PI);
}


void test_restriction_extension() {
  using cell_type = cell::edge;
  using mesh_type = mesh<cell_type>;
  using fe_type = cell::edge::fe::lagrange_p1_bubble;
  using fes_type = finite_element_space<fe_type>;
  
  mesh<cell_type> m(gen_segment_mesh(0.0, 1.0, 40));
  submesh<cell_type, cell_type> center_sm(m.query_elements([] (const double* x) -> bool {
    return x[0] > 0.25 and x[0] < 0.75;
  }));
  mesh<cell_type> center(center_sm);

  fes_type
    m_fes(m),
    center_fes(center);

  const fes_type::element f_h(projector::l2<fe_type, quad::edge::gauss3>(f, m_fes));
  exporter::ascii<mesh_type>("f_h.dat", to_mesh_vertex_data<fe_type>(f_h));

  const fes_type::element f_h_center(f_h.restrict(center_fes, center_sm));
  exporter::ascii<mesh_type>("f_h_center.dat", to_mesh_vertex_data<fe_type>(f_h_center));

  const fes_type::element f_h_extended(f_h_center.extend(m_fes, center_sm));
  exporter::ascii<mesh_type>("f_h_extended.dat", to_mesh_vertex_data<fe_type>(f_h_extended));
}


void test_alucell_import() {
  try {
    using cell_type = cell::tetrahedron;
    using mesh_type = mesh<cell_type>;
  
    mesh<cell_type> cuve(importer::alucell::mesh<cell_type>("/usr/scratch/master/alu-data/trunk/AP32/stat/dbfile_stat", "cuveb"));
    submesh<cell_type, cell_type> electrolyte(cuve.get_submesh_with_reference(2));
    mesh<cell_type> electrolyte_m(electrolyte);

    using fe_type = cell::tetrahedron::fe::lagrange_p1;
    using cfe_type = composite_finite_element<fe_type, fe_type, fe_type>;
    using cfes_type = composite_finite_element_space<cfe_type>;
  
    cfes_type velocity_fes(cuve);
    cfes_type electrolyte_velocity_fes(electrolyte_m);
  
    cfes_type::element
      velocity(to_composite_p1_finite_element_function<3, cell_type>(
		 importer::alucell::variable<cell_type>(
		   "/usr/scratch/master/alu-data/trunk/AP32/stat/dbfile_stat",
		   "cuveb", "vitesse", cuve)));

    cfes_type::element electrolyte_velocity(velocity.restrict(electrolyte_velocity_fes, electrolyte));
  
  
    exporter::ensight6(
      "cuveb",
      to_mesh_vertex_data<fe_type>(velocity.get_component<0>()), "v0",
      to_mesh_vertex_data<fe_type>(velocity.get_component<1>()), "v1",
      to_mesh_vertex_data<fe_type>(velocity.get_component<2>()), "v2");
  
    exporter::ensight6("electrolyte",
      to_mesh_vertex_data<fe_type>(electrolyte_velocity.get_component<0>()), "v0",
      to_mesh_vertex_data<fe_type>(electrolyte_velocity.get_component<1>()), "v1",
      to_mesh_vertex_data<fe_type>(electrolyte_velocity.get_component<2>()), "v2");

    /*exporter::ensight6("electrolyte",
		       velocity.get_component<0>().restrict(electrolyte_velocity_fes.get_finite_element_space<0>(), electrolyte), "v0",
		       velocity.get_component<1>().restrict(electrolyte_velocity_fes.get_finite_element_space<1>(), electrolyte), "v1",
		       velocity.get_component<2>().restrict(electrolyte_velocity_fes.get_finite_element_space<2>(), electrolyte), "v2");*/
  
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
}


int main(int argc, char *argv[]) {
  test_restriction_extension();
  test_alucell_import();

  return 0;
}
