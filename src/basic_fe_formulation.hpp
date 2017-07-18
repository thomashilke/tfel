#ifndef BASIC_FE_FORMULATION_H
#define BASIC_FE_FORMULATION_H

#include "mesh.hpp"
#include "fe.hpp"
#include "fes.hpp"
#include "quadrature.hpp"
#include "projector.hpp"
#include "export.hpp"
#include "timer.hpp"
#include "expression.hpp"

#include "composite_fe.hpp"
#include "composite_fes.hpp"
#include "composite_form.hpp"


template<typename fe>
class basic_fe_formulation {
public:
  typedef fe fe_type;
  typedef typename fe_type::cell_type cell_type;
  typedef finite_element_space<fe_type> fes_type;
  typedef typename fes_type::element element_type;
  typedef typename default_quadrature<cell_type>::type volume_quadrature_type;
  typedef typename default_quadrature<typename cell_type::boundary_cell_type>::type boundary_quadrature_type;

  virtual ~basic_fe_formulation() {}
  
private:
};

template<typename ... fe_pack>
class basic_fe_formulation<composite_finite_element<fe_pack...> > {
public:
  using fe_type = composite_finite_element<fe_pack...>;
  using fe_list = type_list<fe_pack...>;
  typedef typename fe_type::cell_type cell_type;
  typedef composite_finite_element_space<fe_type> fes_type;
  typedef typename fes_type::element element_type;
  typedef typename default_quadrature<cell_type>::type volume_quadrature_type;
  typedef typename default_quadrature<typename cell_type::boundary_cell_type>::type boundary_quadrature_type;

  virtual ~basic_fe_formulation() {}
  
private:
};

#endif /* BASIC_FE_FORMULATION_H */
