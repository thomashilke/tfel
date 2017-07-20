#ifndef _COMPOSITE_FORM_H_
#define _COMPOSITE_FORM_H_

#include "composite_fes.hpp"
#include "form.hpp"


template<typename fe_list, std::size_t n, std::size_t m, typename unique_fe_list>
struct select_function_valuation_impl {
  static void call(const double** ptr, std::size_t q, std::size_t i,
		   const fe_value_manager<unique_fe_list>& values,
		   const fe_value_manager<unique_fe_list>& zvalues) {
    const std::size_t fe_index = get_index_of_element<head_t<fe_list>, unique_fe_list>::value;
    
    *ptr = &(zvalues.template get_values<fe_index>().at(q,0,0));
    select_function_valuation_impl<tail_t<fe_list>, n + 1, m, unique_fe_list>::call(ptr + 1, q, i, values, zvalues);
  }
};

template<typename fe_list, std::size_t m, typename unique_fe_list>
struct select_function_valuation_impl<fe_list, m, m, unique_fe_list> {
  static void call(const double** ptr, std::size_t q, std::size_t i,
		   const fe_value_manager<unique_fe_list>& values,
		   const fe_value_manager<unique_fe_list>& zvalues) {
    const std::size_t fe_index = get_index_of_element<head_t<fe_list>, unique_fe_list>::value;
    
    *ptr = &(values.template get_values<fe_index>().at(q,i,0));
    select_function_valuation_impl<tail_t<fe_list>, m + 1, m, unique_fe_list>::call(ptr + 1, q, i, values, zvalues);
  }
};

template<std::size_t n, std::size_t m, typename unique_fe_list>
struct select_function_valuation_impl<type_list<>, n, m, unique_fe_list> {
  static void call(const double** ptr, std::size_t q, std::size_t i,
		   const fe_value_manager<unique_fe_list>& values,
		   const fe_value_manager<unique_fe_list>& zvalues) {}
};



template<typename fe_list, std::size_t m, typename unique_fe_list>
void select_function_valuation(const double** ptr, std::size_t q, std::size_t i,
			       const fe_value_manager<unique_fe_list>& values,
			       const fe_value_manager<unique_fe_list>& zvalues) {
  select_function_valuation_impl<fe_list, 0, m, unique_fe_list>::call(ptr, q, i, values, zvalues);
}


#include "composite_linear_form.hpp"
#include "composite_bilinear_form.hpp"


#endif /* _COMPOSITE_FORM_H_ */
