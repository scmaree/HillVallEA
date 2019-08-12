#pragma once

/*

Implementation by S.C. Maree, 2018
s.c.maree[at]amc.uva.nl
smaree.com

*/

#include "HillVallEA/fitness.h"
#include "CEC2013_niching_benchmark/cec2013.h"

namespace hillvallea
{

  class cec2013_t : public fitness_t
  {

  public:

    // data members
    CEC2013 *cec2013_function_pointer;
    int problem_index;
    
    cec2013_t(int problem_index)
    {
      this->problem_index = problem_index;
      
      // CEC2013 Niching benchmark setup
      // This is a wrapper around the global pointer cec2013_function_pointer
      //----------------------------------------------------------------
      cec2013_function_pointer = new CEC2013(problem_index);

      number_of_parameters = cec2013_function_pointer->get_dimension();
      maximum_number_of_evaluations = cec2013_function_pointer->get_maxfes();

    }
    cec2013_t() {}

    // any positive value
    void set_number_of_parameters(size_t & number_of_parameters)
    {
      this->number_of_parameters = number_of_parameters;
    }


    void get_param_bounds(vec_t & lower, vec_t & upper) const
    {

      lower.clear();
      lower.resize(number_of_parameters, 0);
      
      upper.clear();
      upper.resize(number_of_parameters, 5);
      
      for (int i = 0; i < number_of_parameters; ++i) {
        lower[i] = cec2013_function_pointer->get_lbound(i);
        upper[i] = cec2013_function_pointer->get_ubound(i);
      }

    }

    void define_problem_evaluation(solution_t & sol)
    {
    
      assert(cec2013_function_pointer->get_dimension() == (int) sol.param.size());
      
      // HillVallEA performs minimization!
      sol.f = -cec2013_function_pointer->evaluate(sol.param);
      sol.penalty = 0.0;

    }


    std::string name() const
    {
      std::stringstream ss;
      ss << "CEC2013_p" << problem_index;
      return ss.str();
    }


  };
}
