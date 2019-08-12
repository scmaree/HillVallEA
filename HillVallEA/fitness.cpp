
#include "fitness.h"
#include "mathfunctions.hpp"



hillvallea::fitness_t::fitness_t()
{
  number_of_evaluations = 0;
  number_of_parameters = 0;
  redefine_random_initialization = false;
}

hillvallea::fitness_t::~fitness_t() {}


void hillvallea::fitness_t::set_number_of_parameters(size_t & number_of_parameters)
{
  std::cout << "fitness_function error 'set_number_of_parameters' not implemented" << std::endl;
  assert(false);
  return;
}

void hillvallea::fitness_t::get_param_bounds(vec_t & lower, vec_t & upper) const
{
  std::cout << "fitness_function error 'get_param_bounds' not implemented" << std::endl;
  assert(false);
  return;
}

size_t hillvallea::fitness_t::get_number_of_parameters() const
{
  return number_of_parameters;
}

void hillvallea::fitness_t::evaluate(solution_t & sol)
{
  assert(sol.param.size() == number_of_parameters);

  define_problem_evaluation(sol);
  
  number_of_evaluations++;
}

void hillvallea::fitness_t::evaluate(solution_pt & sol)
{ 
  evaluate(*sol); 
}

// evaluates the function
// for new functions, set problem_evaluation.
// evaluate covers the evaluation itself and can be set to cover other stuff
// such as counting the number of evaluations or printing

void hillvallea::fitness_t::define_problem_evaluation(solution_t & sol)
{
  std::cout << "fitness_function error 'problem_evaluation' not implemented" << std::endl;
  assert(false);
  return;
}


std::string hillvallea::fitness_t::name() const
{
  std::cout << "fitness_function warning 'name' not implemented" << std::endl;
  return "no name";
}

void hillvallea::fitness_t::init_solutions_randomly(population_pt & population, size_t sample_size, size_t number_of_elites, rng_pt rng)
{
  std::cout << "fitness_function warning 'init_solutions_randomly' not implemented" << std::endl;
  redefine_random_initialization = false;
}
