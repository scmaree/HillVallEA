/*

HillVallEA

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA


*/

#include "solution.hpp"
#include "optimizer.hpp"
#include "amalgam.hpp"
#include "amalgam_univariate.hpp"
#include "iamalgam.hpp"
#include "iamalgam_univariate.hpp"
#include "cmsaes.hpp"

hillvallea::optimizer_pt hillvallea::init_optimizer(const int local_optimizer_index, const size_t number_of_parameters, const vec_t & lower_param_bounds, const vec_t & upper_param_bounds, double init_univariate_bandwidth, fitness_pt fitness_function, rng_pt rng)
{

  // parse settings
  switch (local_optimizer_index)
  {
    case 0: return std::make_shared<amalgam_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng); break;
    case 1: return std::make_shared<amalgam_univariate_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng); break;
    case 10: return std::make_shared<cmsaes_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng); break;
    case 20: return std::make_shared<iamalgam_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng); break;
    case 21: return std::make_shared<iamalgam_univariate_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng); break;
    default: return std::make_shared<amalgam_t>(number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng); break;
  }

}

// initialization of the general parameters of the EDAs
hillvallea::optimizer_t::optimizer_t(const size_t number_of_parameters, const vec_t & lower_param_bounds, const vec_t & upper_param_bounds, double init_univariate_bandwidth, fitness_pt fitness_function, rng_pt rng)
{

  active = true;
  this->number_of_parameters = number_of_parameters;
  this->lower_param_bounds = lower_param_bounds;
  this->upper_param_bounds = upper_param_bounds;
  this->fitness_function = fitness_function;
  number_of_generations = 0;
  this->rng = rng;
  pop = std::make_shared<population_t>();
  // best;
  average_fitness_history.resize(0); 
  selection_fraction = 0; // this will definitely cause weird stuff.
  this->init_univariate_bandwidth = init_univariate_bandwidth;
  maximum_no_improvement_stretch = 1000000;
  param_std_tolerance = 0;
  fitness_std_tolerance = 0;

}

hillvallea::optimizer_t::~optimizer_t() {}

 // Initialization
 //---------------------------------------------------------------------------------
void hillvallea::optimizer_t::initialize_from_population(population_pt pop)
{
  std::cout << "initialize_from_population not implemented" << std::endl;
  assert(false);
  return;
}


size_t hillvallea::optimizer_t::recommended_popsize(const size_t problem_dimension) const
{
  std::cout << "recommended_popsize not implemented" << std::endl;
  assert(false);
  return 0;
}

void hillvallea::optimizer_t::generation()
{
  std::cout << "generation not implemented" << std::endl;
  assert(false);
  return;
}

bool hillvallea::optimizer_t::checkTerminationCondition()
{
  std::cout << "checkTerminiationCriteria not implemented" << std::endl;
  assert(false);
  return true;
}

size_t hillvallea::optimizer_t::sample_new_population(const size_t sample_size)
{
  std::cout << "sample_new_population not implemented" << std::endl;
  assert(false);
  return 0;
}

std::string hillvallea::optimizer_t::name() const
{
  std::cout << "name not implemented" << std::endl;
  assert(false);
  return "name not implemented"; 
}

void hillvallea::optimizer_t::estimate_sample_parameters()
{
  std::cout << "estimate_sample_parameters not implemented" << std::endl;
  assert(false);
}

