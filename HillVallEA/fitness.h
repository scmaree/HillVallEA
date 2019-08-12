#pragma once

/*

hillvallea Multi-objective

By S.C. Maree, 2016
s.c.maree[at]amc.uva.nl
smaree.com

*/


#include <functional>
#include "hillvallea_internal.hpp"
#include "population.hpp"


// Defines the fitness function of our choice
// implementation idea from libcmaes by Emmanuel Benazera.
//------------------------------------------------

namespace hillvallea
{

  class fitness_t
  {

  public:

    fitness_t();
    ~fitness_t();

    size_t number_of_parameters;
    unsigned int number_of_evaluations;
    unsigned int maximum_number_of_evaluations;

    size_t get_number_of_parameters() const;

    // evaluates the function
    // for new functions, define problem_evaluation in "define_problem_evaluation".
    // evaluate covers the evaluation itself and can be set to cover other stuff
    // such as counting the number of evaluations or printing
    void evaluate(solution_t & sol);
    void evaluate(solution_pt & sol);

    // Placeholders for user-defined objective functions
    //----------------------------------------------------------------------------------------
    virtual void set_number_of_parameters(size_t & number_of_parameters);
    virtual void get_param_bounds(vec_t & lower, vec_t & upper) const;
    virtual void define_problem_evaluation(solution_t & sol);

    virtual std::string name() const;
    
    // redefine initialization
    bool redefine_random_initialization;
    
    virtual void init_solutions_randomly(population_pt & population, size_t sample_size, size_t number_of_elites, rng_pt rng);
    
  };


}
