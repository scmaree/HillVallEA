/*

HillVallEA 

Real-valued Multi-Modal Evolutionary Optimization

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA

Example script to demonstrate the usage of HillVallEA
 on the well-known 2D Six Hump Camel Back function

*/

#include "HillVallEA/hillvallea.hpp"


// implementation of the Six Hump Camel Back function
hillvallea::fitness_t sixHumpCamelBack = [](hillvallea::solution_t & sol)
{

  // sanity check to see if the objective function is correctly called.
  assert(sol.param.size() == 2);
  
  // fitness function
  double p0s = sol.param[0]*sol.param[0]; // param 0 squared
  double p1s = sol.param[1]*sol.param[1]; // param 1 squared
  
  sol.f = (4.0-2.1*p0s + p0s*p0s/3.0) * p0s + sol.param[0]*sol.param[1] + (-4.0 + 4.0*p1s)*p1s;
  
  // penalty function, if set to > 0, solution is infeasible
  // (not extensively tested)
  sol.penalty = 0.0;

};


// Main: Run the CEC2013 niching benchmark
//--------------------------------------------------------
int main(int argc, char **argv)
{
  
  // Problem definition
  // Note: define as minimization problem!
  //-----------------------------------------
  size_t number_of_parameters = 2;                               // number of problem variables
  hillvallea::vec_t lower_range_bounds(number_of_parameters);   // lower variables bounds
  hillvallea::vec_t upper_range_bounds(number_of_parameters);   // upper variables bounds
  lower_range_bounds[0] = -3.0;
  lower_range_bounds[1] = -2.0;
  upper_range_bounds[0] = 3.0;
  upper_range_bounds[1] = 2.0;
  
  // initialize HillVallEA uniform on the entire search space
  // decrease if optimum is expected to be in a subspace of the entire search space
  hillvallea::vec_t lower_init_ranges = lower_range_bounds;
  hillvallea::vec_t upper_init_ranges = upper_range_bounds;
  
  // HillVallEA Settings
  //-----------------------------------------
  // Type of local optimizer to be used.
  // 0 = AMaLGaM, 1 = AMaLGaM-Univariate, 20 = iAMaLGaM, 21 = iAMaLGaM-Univariate
  size_t local_optimizer_index = 1; // AMaLGaM-Univariate (1) is suggested
  
  int maximum_number_of_evaluations = 10000; // maximum number of evaluations
  int maximum_number_of_seconds = 60; // maximum runtime in seconds
  
  // if the optimum is known, you can terminate HillVallEA if it found a solution
  // with fitness below the value_to_reach (vtr)
  double value_to_reach = 0;
  bool use_vtr = false;
  
  // random seed initialization for reproducibility
  int random_seed = 100;
  
  // Output to test files
  bool write_generational_solutions = false;
  bool write_generational_statistics = true;
  std::string write_directory = "./";
  std::string file_appendix = ""; // can be used when multiple runs are outputted in the same directory
  
  // Initialization of HillVallEA
  //-----------------------------------------
  hillvallea::hillvallea_t opt(
     &sixHumpCamelBack,
     local_optimizer_index,
     number_of_parameters,
     lower_init_ranges,
     upper_init_ranges,
     lower_range_bounds,
     upper_range_bounds,
     maximum_number_of_evaluations,
     maximum_number_of_seconds,
     value_to_reach,
     use_vtr,
     random_seed,
     write_generational_solutions,
     write_generational_statistics,
     write_directory,
     file_appendix
  );

  // Running HillVallEA
  std::cout << "Running HillVallEA on the Six Hump Camel back function" << std::endl;
  
  
  opt.run();
  
  std::cout << "HillVallEA finished" << std::endl;
  std::cout << "Generation statistics written to " << write_directory << "statistics" << file_appendix << ".dat" << std::endl;
  std::cout << "Elitist archive written to       " << write_directory << "elites" << file_appendix << ".dat" << std::endl;
  
  std::cout << "HillVallEA Obtained " << opt.elitist_archive.size() << " elites: " << std::endl;
  
  std::cout << "    Fitness      Penalty   Params" << std::endl;
  for(size_t i = 0; i < opt.elitist_archive.size(); ++i)
  {
    std::cout << std::setw(11) << std::scientific << std::setprecision(3) << opt.elitist_archive[i]->f << "  ";
    std::cout << std::setw(11) << std::scientific << std::setprecision(3) << opt.elitist_archive[i]->penalty << "  ";
    std::cout << std::setw(11) << std::scientific << std::setprecision(3) << opt.elitist_archive[i]->param << std::endl;
  }
  

  return(0);
}
