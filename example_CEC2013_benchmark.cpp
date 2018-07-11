/*

HillVallEA 

Real-valued Multi-Modal Evolutionary Optimization

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA

Main script to run HillVallEA on the test problems of the CEC2013 niching bechmark

*/

#include "HillVallEA/hillvallea.hpp"
#include "CEC2013_niching_benchmark/cec2013.h"

CEC2013 *cec2013_function_pointer;


hillvallea::fitness_t cec2013 = [](hillvallea::solution_t & sol)
{

  assert(cec2013_function_pointer->get_dimension() == (int) sol.param.size());

  // HillVallEA performs minimization!
  sol.f = -cec2013_function_pointer->evaluate(sol.param);
  sol.penalty = 0.0;

};

void run_CEC2013_niching_problem(int problem_index, int number_of_runs, double & mean_pr, double & mean_f1)
{

  // CEC2013 Niching benchmark setup
  // This is a wrapper around the global pointer cec2013_function_pointer
  //----------------------------------------------------------------
  cec2013_function_pointer = new CEC2013(problem_index); 
  hillvallea::fitness_t * fitness_function = &cec2013;
  int number_of_parameters = cec2013_function_pointer->get_dimension();
  int maximum_number_of_evaluations = cec2013_function_pointer->get_maxfes();

  hillvallea::vec_t lower_range_bounds(number_of_parameters);
  hillvallea::vec_t upper_range_bounds(number_of_parameters);

  for (int i = 0; i < number_of_parameters; ++i) {
    lower_range_bounds[i] = cec2013_function_pointer->get_lbound(i);
    upper_range_bounds[i] = cec2013_function_pointer->get_ubound(i);
  }

  // Run HillVallEA 
  //---------------------------------------------------------------
  hillvallea::vec_t pr(number_of_runs, 0.0);
  hillvallea::vec_t f1(number_of_runs, 0);

  for (int run = 0; run < number_of_runs; ++run)
  {

    // Run HillVallEA
    // Fix random seed for reproducibility
    //-------------------------------------------------------
    int random_seed = 1000 + run;

    hillvallea::hillvallea_t hillvallea(fitness_function, number_of_parameters, lower_range_bounds, upper_range_bounds, maximum_number_of_evaluations, random_seed);
    hillvallea.run();

    // write the elitist archive
    // write_CEC2013_niching_file(problem_index, run + 1, hillvallea.elitist_archive);

    // Compute Peak Ratio using the CEC2013 guidelines
    //-------------------------------------------------------
    std::vector<std::vector<double>> optimum_candidates;
    for (size_t i = 0; i < hillvallea.elitist_archive.size(); ++i) {
      optimum_candidates.push_back(hillvallea.elitist_archive[i]->param);
    }

    // Compute the peak ratio & success rate
    std::vector<std::vector<double>> seeds;
    int peaks_found = how_many_goptima(optimum_candidates, seeds, cec2013_function_pointer, 1e-5, cec2013_function_pointer->get_rho());

    pr[run] = ((double)peaks_found) / cec2013_function_pointer->get_no_goptima();
    f1[run] = ((double)peaks_found) / optimum_candidates.size();

    std::cout << std::fixed << std::setw(7) << std::setprecision(0) << problem_index
      << std::fixed << std::setw(7) << std::setprecision(0) << run
      << std::fixed << std::setw(14) << std::setprecision(3) << pr[run]
      << std::fixed << std::setw(8) << std::setprecision(3) << f1[run]
      << std::endl;

  }

  mean_pr = pr.mean();
  mean_f1 = f1.mean();

  delete cec2013_function_pointer;

}

// Main: Run the CEC2013 niching benchmark
//--------------------------------------------------------
int main(int argc, char **argv)
{

  // experiment settings
  //---------------------------------------------------------------------------------------
  int number_of_runs = 5;

  std::cout << "Running HillVallEA on the problems of the CEC2013 niching benchmark" << std::endl;
  
  std::cout << "Problem    Run    Peak Ratio      F1" << std::endl;
  std::cout << "-------------------------------------" << std::endl;
  double pr, f1;

  // iterate over the problem_list
  for (int problem_index = 1; problem_index <= 20; ++problem_index)
  {

    run_CEC2013_niching_problem(problem_index, number_of_runs, pr, f1);

    std::cout << "-------------------------------------" << std::endl;
    std::cout << std::fixed << std::setw(7) << std::setprecision(0) << problem_index
      << "    avg"
      << std::fixed << std::setw(14) << std::setprecision(3) << pr
      << std::fixed << std::setw(8) << std::setprecision(3) << f1
      << std::endl;
    std::cout << "-------------------------------------" << std::endl << std::endl;
  }

  return(0);
}
