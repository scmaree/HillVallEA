/*

HillVallEA 

Real-valued Multi-Modal Evolutionary Optimization

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA

Main script to run HillVallEA on the test problems of the CEC2013 niching bechmark

*/
// core_search_alg[a], cluster_alg[c],


#include "HillVallEA/hillvallea.hpp"
#include "cec2013_fitness_function.h"
void  write_CEC2013_niching_file(int core_search_alg, int cluster_alg, int problem, int run, std::vector<hillvallea::solution_pt> elitist_archive)
{
	
	std::ofstream file;
    std::string filename;
    std::string write_directory = "../niching_competition_data/";
    
	std::stringstream ss;
	ss << write_directory << "HillVallEA_core" << core_search_alg << "_cluster" << cluster_alg << "/problem" << std::setw(3) << std::setfill('0') << problem << "run" << std::setw(3) << std::setfill('0') << run << ".dat";
	filename = ss.str();
    
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    assert(file.is_open());

    int precision = std::numeric_limits<long double>::digits10 + 1;
    
    for (size_t i = 0; i < elitist_archive.size(); ++i) {
      file
      << std::fixed << std::setw(precision + 5) << std::setprecision(precision + 1) << elitist_archive[i]->param << " = "
      << std::fixed << std::setw(precision + 5) << std::setprecision(precision + 1) << elitist_archive[i]->f << " @ "
      << std::fixed << std::setw(precision + 5) << elitist_archive[i]->feval_obtained << " " 
      << std::fixed << std::setw(precision + 5) << std::setprecision(precision + 1) << elitist_archive[i]->time_obtained;

      file << std::endl;
    }
	
}

void run_CEC2013_niching_problem(int core_search_alg, int cluster_alg, int problem_index, int number_of_runs, double & mean_pr, double & mean_f1)
{

  // Run HillVallEA 
  //---------------------------------------------------------------
  hillvallea::vec_t pr(number_of_runs, 0.0);
  hillvallea::vec_t f1(number_of_runs, 0);
  
  #pragma omp parallel for
  for (int run = 0; run < number_of_runs; ++run)
  {
    
    // init fitness function
    //------------------------------------------------------
    hillvallea::fitness_pt fitness_function = std::make_shared<hillvallea::cec2013_t>(problem_index);
    hillvallea::vec_t lower_range_bounds, upper_range_bounds;
    fitness_function->get_param_bounds(lower_range_bounds, upper_range_bounds);
    

    // Run HillVallEA
    // Fix random seed for reproducibility
    //-------------------------------------------------------
    int random_seed = 23000 + run;

    hillvallea::hillvallea_t hillvallea(fitness_function, (int) fitness_function->number_of_parameters, lower_range_bounds, upper_range_bounds, fitness_function->maximum_number_of_evaluations, random_seed);
    hillvallea.local_optimizer_index = core_search_alg;
    hillvallea.cluster_alg = cluster_alg;
    
    hillvallea.run();

    // write the elitist archive
    // write_CEC2013_niching_file(core_search_alg, cluster_alg, problem_index, run + 1, hillvallea.elitist_archive);

    // Compute Peak Ratio using the CEC2013 guidelines
    //-------------------------------------------------------
    std::vector<std::vector<double>> optimum_candidates;
    for (size_t i = 0; i < hillvallea.elitist_archive.size(); ++i) {
      optimum_candidates.push_back(hillvallea.elitist_archive[i]->param);
    }

    
    // Compute the peak ratio & success rate
    std::vector<std::vector<double>> seeds;
    CEC2013 *cec2013_function_pointer = new CEC2013(problem_index);
    int peaks_found = how_many_goptima(optimum_candidates, seeds, cec2013_function_pointer, 1e-5, cec2013_function_pointer->get_rho());

    pr[run] = ((double)peaks_found) / cec2013_function_pointer->get_no_goptima();
    f1[run] = ((double)peaks_found) / optimum_candidates.size();

    std::cout << std::fixed << std::setw(7) << std::setprecision(0) << problem_index
      << std::fixed << std::setw(7) << std::setprecision(0) << run
      << std::fixed << std::setw(14) << std::setprecision(3) << pr[run]
      << std::fixed << std::setw(8) << std::setprecision(3) << f1[run]
      << std::endl;

    delete cec2013_function_pointer;
  }

  mean_pr = pr.mean();
  mean_f1 = f1.mean();

  

}

// Main: Run the CEC2013 niching benchmark
//--------------------------------------------------------
int main(int argc, char **argv)
{

  // experiment settings
  //---------------------------------------------------------------------------------------
  int number_of_runs = 50;

  std::cout << "Running HillVallEA on the problems of the CEC2013 niching benchmark" << std::endl;
  
  std::cout << "Problem    Run    Peak Ratio      SR" << std::endl;
  std::cout << "-------------------------------------" << std::endl;
  double f1;

  std::vector<int> problem_list;
  for (int i = 1; i <= 20; ++i) {
	  problem_list.push_back(i);
  }
  
  std::vector<int> core_search_alg;
  // core_search_alg.push_back(0);
  core_search_alg.push_back(1); // amalgam-univariate (AMu)
  // core_search_alg.push_back(10);
  // core_search_alg.push_back(20);
  // core_search_alg.push_back(21);
  
  std::vector<int> cluster_alg;
  cluster_alg.push_back(0); // HVC
  // cluster_alg.push_back(1); // HGML
  // cluster_alg.push_back(2); // NBC
  
  // iterate over the problem_list
  
  hillvallea::vec_t prs(problem_list.size(),1.0);
  
  for(size_t a = 0; a < core_search_alg.size(); ++a)
  {
  
    for(size_t c = 0; c < cluster_alg.size(); ++c)
    {
      std::cout << "Core Search alg = " << core_search_alg[a] << " with cluster alg = " << cluster_alg[c] << std::endl;
      
      for (int i = 0; i < problem_list.size(); ++i)
      {

        int problem_index = problem_list[i];
        
        run_CEC2013_niching_problem(core_search_alg[a], cluster_alg[c], problem_index, number_of_runs, prs[i], f1);

        std::cout << "-------------------------------------" << std::endl;
        std::cout << std::fixed << std::setw(7) << std::setprecision(0) << problem_index
          << "    avg"
          << std::fixed << std::setw(14) << std::setprecision(3) << prs[i]
          << std::fixed << std::setw(8) << std::setprecision(3) << f1
          << std::endl;
        std::cout << "-------------------------------------" << std::endl << std::endl;
      }

      std::cout << "SUM PR = " << prs.sum() << std::endl;
    }
  }
  return(0);
}
