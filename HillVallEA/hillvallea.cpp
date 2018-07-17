/*

HillVallEA

By S.C. Maree
s.c.maree[at]amc.uva.nl
github.com/SCMaree/HillVallEA


*/

#include "hillvallea.hpp"
#include "population.hpp"
#include "mathfunctions.hpp"

namespace hillvallea
{
  // Constructor
  hillvallea_t::hillvallea_t(
    fitness_t * fitness_function,
    const int local_optimizer_index,
    const int  number_of_parameters,
    const vec_t & lower_init_ranges,
    const vec_t & upper_init_ranges,
    const vec_t & lower_param_bounds,
    const vec_t & upper_param_bounds,
    const int maximum_number_of_evaluations,
    const int maximum_number_of_seconds,
    const double vtr,
    const bool use_vtr,
    const int random_seed,
    const bool write_generational_solutions,
    const bool write_generational_statistics,
    const std::string write_directory,
    const std::string file_appendix
  )
  {

    // copy all settings
    this->fitness_function = fitness_function;
    this->local_optimizer_index = local_optimizer_index;
    this->number_of_parameters = number_of_parameters;
    this->lower_init_ranges = lower_init_ranges;
    this->upper_init_ranges = upper_init_ranges;
    this->lower_param_bounds = lower_param_bounds;
    this->upper_param_bounds = upper_param_bounds;
    this->maximum_number_of_evaluations = maximum_number_of_evaluations;
    this->maximum_number_of_seconds = maximum_number_of_seconds;
    this->vtr = vtr;
    this->use_vtr = use_vtr;
    this->random_seed = random_seed;
    this->write_generational_solutions = write_generational_solutions;
    this->write_generational_statistics = write_generational_statistics;
    this->write_directory = write_directory;
    this->file_appendix = file_appendix;

    rng = std::make_shared<std::mt19937>((unsigned long)(random_seed));
    std::uniform_real_distribution<double> unif(0, 1);

    // Parameters of the recursion scheme
    //---------------------------------------------
    population_size_initializer = 6.0;
    population_size_incrementer = 2.0;
    cluster_size_initializer = 1.0;
    cluster_size_incrementer = 1.2;
    add_elites_max_trials = 5;

    search_volume = 1.0;
    for (int i = 0; i < number_of_parameters; ++i) {
      search_volume *= (upper_init_ranges[i] - lower_init_ranges[i]);
    }

    clustering_max_number_of_neighbours = (size_t)(number_of_parameters + 1);
    TargetTolFun = 1e-5;
  }


  // Quick Constructor
  hillvallea_t::hillvallea_t(
    fitness_t * fitness_function,
    const int  number_of_parameters,
    const vec_t & lower_param_bounds,
    const vec_t & upper_param_bounds,
    const int maximum_number_of_evaluations,
    const int random_seed
  )
  {

    // copy all settings
    this->fitness_function = fitness_function;
    this->local_optimizer_index = 1; // default: AMu
    this->number_of_parameters = number_of_parameters;
    this->lower_init_ranges = lower_param_bounds; // init range == bounds
    this->upper_init_ranges = upper_param_bounds;  // init range == bounds
    this->lower_param_bounds = lower_param_bounds;
    this->upper_param_bounds = upper_param_bounds;
    this->maximum_number_of_evaluations = maximum_number_of_evaluations;
    this->maximum_number_of_seconds = 0;
    this->vtr = 0;
    this->use_vtr = false;
    this->random_seed = random_seed;
    this->write_generational_solutions = false;
    this->write_generational_statistics = false;
    this->write_directory = "";
    this->file_appendix = "";

    rng = std::make_shared<std::mt19937>((unsigned long)(random_seed));
    std::uniform_real_distribution<double> unif(0, 1);

    // Parameters of the recursion scheme
    //---------------------------------------------
    population_size_initializer = 6.0;
    population_size_incrementer = 2.0;
    cluster_size_initializer = 1.0;
    cluster_size_incrementer = 1.2;
    add_elites_max_trials = 5;

    search_volume = 1.0;
    for (int i = 0; i < number_of_parameters; ++i) {
      search_volume *= (upper_init_ranges[i] - lower_init_ranges[i]);
    }

    clustering_max_number_of_neighbours = (size_t)(number_of_parameters + 1);
    TargetTolFun = 1e-5;
  }


  hillvallea_t::~hillvallea_t()
  {

  }

  // Write statistic Files
  //----------------------------------------------------------------------------------
  void hillvallea_t::new_statistics_file()
  {
    std::string filename;
    if (file_appendix.empty()) {
      filename = write_directory + "statistics.dat";
    }
    else {
      filename = write_directory + "statistics" + file_appendix + ".dat";
    }

    statistics_file.open(filename, std::ofstream::out | std::ofstream::trunc);
    assert(statistics_file.is_open());

    statistics_file << "Pop  Cluster    Gen   Evals        Time No.Elites  Best-elite   Average-obj       Std-obj" << std::endl;

  }

  void hillvallea_t::write_statistics_line_population(const population_t & pop, const std::vector<optimizer_pt> & local_optimizers, const std::vector<solution_pt> & elitist_archive)
  {
    
    
    clock_t current_time = clock();
    double runtime = double(current_time - starting_time) / CLOCKS_PER_SEC;
    
    solution_t best = *pop.first();
    
    for (auto sol = elitist_archive.begin(); sol != elitist_archive.end(); ++sol)
    {
      if (solution_t::better_solution(**sol, best))
      {
        best = **sol;
      }
    }
    
    statistics_file
      << std::setw(3) << number_of_generations
      << std::setw(9) << local_optimizers.size()
      << std::setw(7) << 0
      << std::setw(8) << number_of_evaluations
      << std::setw(12) << std::scientific << std::setprecision(3) << runtime
      << std::setw(10) << elitist_archive.size()
    << std::setw(12) << std::scientific << std::setprecision(3) <<  best.f
      << std::setw(14) << std::scientific << std::setprecision(3) << pop.average_fitness()
      << std::setw(14) << std::scientific << std::setprecision(3) << pop.relative_fitness_std()
      << std::endl;
  }
  
  void hillvallea_t::write_statistics_line_cluster(const population_t & cluster_pop, int cluster_number, int cluster_generation, const std::vector<optimizer_pt> & local_optimizers, const std::vector<solution_pt> & elitist_archive)
  {
    
    
    clock_t current_time = clock();
    double runtime = double(current_time - starting_time) / CLOCKS_PER_SEC;
    
    solution_t best = *cluster_pop.first();
    
    for (auto sol = elitist_archive.begin(); sol != elitist_archive.end(); ++sol)
    {
      if (solution_t::better_solution(**sol, best))
      {
        best = **sol;
      }
    }
    
    statistics_file
    << std::setw(3) << ""
    << std::setw(9) << cluster_number
    << std::setw(7) << cluster_generation
    << std::setw(8) << number_of_evaluations
    << std::setw(12) << std::scientific << std::setprecision(3) << runtime
    << std::setw(10) << elitist_archive.size()
    << std::setw(12) << std::scientific << std::setprecision(3) << best.f
    << std::setw(14) << std::scientific << std::setprecision(3) << cluster_pop.average_fitness()
    << std::setw(14) << std::scientific << std::setprecision(3) << cluster_pop.relative_fitness_std()
    << std::endl;
  }

  
  void hillvallea_t::close_statistics_file()
  {
    statistics_file.close();
  }


  // Write population
  //-------------------------------------------------------------------------------
  void hillvallea_t::write_population_file(population_pt pop, std::vector<optimizer_pt> & local_optimizers) const
  {
    std::ofstream file;
    std::stringstream ss;
    ss << write_directory << "population" << std::setw(5) << std::setfill('0') << number_of_generations << "_inital_population" << file_appendix << ".dat";
    file.open(ss.str(), std::ofstream::out | std::ofstream::trunc);
    assert(file.is_open());

    for (auto sol = pop->sols.begin(); sol != pop->sols.end(); ++sol)
    {
      file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << (*sol)->param << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << (*sol)->f << " " << (*sol)->penalty << std::endl;
    }

    file.close();
  }

  // Write selection
  //-------------------------------------------------------------------------------
  void hillvallea_t::write_selection_file(population_pt pop, std::vector<optimizer_pt> & local_optimizers) const
  {

    // print clusters to file so that i can check them in matlab.
    std::ofstream file;
    std::stringstream ss;
    ss << write_directory << "population" << std::setw(5) << std::setfill('0') << number_of_generations << "_selection" << file_appendix << ".dat";
    file.open(ss.str(), std::ofstream::out | std::ofstream::trunc);
    assert(file.is_open());

    for (auto sol = pop->sols.begin(); sol != pop->sols.end(); ++sol)
    {
      file << (*sol)->param << (*sol)->f << " " << (*sol)->penalty << std::endl;
    }

    file.close();
  }


  // Write clusters
  //-------------------------------------------------------------------------------
  void hillvallea_t::write_cluster_population(int generation_nuber, size_t cluster_number, int cluster_generation, population_pt pop) const
  {
    std::ofstream file;
    std::stringstream ss;
    ss << write_directory << "population" << std::setw(5) << std::setfill('0') << generation_nuber << "_cluster" << std::setw(5) << std::setfill('0') << cluster_number << "_generation" << cluster_generation << file_appendix << ".dat";
    file.open(ss.str(), std::ofstream::out | std::ofstream::trunc);
    assert(file.is_open());

    for (auto sol = pop->sols.begin(); sol != pop->sols.end(); ++sol)
    {
      file << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << (*sol)->param << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << (*sol)->f << " " << (*sol)->penalty << std::endl;
    }

    file.close();

  }
  
  void hillvallea_t::write_elitist_archive_file(const std::vector<solution_pt> & elitist_archive, bool final) const
  {
    std::ofstream file;
    std::string filename;
    
    if(final) {
      filename = write_directory + "elites" + file_appendix + ".dat";
    }
    else
    {
      std::stringstream ss;
      ss << write_directory << "population" << std::setw(5) << std::setfill('0') << number_of_generations << "_elites" << file_appendix << ".dat";
      filename = ss.str();
    }
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    assert(file.is_open());

    int prescision = std::numeric_limits<long double>::digits10 + 1;
    
    for (size_t i = 0; i < elitist_archive.size(); ++i) {
      file
      << std::fixed << std::setw(prescision + 5) << std::setprecision(prescision + 1) << elitist_archive[i]->param
      << std::fixed << std::setw(prescision + 5) << std::setprecision(prescision + 1) << elitist_archive[i]->f
      << std::fixed << std::setw(prescision + 5) << std::setprecision(prescision + 1) << elitist_archive[i]->penalty;

      if (i == 0 && use_vtr)
        file << " " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << vtr << " " << success;

      file << std::endl;
    }
  }

  // Termination Criteria
  //-------------------------------------------------------------------------------
  bool hillvallea_t::terminate_on_runtime() const
  {
    // stop if we run out of time.
    if (maximum_number_of_seconds > 0)
    {
      clock_t current_time = clock();
      double runtime = double(current_time - starting_time) / CLOCKS_PER_SEC;
      
      if (runtime > maximum_number_of_seconds) {
        return true;
      }
    }
    
    return false;
  }

  bool hillvallea_t::terminate_on_approaching_elite(optimizer_t & local_optimizer, std::vector<solution_pt> & elite_candidates)
  {
    // find the nearest (candidate) elite that has similar or better fitness
    solution_pt nearest_elite;
    double distance_to_nearest_elite = 1e300;
    double distance;
    double best_fitness_so_far = local_optimizer.pop->sols[0]->f;
    double TargetTolFun = 1e-5;
    if (elitist_archive.size() > 0)
    {
      // find the nearest elite
      for (size_t j = 0; j < elitist_archive.size(); ++j)
      {

        // only consider elites that have better fitness
        if (elitist_archive[j]->f < (local_optimizer.pop->sols[0]->f + TargetTolFun))
        {

          distance = elitist_archive[j]->param_distance(*local_optimizer.pop->sols[0]);
          if (distance < distance_to_nearest_elite)
          {
            distance_to_nearest_elite = distance;
            nearest_elite = elitist_archive[j];
          }

          // also find the best elite in the archive for later
          if (elitist_archive[j]->f < best_fitness_so_far) {
            best_fitness_so_far = elitist_archive[j]->f;
          }
        }

      }
    }

    for (size_t j = 0; j < elite_candidates.size(); ++j)
    {
      // only consider good elite candidates
      if (elite_candidates[j]->f < (best_fitness_so_far + TargetTolFun))
      {
        distance = elite_candidates[j]->param_distance(*local_optimizer.pop->sols[0]);

        if (distance < distance_to_nearest_elite) {
          distance_to_nearest_elite = distance;
          nearest_elite = elite_candidates[j];
        }
      }
    }

    if (nearest_elite != nullptr)
    {
      if (check_edge(*nearest_elite, *local_optimizer.pop->sols[0], 5)) {
        return true;
      }
    }

    return false;
  }

  bool hillvallea_t::terminate_on_converging_to_local_optimum(optimizer_t & local_optimizer, std::vector<solution_pt> & elite_candidates)
  {

    // if this is the first, never terminate.
    if (elite_candidates.size() == 0 && elitist_archive.size() == 0) {
      return false;
    }

    // find best solution in the archive and candidates
    //----------------------------------------------------------
    double best = 1e308;
    double max_number_of_generations_to_obtain_elite = 0;

    for (size_t i = 0; i < elitist_archive.size(); ++i) {
      if (elitist_archive[i]->f < best) {
        best = elitist_archive[i]->f;
      }

      if (elitist_archive[i]->generation_obtained > max_number_of_generations_to_obtain_elite) {
        max_number_of_generations_to_obtain_elite = elitist_archive[i]->generation_obtained;
      }
    }

    for (size_t i = 0; i < elite_candidates.size(); ++i) {
      if (elite_candidates[i]->f < best) {
        best = elite_candidates[i]->f;
      }

      if (elite_candidates[i]->generation_obtained > max_number_of_generations_to_obtain_elite) {
        max_number_of_generations_to_obtain_elite = elite_candidates[i]->generation_obtained;
      }
    }
    
    // Only terminate if the local optimizer if its best is significantly worse than the best elite
    if (local_optimizer.pop->sols[0]->f > best + TargetTolFun)
    {

      // compute the time to optimum
      //----------------------------------------------------------
      int lookback_window = std::min((int)local_optimizer.average_fitness_history.size(), 5);

      if (lookback_window < 5) {
        return false;
      }

      // terminate only if the recent averages are all decreasing
      // note, if this passes, it implies curr_dto < prev_dto,
      // and thus tto > 0
      for (size_t i = 0; i < lookback_window-1; ++i)
      {
        if (local_optimizer.average_fitness_history[local_optimizer.average_fitness_history.size() - lookback_window + i] 
            < local_optimizer.average_fitness_history[local_optimizer.average_fitness_history.size() - lookback_window + i + 1]) {
          return false;
        }
      }

      double curr_fitness = local_optimizer.average_fitness_history.back();
      double prev_fitness = local_optimizer.average_fitness_history[local_optimizer.average_fitness_history.size() - lookback_window];

      double curr_dto = curr_fitness - best;
      double prev_dto = prev_fitness - best;

      // this is the best so far, let it run
      if (curr_dto <= 0.0) {
        return false;
      }

      // normal case     
      double dto_vtr = 1e-12;
      double cr5 = (curr_dto - prev_dto) / prev_dto;
      double cr = cr5; // pow((1 + cr5), 1.0 / lookback_window) - 1;

      double tto = lookback_window * log(dto_vtr / curr_dto) / log(1.0 + cr);
      
      if ((local_optimizer.number_of_generations + tto) > 50 * max_number_of_generations_to_obtain_elite) {
        return true;
      }

    }

    return false;

  }

  //----------------------------------------------------------------------------------------------
  // samples an initial population uniformly random, clusters it into a set of local_optimizers
  void hillvallea_t::initialize(population_pt pop, size_t population_size, std::vector<optimizer_pt> & local_optimizers, const std::vector<solution_pt> & elitist_archive)
  {

    // Initialize running parameters of hillvallea
    //-------------------------------------------------
    local_optimizers.clear();

    // initially, we create a single cluster that we initialize by uniform sampling
    //-------------------------------------------------------------------------------------------------------------
    std::vector<solution_pt> backup_sols = pop->sols;
    pop->fill_uniform(population_size, number_of_parameters, lower_init_ranges, upper_init_ranges, rng);

    {
      int fevals = pop->evaluate(this->fitness_function, 0); // no elite yet.
      number_of_evaluations += fevals;
      number_of_evaluations_init += fevals;
    }

    pop->sort_on_fitness();

    // create a dummy local_optimizer for the initial population so that we can perform selection and we can write it down.
    double init_univariate_bandwidth = 1.0;
    optimizer_pt local_optimizer = init_optimizer(local_optimizer_index, number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng);
    local_optimizer->initialize_from_population(pop);

    population_pt selection = std::make_shared<population_t>();
    local_optimizer->pop->truncation_percentage(*selection, local_optimizer->selection_fraction);

    // Hill-Valley Clustering
    // note: i init the univariate bandwidth base don the average edge length.. so this is correct here :)
    //------------------------------------------------
    std::vector<population_pt> clusters;

    clustering(*selection, clusters, elitist_archive, init_univariate_bandwidth);

    // Init local optimizers
    //---------------------------------------------------------------------------
    for (auto cluster = clusters.begin(); cluster != clusters.end(); ++cluster)
    {
      if ((*cluster)->sols.size() > 0)
      {
        optimizer_pt opt = init_optimizer(local_optimizer_index, number_of_parameters, lower_param_bounds, upper_param_bounds, init_univariate_bandwidth, fitness_function, rng);
        opt->initialize_from_population(*cluster);
        opt->average_fitness_history.push_back(opt->pop->average_fitness());
        local_optimizers.push_back(opt);
      }
    }

    if (write_generational_statistics) {
      write_statistics_line_population(*pop, local_optimizers, elitist_archive);
    }

    if (write_generational_solutions)
    {
      write_population_file(pop, local_optimizers);
      write_selection_file(selection, local_optimizers);
    }

  }

  void hillvallea_t::run()
  {

    //---------------------------------------------
    // reset all runlogs (in case hillvallea is run multiple time)
    starting_time = clock();
    success = false;
    terminated = false;
    number_of_evaluations = 0;
    number_of_evaluations_init = 0;
    number_of_evaluations_clustering = 0;
    number_of_generations = 0;
    bool restart = true;
    int number_of_generations_without_new_clusters = 0;
    elitist_archive.clear();


    // allocate population
    //---------------------------------------------
    pop = std::make_shared<population_t>();

    // Init population sizes
    //---------------------------------------------
    double current_population_size = number_of_parameters * pow(2.0, population_size_initializer);
    double current_cluster_size;
    {
      hillvallea::optimizer_pt dummy_optimizer = init_optimizer(local_optimizer_index, number_of_parameters, lower_param_bounds, upper_param_bounds, 1.0, fitness_function, rng);
      current_cluster_size = cluster_size_initializer *dummy_optimizer->recommended_popsize(number_of_parameters);
    }

    if(write_generational_statistics) {
      new_statistics_file();
    }
    
    // The big restart scheme
    //--------------------------------------------
    while (restart)
    {

      // each restart, collect the elite_candidates and add them to the archive at the end of the run
      // directly adding them to the archive is more expensive in terms of fevals, as you might find a bunch of local optima first
      // and adding solutions to the elitist archive costs fevals because we check them using the HillVallyTest
      std::vector<solution_pt> elite_candidates;
      std::vector<optimizer_pt> local_optimizers;

      // stop if the init popsize is too large. 
      // there is a bit to gain here by decreasing the popsize so that we can run some more local opts, but I don't think its worth it. 
      if (maximum_number_of_evaluations > 0 && current_population_size > (maximum_number_of_evaluations - number_of_evaluations)) {
        restart = false;
        break;
      }

      // compute initial population
      initialize(pop, (size_t) current_population_size, local_optimizers, elitist_archive);
      
      // we only create local optimizers from the global opts
      // so the local optimizer still inits new global opts
      // therefore, this is basically never hit.
      if (local_optimizers.size() == 0) 
      {
        current_population_size *= population_size_incrementer*population_size_incrementer;
        number_of_generations_without_new_clusters++;
        
        if (number_of_generations_without_new_clusters == 3) {
          restart = false;
          break;
        }
        else {
          continue;
        }
      }
      else {
        number_of_generations_without_new_clusters = 0;
      }

      // Run each of the local optimizers until convergence
      for (size_t i = 0; i < local_optimizers.size(); ++i)
      {

        // current local optimizer
        while (true)
        {

          // stop if the feval budget is reached
          size_t fevals_needed_to_check_elites = 1 + (size_t)(elite_candidates.size() * add_elites_max_trials * (elitist_archive.size() + elite_candidates.size() * 0.5));

          if (maximum_number_of_evaluations > 0 && number_of_evaluations + fevals_needed_to_check_elites + current_cluster_size >= maximum_number_of_evaluations) {
            restart = false;
            if (local_optimizers[i]->pop->size() > 0) {
              elite_candidates.push_back(local_optimizers[i]->pop->sols[0]);
            }
            break;
          }

          // stop if we run out of time.
          if (terminate_on_runtime()) {
            restart = false;
            if (local_optimizers[i]->pop->size() > 0) {
              elite_candidates.push_back(local_optimizers[i]->pop->sols[0]);
            }
            break;
          }

          // stop if the vtr is hit
          if (use_vtr && local_optimizers[i]->pop->sols[0]->f < vtr)
          {
            restart = false;
            best = *local_optimizers[i]->pop->sols[0];
            success = true;
            break;
          }

          // stop this local optimizer if it approaches a previously obtained elite (candidate) 
          if ((1 + local_optimizers[i]->number_of_generations) % 5 == 0)
          {
            if (terminate_on_approaching_elite(*local_optimizers[i], elite_candidates)) {
              local_optimizers[i]->active = false;
              break;
            }
          }

          // stop this local optimizer if it converges to a local optimum
          if (terminate_on_converging_to_local_optimum(*local_optimizers[i], elite_candidates)) {
            local_optimizers[i]->active = false;
            break;
          }

          // if the cluster is active, and after checking it, it is terminated, 
          // we add the best solution to the elitist archive
          if (local_optimizers[i]->active && local_optimizers[i]->checkTerminationCondition()) 
          {
            if (local_optimizers[i]->pop->size() > 0)
            {
              if (elitist_archive.size() == 0 || local_optimizers[i]->pop->sols[0]->f < elitist_archive[0]->f + TargetTolFun) {
                elite_candidates.push_back(local_optimizers[i]->pop->sols[0]);
                elite_candidates.back()->generation_obtained = local_optimizers[i]->number_of_generations;
              }
              
              break;
            }
          }

          // if it is still active, run a generation of the local optimizer
          if (local_optimizers[i]->active)
          {

            local_optimizers[i]->estimate_sample_parameters();

            int local_number_of_evaluations = (int)local_optimizers[i]->sample_new_population((size_t) current_cluster_size);
            number_of_evaluations += local_number_of_evaluations;

            if (write_generational_solutions) {
              write_cluster_population(number_of_generations, i, local_optimizers[i]->number_of_generations, local_optimizers[i]->pop);
            }

            local_optimizers[i]->pop->truncation_percentage(*local_optimizers[i]->pop, local_optimizers[i]->selection_fraction);
            local_optimizers[i]->average_fitness_history.push_back(local_optimizers[i]->pop->average_fitness());

            //if (write_generational_solutions) {
            //  write_cluster_selection(number_of_generations, i, local_optimizers[i]->number_of_generations, local_optimizers[i]->pop);
            //}
            
            if (write_generational_statistics) {
              write_statistics_line_cluster(*local_optimizers[i]->pop, i, local_optimizers[i]->number_of_generations, local_optimizers, elitist_archive);
            }
          }
        }
      }


      // write elitist_archive of this generation.
      if (write_generational_solutions) {
        write_elitist_archive_file(elitist_archive, false);
      }

      // check if the elites are novel and add the to the archive. 
      int number_of_new_global_opts_found = -1;
      int number_of_global_opts_found = -1;
      add_elites_to_archive(elitist_archive, elite_candidates, number_of_global_opts_found, number_of_new_global_opts_found);

      // if we found no new global opt, this is either due to the fact that there are no new basins found, 
      // or cuz the cluser size is too small.  increase both
      if (number_of_new_global_opts_found == 0) {
        current_cluster_size *= cluster_size_incrementer;
        current_population_size *= population_size_incrementer;
      }
      number_of_generations++;

    }

    // sort the archive s.t. the best is first
    if(elitist_archive.size() > 0) {
      best = *elitist_archive[0]; 
    }

    // write the final solution(s)
    // only if we care to output anything. Else, write nothing 
    if (write_generational_statistics || write_generational_solutions) {
      write_elitist_archive_file(elitist_archive, true);
    }

    if (write_generational_statistics) {
      close_statistics_file();
    }


  }


  void hillvallea_t::add_elites_to_archive(std::vector<solution_pt> & elitist_archive, const std::vector<solution_pt> & elite_candidates, int & number_of_global_opts_found, int & number_of_new_global_opts_found)
  {

    number_of_global_opts_found = 0;
    number_of_new_global_opts_found = 0;

    // find best solution in the archive
    double best_archive = 1e308;
    for (size_t i = 0; i < elitist_archive.size(); ++i)
    {
      if (elitist_archive[i]->f < best_archive) {
        best_archive = elitist_archive[i]->f;
      }
    }

    // find best candidate
    double best_candidate = 1e308;
    for (size_t i = 0; i < elite_candidates.size(); ++i)
    {
      if (elite_candidates[i]->f < best_candidate) {
        best_candidate = elite_candidates[i]->f;
      }
    }

    // clear the achive if the best candidate is better
    // than the best in the archive
    double best = std::min(best_candidate, best_archive);
    if ((best_candidate + TargetTolFun) < best_archive) {
      elitist_archive.clear();
    }

    // potential candidates are only those that are global optima 
    std::vector<solution_pt> potential_candidates;
    for (size_t i = 0; i < elite_candidates.size(); ++i)
    {
      if (elite_candidates[i]->f < (best + TargetTolFun)) {
        potential_candidates.push_back(elite_candidates[i]);
      }
    }

    // for each candidate, check if it is a novel local optimum
    for (size_t i = 0; i < potential_candidates.size(); ++i)
    {

      // check if the potential global optima is novel
      bool novel = true;

      for (size_t j = 0; j < elitist_archive.size(); ++j)
      {

        // a valid edge (check_edge returns true) suggest that the two optima are the same. 
        if (check_edge(*elitist_archive[j], *potential_candidates[i], add_elites_max_trials))
        {
          novel = false;
          number_of_global_opts_found++;

          // replace the elite with the candidate if it is better
          if (solution_t::better_solution_via_pointers(potential_candidates[i], elitist_archive[j])) {
            elitist_archive[j] = std::make_shared<solution_t>(*potential_candidates[i]);
            elitist_archive[j]->elite = true;
            elitist_archive[j]->time_obtained = ((double) (clock() - starting_time)) / CLOCKS_PER_SEC * 1000.0;
            elitist_archive[j]->feval_obtained = number_of_evaluations;
            elitist_archive[j]->generation_obtained = potential_candidates[i]->generation_obtained;
          }

          break;
        }
      }

      // it is novel, add it to the archive
      if (novel) {
        elitist_archive.push_back(std::make_shared<solution_t>(*potential_candidates[i]));
        elitist_archive.back()->elite = true;
        elitist_archive.back()->time_obtained = ((double)(clock() - starting_time)) / CLOCKS_PER_SEC * 1000.0;
        elitist_archive.back()->feval_obtained = number_of_evaluations;
        elitist_archive.back()->generation_obtained = potential_candidates[i]->generation_obtained;
        number_of_new_global_opts_found++;
      }

    }

  }

}



// returns true if it is a valid edge. (if the solutions belong to the same basin)
bool hillvallea::hillvallea_t::check_edge(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, int max_trials)
{
  std::vector<solution_pt> test_points;
  return check_edge(sol1, sol2, max_trials, test_points);
}

bool hillvallea::hillvallea_t::check_edge(const hillvallea::solution_t &sol1, const hillvallea::solution_t &sol2, int max_trials, std::vector<solution_pt> & test_points)
{

  // elites are already checked to be in different basins
  if (sol1.elite && sol2.elite) {
    return false;
  }

  // check max_trials with number_of_evaluations remaining
  if (maximum_number_of_evaluations > 0 && max_trials > maximum_number_of_evaluations - number_of_evaluations) {
    max_trials = maximum_number_of_evaluations - number_of_evaluations;
  }

  // find the worst solution of the two. 
  solution_t worst;
  if (solution_t::better_solution(sol1, sol2)) {
    worst = sol2;
  }
  else {
    worst = sol1;
  }

  for (size_t k = 0; k < max_trials; k++)
  {

    solution_pt x_test = std::make_shared<solution_t>();

    x_test->param = sol1.param + ((k + 1.0) / (max_trials + 1.0)) * (sol2.param - sol1.param);

    (*fitness_function)(*x_test);
    number_of_evaluations++;
    number_of_evaluations_clustering++;

    test_points.push_back(x_test);

    // if f[i] is better than f_test, we don't like the connection. So we stop.
    if (solution_t::better_solution(worst, *x_test)) {
      return false;
    }
  }

  return true;

}

void hillvallea::hillvallea_t::clustering(population_t & pop, std::vector<population_pt> & clusters, const std::vector<solution_pt> & elitist_archive, double & average_edge_length)
{

  // reset all clusters
  clusters.clear();

  // exclude the trivial case
  if (pop.size() == 0) {
    return;
  }

  // add elites to the archive, and mark them as elite
  if (elitist_archive.size() > 0)
  {
    for (size_t i = 0; i < elitist_archive.size(); ++i)
    {
      elitist_archive[i]->elite = true;
      pop.sols.push_back(std::make_shared<solution_t>(*elitist_archive[i]));
    }
    pop.sort_on_fitness();
  }


  // Initialize the first cluster as the cluster containing the best (=first) solution
  //-------------------------------------------------------------------------------------
  size_t number_of_clusters = 1;
  std::vector<size_t> cluster_index(pop.size(), -1);
  cluster_index[0] = 0;

  // remember how many solutions in the population created, so that we can allocate them later. 
  std::vector<solution_pt> test_points;
  std::vector<size_t> cluster_index_of_test_points;
  average_edge_length = pow(search_volume / pop.size(), 1.0 / number_of_parameters);
  double* dist = (double *)Malloc((long)pop.size() * sizeof(double));

  for (size_t i = 1; i < pop.size(); i++)
  {

    // compute the distance to all better solutions. 
    dist[i] = 0.0;
    size_t nearest_better_index = 0, worst_better_index = 0;
    for (size_t j = 0; j < i; j++) {
      dist[j] = pop.sols[i]->param_distance(*pop.sols[j]);

      if (dist[j] < dist[nearest_better_index]) {
        nearest_better_index = j;
      }

      if (dist[j] > dist[worst_better_index]) {
        worst_better_index = j;
      }
    }

    // Check neighbours
    bool edge_added = false;
    std::vector<size_t> does_not_belong_to(clustering_max_number_of_neighbours, -1);
    size_t old_nearest_better_index;
    std::vector<solution_pt> new_test_points_for_this_sol;

    for (size_t j = 0; j < std::min(i, clustering_max_number_of_neighbours); j++)
    {

      // find the next-to nearest index
      if (j > 0) 
      {
        old_nearest_better_index = nearest_better_index;
        nearest_better_index = worst_better_index;

        for (size_t k = 0; k < i; k++) {

          if (dist[k] > dist[old_nearest_better_index] && dist[k] < dist[nearest_better_index]) {
            nearest_better_index = k;
          }
        }
      }

      bool skip_neighbour = false;
      for (size_t k = 0; k < does_not_belong_to.size(); ++k)
      {
        if (does_not_belong_to[k] == cluster_index[nearest_better_index])
        {
          skip_neighbour = true;
          break;
        }
      }

      if (skip_neighbour) {
        continue;
      }

      int max_number_of_trial_solutions = 1 + ((int)(dist[nearest_better_index] / average_edge_length));
      std::vector<solution_pt> new_test_points;

      if (check_edge(*pop.sols[i], *pop.sols[nearest_better_index], max_number_of_trial_solutions, new_test_points))
      {
        cluster_index[i] = cluster_index[nearest_better_index];
        edge_added = true;

        // if the edge is accepted, add all test_poitns to their cluster
        for (size_t k = 0; k < new_test_points.size(); ++k) {
          test_points.push_back(new_test_points[k]);
          cluster_index_of_test_points.push_back(cluster_index[nearest_better_index]);
        }

        break;
      }
      else
      {
        does_not_belong_to[j] = cluster_index[nearest_better_index];

        // if the edge is not accepted, add all solutions to that cluster
        // all but the last because that one caused the rejection
        if (new_test_points.size() > 0) {
          for (size_t k = 0; k < new_test_points.size() - 1; ++k) {
            new_test_points_for_this_sol.push_back(new_test_points[k]);
          }
        }

      }
    }

    // its a new clusters, label it like that. 
    if (!edge_added)
    {
      cluster_index[i] = number_of_clusters;
      number_of_clusters++;

      // if its a new cluster, add all its testpoints as well. 
      for (size_t k = 0; k < new_test_points_for_this_sol.size(); ++k)
      {
        test_points.push_back(new_test_points_for_this_sol[k]);
        cluster_index_of_test_points.push_back(cluster_index[i]);
      }

    }

  }

  // create & fill the clusters
  //---------------------------------------------------------------------------
  std::vector<population_pt> candidate_clusters(number_of_clusters);
  std::vector<bool> cluster_active(number_of_clusters, true);

  for (size_t i = 0; i < number_of_clusters; ++i) {
    candidate_clusters[i] = std::make_shared<population_t>();
  }

  for (size_t i = 0; i < cluster_index.size(); ++i) {
    candidate_clusters[cluster_index[i]]->sols.push_back(pop.sols[i]);

    // only if the elite is the best of that population, we do not run it again.
    // don't get confused by this. We disable the cluster as soon as the first solution is an elite. 
    if (candidate_clusters[cluster_index[i]]->sols.size() == 1 && pop.sols[i]->elite) {
      cluster_active[cluster_index[i]] = false;
    }
  }

  for (size_t i = 0; i < test_points.size(); ++i) {
    candidate_clusters[cluster_index_of_test_points[i]]->sols.push_back(test_points[i]);
  }

  for (size_t i = 0; i < candidate_clusters.size(); ++i)
  {
    if (cluster_active[i]) {
      clusters.push_back(candidate_clusters[i]);
    }
  }

}



