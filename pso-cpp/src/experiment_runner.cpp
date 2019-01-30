#include <cmath>
#include <iostream>
#include <sstream>

#include <omp.h>

#include "experiment.h"
#include "particles.h"
#include "particles_functional.h"
#include "particle_update_mode.h"

#define N_PARTICLES 500
#define DIMENSION 10
#define MAX_STEPS 1000000
#define EPSILON 5

enum class RunType {
  Sequential,
  Parallel,
};

enum class TestFunction {
  Rastrigin,
  Ackley,
  Quadratic,
};

double average(const std::vector<double>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double std_dev(const std::vector<double>& v) {
  auto sum = 0.0;
  const auto avg = average(v);
  for (const auto x : v) {
    const auto d = avg - x;
    sum += d * d;
  }
  return std::sqrt(sum/(v.size() - 1));
}

struct SummaryStatistics {
  SummaryStatistics(const std::vector<Result>& results)
    : avg_time(0.0), dev_time(0.0), avg_fitness(0.0),
      dev_fitness(0.0), avg_steps(0.0), dev_steps(0.0)
  {
    std::vector<double> times(results.size());
    std::vector<double> fitnesses(results.size());
    std::vector<double> steps(results.size());
    for (std::size_t i = 0; i < results.size(); ++i) {
      times[i] = results[i].time;
      fitnesses[i] = results[i].fitness;
      steps[i] = results[i].n_steps;
    }
    avg_time = average(times);
    dev_time = std_dev(times);
    avg_fitness = average(fitnesses);
    dev_fitness = std_dev(fitnesses);
    avg_steps = average(steps);
    dev_steps = std_dev(steps);
  }
	       
  double avg_time = 0.0;
  double dev_time = 0.0;
  double avg_fitness = 0.0;
  double dev_fitness = 0.0;
  double avg_steps = 0.0;
  double dev_steps = 0.0;
};

void print_summary(const std::string& title, const SummaryStatistics& summary) {
  const std::string pm = " +/- ";
  std::cout << title;
  std::cout << "\n\ttime (s): " << summary.avg_time << pm << summary.dev_time
	    << "\n\tfitness: " << summary.avg_fitness << pm << summary.dev_fitness
	    << "\n\tsteps: " << summary.avg_steps << pm << summary.dev_steps
	    << std::endl;
}

int main() {
  using type_t = double;
  const type_t chi = 2.0 / std::abs(2 - 4.1 - std::sqrt(4.1 * 4.1 - 4 * 4.1));
  const type_t omega = chi;
  const type_t c1 = 2.05;
  const type_t c2 = 2.05;
  const type_t x_min = -5.12;
  const type_t x_max = 5.12;
  const type_t v_min = 2 * x_min;
  const type_t v_max = 2 * x_max;


  // Params
  std::size_t n_particles = N_PARTICLES;
  std::size_t dim = DIMENSION;
  RunType run_type = RunType::Parallel;
  bool verbose = true;
  const auto max_steps = MAX_STEPS;
  auto test_function = TestFunction::Rastrigin;
  
  ParticleUpdateMode particle_update_mode;
  RunMode run_mode;
  switch (run_type) {
  case RunType::Sequential:
    std::cout << "Running in Sequential mode\n";
    particle_update_mode = ParticleUpdateMode::Sequential;
    run_mode = RunMode::Sequential;
    break;
  case RunType::Parallel:
    std::cout << "Running in Parallel mode\n";
    particle_update_mode = ParticleUpdateMode::Parallel;
    run_mode = RunMode::Parallel;
    break;
  }

  omp_set_num_threads(6);

  auto quadratic = [](const std::vector<type_t>& p) -> type_t {
    double v = 0.0;
    for (const auto& x : p) {
      v += x * x;
    }
    return v;
  };

  auto rastrigin = [](const std::vector<type_t>& p) -> type_t {
    const auto A = 10.0;
    static const auto pi = std::acos(-1);
    auto v = A * p.size();
    for (const auto& x : p) {
      v += x * x - A * std::cos(2 * pi * x);
    }
    return v;
  };

  auto ackley = [](const std::vector<type_t>& p) -> type_t {
    static const auto pi = std::acos(-1);
    auto sum1 = 0.0;
    auto sum2 = 0.0;
    for (const auto& x : p) {
      sum1 += x * x;
      sum2 += std::cos(pi * 2 * x);
    }
    return -20 * std::exp(-0.2 * std::sqrt(sum1/p.size())) - std::exp(sum2/p.size()) + std::exp(1) + 20;
  };
  // Silence unused warning...
  (void)quadratic;
  (void)rastrigin;
  (void)ackley;

  std::function<type_t(const std::vector<type_t>&)> fitness;
  std::string func_name;
  if (test_function == TestFunction::Quadratic) {
    func_name == "Quadratic";
    fitness = quadratic;
  } else if (test_function == TestFunction::Rastrigin) {
    func_name == "Rastrigin";
    fitness = rastrigin;
  } else if (test_function == TestFunction::Ackley) {
    func_name == "Ackley";
    fitness = ackley;
  }

  type_t epsilon = EPSILON;
  int n_experiments = 20;

  std::vector<Result> soa_results(n_experiments);
  std::vector<Result> soaf_results(n_experiments);
  std::vector<Result> aos_results(n_experiments);
  int completed = 0;
  std::vector<double> round_times;
  for (int i = 0; i < n_experiments; ++i) {
    SoAExperiment<Particles, double>
      soa_experiment(fitness, particle_update_mode,
    		     x_min, x_max, v_min, v_max, n_particles,
    		     dim, omega, c1, c2);
    soa_results[i] = soa_experiment.Run(max_steps, epsilon);
    SoAExperiment<ParticlesFunctional, double>
      soa_functional_experiment(fitness, particle_update_mode,
				x_min, x_max, v_min, v_max, n_particles,
				dim, omega, c1, c2);
    soaf_results[i] = soa_functional_experiment.Run(max_steps, epsilon);
    AoSExperiment<double> aos_experiment(fitness, run_mode, x_min, x_max,
					 v_min, v_max, n_particles, dim, omega, c1,
					 c2);
    aos_results[i] = aos_experiment.Run(max_steps, epsilon);
    if (verbose) {
      std::cout << "SoA: " << soa_results[i] << '\n'
		<< "SoAf: " << soaf_results[i] << '\n'
		<< "AoS: " << aos_results[i] << '\n';
    }
    ++completed;
    round_times.push_back(soa_results[i].time + soaf_results[i].time + aos_results[i].time);
    auto avg = average(round_times);
    std::cout << "Completed " << completed << " of " << n_experiments << " rounds.\n"
	      << "Average round time: " << avg << "s\n"
	      << "Approximate time remaining: " << avg * (n_experiments - completed)
	      << "s\n\n";
  }

  SummaryStatistics soa_summary(soa_results);
  SummaryStatistics soaf_summary(soaf_results);
  SummaryStatistics aos_summary(aos_results);
  const std::string pm = " +/- ";
  std::cout << "Done.\n\n";
  std::cout << "Number of particles: " << n_particles << '\n';
  std::cout << "Dimension: " << dim << '\n';
  std::cout << "Function: " << func_name << '\n';
  print_summary("AoS", aos_summary);
  print_summary("SoA", soa_summary);
  print_summary("SoAf", soaf_summary);
  return 0;
}
