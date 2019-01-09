#include <cmath>
#include <iostream>

#include <omp.h>

#include "experiment.h"
#include "particles.h"
#include "particles_functional.h"
#include "particle_update_mode.h"

enum class RunType {
  Sequential,
  Parallel,
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

int main() {
  using type_t = double;
  type_t omega = 1.8;
  type_t c1 = 1.0;
  type_t c2 = 1.2;
  type_t x_min = -50.0;
  type_t x_max = 50.0;
  type_t v_min = -50.0;
  type_t v_max = 50.0;

  // Params
  std::size_t n_particles = 10000;
  std::size_t dim = 500;
  RunType run_type = RunType::Sequential;
  bool verbose = false;

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
      sum1 += 0.5 * x * x;
      sum2 += 0.5 * std::cos(pi * 2 * x);
    }
    return -20 * std::exp(-0.2 * std::sqrt(sum1)) - std::exp(sum2) + std::exp(1) + 20;
  };

  auto fitness = rastrigin;
  auto max_steps = 10000;
  type_t epsilon = 1e-8;
  int n_experiments = 20;
  std::vector<type_t> solution(dim, 0.0);
  std::vector<double> soa_times(n_experiments);
  std::vector<double> soa_functional_times(n_experiments);
  std::vector<double> aos_times(n_experiments);
  int completed = 0;
  std::vector<double> round_times;
  for (int i = 0; i < n_experiments; ++i) {
    SoAExperiment<Particles, double>
      soa_experiment(fitness, particle_update_mode,
		     x_min, x_max, v_min, v_max, n_particles,
		     dim, omega, c1, c2);
    soa_times[i] = soa_experiment.Run(max_steps, solution, epsilon);
    SoAExperiment<ParticlesFunctional, double>
      soa_functional_experiment(fitness, particle_update_mode,
				x_min, x_max, v_min, v_max, n_particles,
				dim, omega, c1, c2);
    soa_functional_times[i] = soa_functional_experiment.Run(max_steps, solution,
							    epsilon);
    AoSExperiment<double> aos_experiment(fitness, run_mode, x_min, x_max,
					 v_min, v_max, n_particles, dim, omega, c1,
					 c2);
    aos_times[i] = aos_experiment.Run(max_steps, solution, epsilon);
    if (verbose) {
      std::cout << "SoA: " << soa_times[i] << '\n'
		<< "SoAf: " << soa_functional_times[i] << '\n'
		<< "AoS: " << aos_times[i] << '\n';
    }
    ++completed;
    round_times.push_back(soa_times[i] + soa_functional_times[i] + aos_times[i]);
    auto avg = average(round_times);
    std::cout << "Completed " << completed << " of " << n_experiments << " rounds.\n"
	      << "Average round time: " << avg << "s\n"
	      << "Approximate time remaining: " << avg * (n_experiments - completed)
	      << "s\n\n";
  }

  std::cout << "Done.\nStruct of Arrays\n\tavg (s): " << average(soa_times)
	    << "\n\tstd dev (s): " << std_dev(soa_times)
	    << "\nStruct of Arrays - Functional\n\tavg (s): " << average(soa_functional_times)
	    << "\n\tstd dev (s): " << std_dev(soa_functional_times)
	    << "\nArray of Structs\n\tavg (s): " << average(aos_times)
	    << "\n\tstd dev (s): " << std_dev(aos_times) << '\n';
  return 0;
}
