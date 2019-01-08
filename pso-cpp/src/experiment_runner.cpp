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

int main() {
  using type_t = double;
  type_t omega = 1.8;
  type_t c1 = 1.0;
  type_t c2 = 1.2;
  type_t x_min = -50.0;
  type_t x_max = 50.0;
  type_t v_min = -50.0;
  type_t v_max = 50.0;
  std::size_t n_particles = 1000;
  std::size_t dim = 250;
  RunType run_type = RunType::Sequential;
  
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
  for (int i = 0; i < n_experiments; ++i) {
    std::cout << "SoA... ";
    SoAExperiment<Particles, double>
      soa_experiment(fitness, particle_update_mode,
		     x_min, x_max, v_min, v_max, n_particles,
		     dim, omega, c1, c2);
    soa_times[i] = soa_experiment.Run(max_steps, solution, epsilon);
    std::cout << soa_times[i] << "s\n";
    std::cout << "SoA functional...";
    SoAExperiment<ParticlesFunctional, double>
      soa_functional_experiment(fitness, particle_update_mode,
				x_min, x_max, v_min, v_max, n_particles,
				dim, omega, c1, c2);
    soa_functional_times[i] = soa_functional_experiment.Run(max_steps, solution,
							    epsilon);
    std::cout << soa_functional_times[i] << "s\n";
    std::cout << "AoS...";
    AoSExperiment<double> aos_experiment(fitness, run_mode, x_min, x_max,
					 v_min, v_max, n_particles, dim, omega, c1,
					 c2);
    aos_times[i] = aos_experiment.Run(max_steps, solution, epsilon);
    std::cout << aos_times[i] << "s\n";
  }

  const auto soa_avg = std::accumulate(soa_times.begin(),
				       soa_times.end(), 0.0) / n_experiments;
  const auto soa_func_avg = std::accumulate(soa_functional_times.begin(),
					    soa_functional_times.end(), 0.0) / n_experiments;
  const auto aos_avg = std::accumulate(aos_times.begin(),
				       aos_times.end(), 0.0) / n_experiments;
  std::cout << "Done.\nStruct of Arrays\n\ttime (s): " << soa_avg
	    << "\nStruct of Arrays - Functional\n\ttime (s): " << soa_func_avg
	    << "\nArray of Structs\n\ttime (s): " << aos_avg
	    << "\n\nAoS / SoA = " << aos_avg / soa_avg << '\n';
  return 0;
}
