#include <iostream>

#include <omp.h>

#include "experiment.h"
#include "particles.h"

int main() {
  using type_t = double;
  type_t omega = 1.8;
  type_t c1 = 1.0;
  type_t c2 = 1.0;
  type_t x_min = -10.0;
  type_t x_max = 10.0;
  type_t v_min = -10.0;
  type_t v_max = 10.0;
  std::size_t n_particles = 1000;
  std::size_t dim = 100;

  omp_set_num_threads(12);
  
  auto fitness = [](const std::vector<double>& p)->type_t {
		   double v = 0.0;
		   for (const auto& x : p) {
		     v += x*x;
		   }
		   return v;
		};
  auto max_steps = 10000;

  std::cout << "Constructing a Particles object...\n";
  auto p = Particles<type_t>(fitness, ParticleUpdateMode::Sequential, n_particles,
			dim, x_min, x_max, v_min, v_max, omega, c1, c2);
  std::cout << "Success!\n";
  SoAExperiment<double> soa_experiment(fitness,
				       ParticleUpdateMode::Sequential,
				       -10.0,
				       10.0,
				       -10.0,
				       10.0,
				       1000,
				       100,
				       omega,
				       c1,
				       c2);
  soa_experiment.Run(max_steps);  
  AoSExperiment<double> experiment(fitness,
				   RunMode::Sequential,
				   -10.0,
				   10.0,
				   -10.0,
				   10.0,
				   1000,
				   100,
				   omega,
				   c1,
				   c2);
  experiment.Run(max_steps);

  return 0;
}
