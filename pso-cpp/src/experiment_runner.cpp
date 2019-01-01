#include <iostream>

#include <omp.h>

#include "experiment.h"

int main() {
  double omega = 1.8;
  double c1 = 1.0;
  double c2 = 1.0;

  omp_set_num_threads(12);
  
  auto fitness = [](const std::vector<double>& p) {
		   double v = 0.0;
		   for (const auto& x : p) {
		     v += x*x;
		   }
		   return v;
		};
  AoSExperiment<double> experiment(fitness,
				   RunMode::Parallel,
				   -10.0,
				   10.0,
				   -10.0,
				   10.0,
				   1000,
				   100,
				   omega,
				   c1,
				   c2);
  auto max_steps = 10000;
  experiment.Run(max_steps);
  return 0;
}
