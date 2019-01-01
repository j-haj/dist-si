#include <iostream>

#include "experiment.h"

int main() {
  double omega = 1.8;
  double c1 = 1.0;
  double c2 = 1.0;

  auto fitness = [](const Particle<double>& p) {
		   double v = 0.0;
		   for (const auto& x : p.Coordinates()) {
		     v += x*x;
		   }
		   return v;
		};
  AoSExperiment<double> experiment(fitness,
				   RunMode::Sequential,
				   -10.0,
				   10.0,
				   -1.0,
				   1.0,
				   10,
				   5,
				   omega,
				   c1,
				   c2);
  auto max_steps = 10000;
  experiment.Run(max_steps);
  return 0;
}
