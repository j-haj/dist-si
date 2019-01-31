// Inside a for-loop with index variable i
const double gbest_fitness = fitness_(gbest_);
const double fitness = fitness_(positions_[i]);
if (fitness > fitness_(best_positions_[i])) {
  best_positions_[i] = positions_[i];
  if (fitness > gbest_fitness) {
#pragma omp critical
    if (fitness > gbest_fitness) {
      gbest_ = positions_[i];
    }
  }
 }
