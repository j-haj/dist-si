#pragma omp parallel for collapse(2)
for (std::size_t i  = 0; i < n_particles; ++i) {
  for (std::size_t j = 0; j < dimension; ++j) {
    const auto r1 = random_uniform(0, 1);
    const auto r2 = random_uniform(0, 1);
    velocities[i][j] = omega * velocities[i][j] +
      c1 * r1 * (positions[i][j] - best_positions[i][j]) +
      c2 * r2 * (positions[i][j] - gbest[j]);
  }
}
