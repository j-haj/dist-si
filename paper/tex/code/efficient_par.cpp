template <typename T>
void UpdateVelocities() {
  const auto r1s = util::rand_unit_vec(n_particles);
  const auto r2s = util::rand_unit_vec(n_particles);
#pragma omp parallel for collapse(2)
  for (std::size_t i = 0; i < n_particles; ++i) {
    for (std::size_t j = 0; j < dimension; ++j) {
      velocities[i][j] = omega * velocities[i][j] +
	c1 * r1s[i] * (positions[i][j] - best_pos[i][j]) +
	c2 * r2s[i] * (positions[i][j] - gbest[j]);
    }
  }
}
