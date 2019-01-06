template <typename T>
void UpdateVelocity(std::vector<T>& v,
		    const std::vector<T>& pos,
		    const std::vector<T>& bpos,
		    const std::vector<T>& gbest)
{
  const T r1 = random_uniform(0, 1);
  const T r2 = random_uniform(0, 1);
  for (std::size_t i = 0; i < v.size(); ++i) {
    v[i] = omega * v[i] + c1 * r1 * (pos[i] - bpos[i])
      + c2 * r2 * (pos[i] - gbest[i]);
  }
}

template <typename T>
void UpdateVelocities() {
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < velocities.size(); ++i) {
    UpdateVelocity(v[i], pos[i], bpos[i], gbest);
  }
}