template <typename T>
void UpdateVelocity(std::vector<T>& v,
		    const std::vector<T>& pos,
		    const std::vector<T>& bpos,
		    const std::vector<T>& gbest,
		    T r1, T r2)
{
  for (std::size_t i = 0; i < v.size(); ++i) {
    v[i] = omega * v[i] + c1 * r1 * (pos[i] - bpos[i])
      + c2 * r2 * (pos[i] - gbest[i]);
  }
}

template <typename T>
void UpdateVelocities() {
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < velocities.size(); ++i) {
    UpdateVelocity(v[i], pos[i], bpos[i], gbest, r1s[i], r2s[i]);
  }
}
