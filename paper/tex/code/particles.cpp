template <typename T>
struct  Particles {
  using vecT = std::vector<T>;
  std::vector<vecT> positions;
  std::vector<vecT> best_positions;
  std::vector<vecT> velocities;
  vecT global_best;
  T c1;
  T c2;
  T omega;
};
