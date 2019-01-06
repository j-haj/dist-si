template <typename T>
struct Particle {
  std::vector<T> position;
  std::vector<T> velocity;
  std::vector<T> local_best;
  T c1;
  T c2;
  T omega;
  std::function<T(const std::vector<T>&)> fitness;
  void UpdateVelocity(const std::vector<T>& gbest);
  void UpdatePosition();
};
