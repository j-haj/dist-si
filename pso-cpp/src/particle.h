#ifndef PARTICLE_H__
#define PARTiCLE_H__

#include <vector>

template <typename T>
class Particle {
public:
  // TODO: properly initialize vectors with random values for position
  // and zero values for velocity.
  Particle(std::size_t n, T omega, T c1, T c2)
    : current_position_(std::vector<T>(n)),
    local_best_pos_(std::vector<T>(n)),
    velocity_(std::vector<T>(n)),
    omega_(omega),
    c1_(c1),
    c2_(c2) {}

  void update_velocity(const std::vector<T>& gbest) noexcept;
  void update_position() noexcept;
  
private:
  std::vector<T> current_position_;
  std::vector<T> local_best_position_;
  std::vector<T> velocity_;
  T omega_;
  T c1_;
  T c2_;

}; // class Particle

template <typename T>
void Particle::pdateVelocity(const std::vector<T>& gbest) noexcept {
  // TODO
}

template <typename T>
void Particle::UpdatePosition() noexcept {
  for (std::size_t i = 0; i < velocity_.size(); ++i) {
    current_position_[i] += velocity_[i];
}
#endif
