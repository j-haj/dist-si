#ifndef PARTICLE_H__
#define PARTiCLE_H__

#include <random>
#include <vector>

#include "util.h"

template <typename T>
class Particle {
public:
  Particle(T x_min, T x_max, T v_min, T v_max,
	   std::size_t n, T omega, T c1, T c2)
    : current_pos_(util::uniform_rand_vec(n, x_min, x_max)),
      local_best_pos_(current_pos_),
      velocity_(std::vector<T>(n, 0)),
      x_min_(x_min),
      x_max_(x_max),
      v_min_(v_min),
      v_max_(v_max),
      omega_(omega),
      c1_(c1),
      c2_(c2) {
    rd_ = std::random_device();
    gen_ = std::mt19937_64(rd_);
    dist_ = std::uniform_real_distribution<T>(0, 1);
  }

  void UpdateVelocity(const Particle<T>& gbest) noexcept;
  void UpdatePosition() noexcept;
  
  T Coordinate(std::size_t i) const noexcept {
    return current_pos_[i];
  }
  
  T Velocity(std::size_t i) const noexcept {
    return velocity_[i];
  }
  
private:
  std::vector<T> current_pos_;
  std::vector<T> local_best_pos_;
  std::vector<T> velocity_;
  T x_min_;
  T x_max_;
  T v_min_;
  T v_max_;
  T omega_;
  T c1_;
  T c2_;
  std::random_device rd_;
  std::mt19937_64 gen_;
  std::uniform_real_distribution<T> dist_;
  

}; // class Particle

template <typename T>
void Particle<T>::UpdateVelocity(const Particle<T>& gbest) noexcept {
  const auto r1 = dist_(gen_);
  const auto r2 = dist_(gen_);
  for (std::size_t i = 0; i < velocity_.size(); ++i) {
    velocity_[i] = velocity_[i] * omega_
      + c1_ * r1 (local_best_pos_[i] - current_pos_[i])
      + c2_ * r2 (gbest[i] - current_pos_[i]);
    // Clamp velocity at min/max
    velocity_[i] = min(max(v_min_, velocity_[i]), v_max_);
  }
}

template <typename T>
void Particle<T>::UpdatePosition() noexcept {
  for (std::size_t i = 0; i < velocity_.size(); ++i) {
    current_pos_[i] += velocity_[i];
    // Reflect positions beyond boundary
    if (current_pos_[i] > x_max_) {
      current_pos_[i] -= current_pos_[i] - x_max_;
      velocity_[i] *= -1.0;
    }
  }
}

#endif
