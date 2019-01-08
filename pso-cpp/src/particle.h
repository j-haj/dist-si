#ifndef PARTICLE_H__
#define PARTICLE_H__

#include <cctype>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "util.h"

template <typename T>
class Particle {
 public:
  Particle(T x_min, T x_max, T v_min, T v_max, std::size_t n, T omega, T c1,
           T c2, std::function<T(const std::vector<T>&)> f)
      : current_pos_(util::uniform_rand_vec<T>(n, x_min, x_max)),
        local_best_pos_(current_pos_),
        velocity_(std::vector<T>(n, 0)),
        fitness_(f),
        x_min_(x_min),
        x_max_(x_max),
        v_min_(v_min),
        v_max_(v_max),
        omega_(omega),
        c1_(c1),
        c2_(c2) {}

  Particle()
      : current_pos_(util::uniform_rand_vec<T>(1, 0, 1)),
        local_best_pos_(current_pos_),
        velocity_(std::vector<T>()),
	fitness_([](const std::vector<T>& p) { return static_cast<T>(p.size()); }),
        x_min_(0),
        x_max_(0),
        v_min_(0),
        v_max_(0),
        omega_(0),
        c1_(0),
        c2_(0) {}

  Particle<T>& operator=(const Particle<T>& p) {
    current_pos_ = p.current_pos_;
    local_best_pos_ = p.local_best_pos_;
    velocity_ = p.velocity_;
    fitness_ = p.fitness_;
    x_min_ = p.x_min_;
    x_max_ = p.x_max_;
    v_min_ = p.v_min_;
    v_max_ = p.v_max_;
    omega_ = p.omega_;
    c1_ = p.c1_;
    c2_ = p.c2_;
    return *this;
  }

  void UpdateVelocity(const Particle<T>& gbest) noexcept;
  void UpdatePosition() noexcept;

  const std::vector<T>& Coordinates() const noexcept { return current_pos_; }

  const std::vector<T>& Velocity() const noexcept { return velocity_; }

  T Fitness() const noexcept { return fitness_(current_pos_); }

 private:
  std::vector<T> current_pos_;
  std::vector<T> local_best_pos_;
  std::vector<T> velocity_;
  std::function<T(const std::vector<T>&)> fitness_;
  T x_min_;
  T x_max_;
  T v_min_;
  T v_max_;
  T omega_;
  T c1_;
  T c2_;

};  // class Particle

template <typename T>
void Particle<T>::UpdateVelocity(const Particle<T>& gbest) noexcept {
  const auto r1 = util::uniform_unit<T>();
  const auto r2 = util::uniform_unit<T>();
  for (std::size_t i = 0; i < velocity_.size(); ++i) {
    velocity_[i] = velocity_[i] * omega_ +
                   c1_ * r1 * (local_best_pos_[i] - current_pos_[i]) +
                   c2_ * r2 * (gbest.Coordinates()[i] - current_pos_[i]);
    // Clamp velocity at min/max
    velocity_[i] = std::min(std::max(v_min_, velocity_[i]), v_max_);
  }
}

template <typename T>
void Particle<T>::UpdatePosition() noexcept {
  for (std::size_t i = 0; i < velocity_.size(); ++i) {
    current_pos_[i] += velocity_[i];
  }
  if (fitness_(current_pos_) < fitness_(local_best_pos_)) {
    local_best_pos_ = current_pos_;
  }
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Particle<T>& p) {
  os << "Particle{ ";
  for (const auto x : p.Coordinates()) {
    os << x << " ";
  }
  os << "}";
  return os;
}

#endif  // PARTICLE_H__
