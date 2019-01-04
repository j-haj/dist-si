#ifndef PARTICLES_FUNCTIONAL_H__
#define PARTICLES_FUNCTIONAL_H__

#include <cctype>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "particle_update_mode.h"
#include "util.h"

template <typename T>
class ParticlesFunctional {
 public:
  ParticlesFunctional(std::function<T(const std::vector<T>&)> f, ParticleUpdateMode mode,
		      std::size_t n_particles, std::size_t dim, T x_min, T x_max, T v_min,
		      T v_max, T omega, T c1, T c2)
      : mode_(mode),
        n_particles_(n_particles),
        dim_(dim),
        fitness_(f),
        x_min_(x_min),
        x_max_(x_max),
        v_min_(v_min),
        v_max_(v_max),
        omega_(omega),
        c1_(c1),
        c2_(c2),
        positions_(std::vector<std::vector<T>>(n_particles)),
        best_positions_(positions_),
        velocities_(positions_),
        gbest_(std::vector<T>(dim_)) {
    for (std::size_t i = 0; i < n_particles; ++i) {
      positions_[i] = util::uniform_rand_vec(dim, x_min, x_max);
    }
    best_positions_ = positions_;
    velocities_ = std::vector<std::vector<T>>(n_particles);
    for (auto& v : velocities_) {
      v = std::vector<T>(dim);
    }
    gbest_ = positions_[0];
  }

  ParticlesFunctional()
      : ParticlesFunctional([](const std::vector<T>& v) { return static_cast<T>(0); },
                  ParticleUpdateMode::Sequential, 1, 1, 0, 0, 0, 0, 0, 0, 0) {}

  ParticlesFunctional<T>& operator=(Particles<T>& p) {
    mode_ = p.mode_;
    n_particles_ = p.n_particles_;
    dim_ = p.dim_;
    fitness_ = p.fitness_;
    x_min_ = p.x_min_;
    x_max_ = p.x_max_;
    v_min_ = p.v_min_;
    v_max_ = p.v_max_;
    omega_ = p.omega_;
    c1_ = p.c1_;
    c2_ = p.c2_;
    positions_ = p.positions_;
    best_positions_ = p.best_positions_;
    velocities_ = p.velocities_;
    gbest_ = p.gbest_;
    return *this;
  }

  void UpdatePositions() noexcept {
    switch (mode_) {
      case ParticleUpdateMode::Sequential:
        UpdatePositionsSequential();
        break;
      case ParticleUpdateMode::Parallel:
        UpdatePositionsParallel();
        break;

    }
  }

  void UpdateVelocities() noexcept {
    switch (mode_) {
      case ParticleUpdateMode::Sequential:
        UpdateVelocitiesSequential();
        break;
      case ParticleUpdateMode::Parallel:
        UpdateVelocitiesParallel();
        break;
    }
  }

  const std::vector<T>& gbest() const noexcept { return gbest_; }
  std::size_t size() const noexcept { return positions_.size(); }

 private:
  void UpdateVelocitiesSequential() noexcept;
  void UpdateVelocitiesParallel() noexcept;
  void UpdatePositionsSequential() noexcept;
  void UpdatePositionsParallel() noexcept;
  void UpdatePosition(std::vector<T>& pos, std::vector<T>& v) noexcept;
  void UpdateVelocity(std::vector<T>& velocity,
		      const std::vector<T>& best_pos,
		      const std::vector<T>& current_pos) noexcept;

  ParticleUpdateMode mode_;
  std::size_t n_particles_;
  std::size_t dim_;
  std::function<T(const std::vector<T>&)> fitness_;
  T x_min_;
  T x_max_;
  T v_min_;
  T v_max_;
  T omega_;
  T c1_;
  T c2_;

  std::vector<std::vector<T>> positions_;
  std::vector<std::vector<T>> best_positions_;
  std::vector<std::vector<T>> velocities_;
  std::vector<T> gbest_;

};  // class Particles

template <typename T>
void ParticlesFunctional<T>::UpdateVelocity(std::vector<T>& velocity,
					    const std::vector<T>& best_pos,
					    const std::vector<T>& current_pos) noexcept {
  const auto r1 = util::uniform_unit<T>();
  const auto r2 = util::uniform_unit<T>();
  for (std::size_t i = 0; i < dim_; ++i) {
    velocity[i] = velocity[i] * omega_ +
      c1_ * r1 * (best_pos[i] - current_pos[i]) +
      c2_ * r2 * (gbest_[i] - current_pos[i]);
    // Clamp velocity
    velocity[i] = std::min(std::max(v_min_, velocity[i]), v_max_);
  }
}

template <typename T>
void ParticlesFunctional<T>::UpdateVelocitiesSequential() noexcept {
  for (std::size_t i = 0; i < velocities_.size(); ++i) {
    UpdateVelocity(velocities_[i], best_positions_[i], positions_[i]);
  }
}

template <typename T>
void ParticlesFunctional<T>::UpdateVelocitiesParallel() noexcept {
#pragma omp parallel for schedule(runtime)
  for (std::size_t i = 0; i < velocities_.size(); ++i) {
    UpdateVelocity(velocities_[i], best_positions_[i], positions_[i]);
  }
}

template <typename T>
void ParticlesFunctional<T>::UpdatePosition(std::vector<T>& pos,
					    std::vector<T>& v) noexcept {
  for (std::size_t i = 0; i < dim_; ++i) {
    pos[i] += v[i];
    if (pos[i] > x_max_) {
      pos[i] -= (pos[i] - x_max_);
      v[i] *= -1.0;
    }
  }
}

template <typename T>
void ParticlesFunctional<T>::UpdatePositionsSequential() noexcept {
  T best_fitness = fitness_(gbest_);
  for (std::size_t i = 0; i < positions_.size(); ++i) {
    UpdatePosition(positions_[i], velocities_[i]);

    // Maintain local best position and gbest
    const T fitness = fitness_(positions_[i]);
    if (fitness < fitness_(best_positions_[i])) {
      best_positions_[i] = positions_[i];
      if (fitness < best_fitness) {
	gbest_ = positions_[i];
	best_fitness = fitness;
      }
    }
  }
}

template <typename T>
void ParticlesFunctional<T>::UpdatePositionsParallel() noexcept {
  T best_fitness = fitness_(gbest_);
#pragma omp parallel for schedule(runtime)
  for (std::size_t i = 0; i < positions_.size(); ++i) {
    UpdatePosition(positions_[i], velocities_[i]);
  
    // Maintain local best position and gbest
    const T fitness = fitness_(positions_[i]);
    if (fitness < fitness_(best_positions_[i])) {
      best_positions_[i] = positions_[i];
    }
    if (fitness < best_fitness) {
#pragma omp critical
      if (fitness < best_fitness) {
	gbest_ = positions_[i];
	best_fitness = fitness;
      }
    }
  }
}


#endif  // PARTICLES_FUNCTIONAL_H__
