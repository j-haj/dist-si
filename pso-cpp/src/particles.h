#ifndef PARTICLES_H__
#define PARTICLES_H__

#include <cctype>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "util.h"

enum class ParticleUpdateMode {
  Sequential,
  Parallel,
  AVX,
};  // enum class ParticleUpdateMode

template <typename T>
class Particles {
 public:
  Particles(std::function<T(const std::vector<T>&)> f, ParticleUpdateMode mode,
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
  }

  Particles()
      : Particles([](const std::vector<T>& v) { return static_cast<T>(0); },
                  ParticleUpdateMode::Sequential, 1, 1, 0, 0, 0, 0, 0, 0, 0) {}

  Particles<T>& operator=(Particles<T>& p) {
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
      case ParticleUpdateMode::AVX:
        UpdatePositionsAVX();
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
      case ParticleUpdateMode::AVX:
        UpdateVelocitiesAVX();
        break;
    }
  }

  const std::vector<T>& gbest() const noexcept { return gbest_; }
  std::size_t size() const noexcept { return positions_.size(); }

 private:
  void UpdateVelocitiesSequential() noexcept;
  void UpdateVelocitiesParallel() noexcept;
  void UpdateVelocitiesAVX() noexcept;
  void UpdatePositionsSequential() noexcept;
  void UpdatePositionsParallel() noexcept;
  void UpdatePositionsAVX() noexcept;

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
void Particles<T>::UpdateVelocitiesSequential() noexcept {
  for (std::size_t i = 0; i < velocities_.size(); ++i) {
    const auto r1 = util::uniform_unit<T>();
    const auto r2 = util::uniform_unit<T>();
    for (std::size_t j = 0; j < dim_; ++j) {
      velocities_[i][j] =
          velocities_[i][j] * omega_ +
          c1_ * r1 * (best_positions_[i][j] - positions_[i][j]) +
          c2_ * r2 * (gbest_[j] - positions_[i][j]);
      // Clamp velocities
      velocities_[i][j] = std::min(std::max(v_min_, velocities_[i][j]), v_max_);
    }
  }
}

template <typename T>
void Particles<T>::UpdateVelocitiesParallel() noexcept {
#pragma omp parallel for
  for (std::size_t i = 0; i < velocities_.size(); ++i) {
    const auto r1 = util::uniform_unit<T>();
    const auto r2 = util::uniform_unit<T>();
    for (std::size_t j = 0; j < dim_; ++j) {
      velocities_[i][j] =
          velocities_[i][j] * omega_ +
          c1_ * r1 * (best_positions_[i][j] - positions_[i][j]) +
          c2_ * r2 * (gbest_[j] - positions_[i][j]);
      // Clamp velocities
      velocities_[i][j] = std::min(std::max(v_min_, velocities_[i][j]), v_max_);
    }
  }
}

template <typename T>
void Particles<T>::UpdateVelocitiesAVX() noexcept {}

template <typename T>
void Particles<T>::UpdatePositionsSequential() noexcept {
  T best_fitness = fitness_(gbest_);
  for (std::size_t i = 0; i < positions_.size(); ++i) {
    for (std::size_t j = 0; j < dim_; ++j) {
      positions_[i][j] += velocities_[i][j];
      if (positions_[i][j] > x_max_) {
        positions_[i][j] -= positions_[i][j] - x_max_;
        velocities_[i][j] *= -1.0;
      }
    }

    // Maintain local best position and gbest
    const T fitness = fitness_(positions_[i]);
    if (fitness < fitness_(best_positions_[i])) {
      best_positions_[i] = positions_[i];
    }
    if (fitness < best_fitness) {
      gbest_ = positions_[i];
      best_fitness = fitness;
    }
  }
}

template <typename T>
void Particles<T>::UpdatePositionsParallel() noexcept {
  T best_fitness = fitness_(gbest_);
#pragma omp parallel for collapse(2)
  for (std::size_t i = 0; i < positions_.size(); ++i) {
    for (std::size_t j = 0; j < dim_; ++j) {
      positions_[i][j] += velocities_[i][j];
      if (positions_[i][j] > x_max_) {
        positions_[i][j] -= positions_[i][j] - x_max_;
        velocities_[i][j] *= -1.0;
      }
    }
  }
#pragma omp parallel for
  for (std::size_t i = 0; i < positions_.size(); ++i) {
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

template <typename T>
void Particles<T>::UpdatePositionsAVX() noexcept {}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  os << "vector { ";
  for (const auto& x : v) {
    os << x << ' ';
  }
  os << "}";
  return os;
}

#endif  // PARTICLES_H__
