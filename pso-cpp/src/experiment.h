#ifndef EXPERIMENT_H__
#define EXPERIMENT_H__

#include <cctype>
#include <chrono>
#include <functional>
#include <vector>

#include "particle.h"

enum class RunMode {
  Sequential,
  Parallel
};


template <typename T>
class AoSExperiment {

public:
  AoSExperiment(std::function<T(const Particle<T>&)> f,
		RunMode run_mode,
		T x_min, T x_max, T v_min, T v_max,
		std::size_t n_particles, std::size_t dim,
		T omega, T c1, T c2)
    : fitness_(f),
      particles_(std::vector<Particle<T>>(n_particles)),
      run_mode_(run_mode),
      n_particles_(n_particles),
      dim_(dim) {

    // Generate particles for simulation
    for (std::size_t i = 0; i < n_particles; ++i) {
      particles_[i] = Particle<T>(x_min, x_max, v_min, v_max,
				  dim, omega, c1, c2);
    }

    gbest_ = particles_[0];
    FindMinParticle();
  }

  void Run(std::size_t max_steps, T epsilon) {
    
    
  }

  
  
private:

  void FindMinParticle() noexcept {
    switch (run_mode_) {
    case RunMode::Sequential:
      FindMinParticleSequential();
      break;
    case RunMode::Parallel:
      FindMinParticleParallel();
      break;
    }
  }

  void FindMinParticleSequential() noexcept {
    T best_fitness = fitness_(gbest_);
    for (std::size_t i = 1; i < n_particles_; ++i) {
      const auto f = fitness_(particles_[i]);
      if (f < best_fitness) {
	best_fitness = f;
	gbest_ = particles_[i];
      }
    }
    
  }

  void FindMinParticleParallel() noexcept {

  }
  
  std::function<T(const std::vector<T>&)> fitness_;
  std::vector<Particle<T>> particles_;
  RunMode run_mode_;
  std::size_t n_particles_;
  std::size_t dim_;
  Particle<T> gbest_;
  
}; // class AoSExperiment

#endif // EXPERIMENT_H__
