#ifndef EXPERIMENT_H__
#define EXPERIMENT_H__

#include <cctype>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
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
      dim_(dim),
      gbest_(Particle<T>(x_min, x_max, v_min, v_max,
			 dim, omega, c1, c2)) {

    // Generate particles for simulation
    for (std::size_t i = 0; i < n_particles; ++i) {
      particles_[i] = Particle<T>(x_min, x_max, v_min, v_max,
				  dim, omega, c1, c2);
    }

    gbest_ = particles_[0];
    FindMinParticle();
  }

  void Run(std::size_t max_steps) {
    std::cout << "Beginning simulation with "
	      << n_particles_
	      << " particles\n";
    auto start = std::chrono::steady_clock::now();
    switch (run_mode_) {
    case RunMode::Sequential:
      std::cout << "RunMode: Sequential\n";
      RunSequential(max_steps);
      break;
    case RunMode::Parallel:
      std::cout << "RunMode: Parallel\n";
      RunParallel(max_steps);
      break;
    }
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    std::cout << "in " << elapsed.count() << " seconds\n";
    //    std::cout << "\tsteps taken: " << n_steps << '\n';
    std::cout << "\tgbest: " << gbest_ << '\n';
  }

  void UpdateParticleVelocities() noexcept {
    switch (run_mode_) {
    case RunMode::Sequential:
      UpdateParticleVelocitiesSequential();
      break;
    case RunMode::Parallel:
      UpdateParticleVelocitiesParallel();
      break;
    }
  }

  void UpdateParticlePositions() noexcept {
    switch (run_mode_) {
    case RunMode::Sequential:
      UpdateParticlePositionsSequential();
      break;
    case RunMode::Parallel:
      UpdateParticlePositionsParallel();
      break;
    }
  }
  
  
private:

  void RunSequential(std::size_t max_steps) {
    std::size_t n_steps = 0;
    T old_fitness = static_cast<T>(0);
    T current_fitness = 1;

    while (n_steps < max_steps) {
      // Get global min
      FindMinParticle();
      
      old_fitness = current_fitness;
      current_fitness = fitness_(gbest_);

      // Update particle velocities
      UpdateParticleVelocities();

      // Update particle positions
      UpdateParticlePositions();

      ++n_steps;
    }
    std::cout << "Simulation done after " << n_steps << " steps ";
  }

  void RunParallel(std::size_t max_steps) {
    
  }

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
    for (const auto& p : particles_) {
      const auto f = fitness_(p);
      if (f < best_fitness) {
	best_fitness = f;
	gbest_ = p;
      }
    }
    
  }

  void FindMinParticleParallel() noexcept {
    T best_fitness = fitness_(gbest_);
    T f;
#pragma omp for local f atomic best_fitness, gbest_
    for (std::size_t i = 0; i < particles_.size(); ++i) {
      f = fitness_(particles_[i]);
      if (f < best_fitness) {
	best_fitness = f;
	gbest_ = p;
    }
  }

  void UpdateParticlePositionsSequential() noexcept {
    for (auto& p : particles_) {
      p.UpdatePosition();
    }
  }

  void UpdateParticlePositionsParallel() noexcept {

  }

  void UpdateParticleVelocitiesSequential() noexcept {
    for (auto& p : particles_) {
      p.UpdateVelocity(gbest_);
    }
  }

  void UpdateParticleVelocitiesParallel() noexcept {

  }
  
  std::function<T(const Particle<T>&)> fitness_;
  std::vector<Particle<T>> particles_;
  RunMode run_mode_;
  std::size_t n_particles_;
  std::size_t dim_;
  Particle<T> gbest_;
  
}; // class AoSExperiment

#endif // EXPERIMENT_H__
