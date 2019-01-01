#ifndef EXPERIMENT_H__
#define EXPERIMENT_H__

#include <cctype>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "particle.h"
#include "particles.h"

enum class RunMode {
  Sequential,
  Parallel
};


template <typename T>
class AoSExperiment {

public:
  AoSExperiment(std::function<T(const std::vector<T>&)> f,
		RunMode run_mode,
		T x_min, T x_max, T v_min, T v_max,
		std::size_t n_particles, std::size_t dim,
		T omega, T c1, T c2)
    : particles_(std::vector<Particle<T>>(n_particles)),
      run_mode_(run_mode),
      n_particles_(n_particles),
      dim_(dim),
      gbest_(Particle<T>(x_min, x_max, v_min, v_max,
			 dim, omega, c1, c2, f)) {

    // Generate particles for simulation
    for (std::size_t i = 0; i < n_particles; ++i) {
      particles_[i] = Particle<T>(x_min, x_max, v_min, v_max,
				  dim, omega, c1, c2, f);
    }

    gbest_ = particles_[0];
    FindMinParticle();
  }

  void Run(std::size_t max_steps) {
    std::size_t n_steps = 0;

    std::cout << "Beginning simulation with "
	      << n_particles_
	      << " particles\n";
    auto start = std::chrono::steady_clock::now();
    
    while (n_steps < max_steps) {
      // Get global min
      FindMinParticle();
      
      // Update particle velocities
      UpdateParticleVelocities();

      // Update particle positions
      UpdateParticlePositions();

      ++n_steps;
    }
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    std::cout << "Simulation done after " << n_steps << " steps ";
    std::cout << "in " << elapsed.count() << " seconds\n";
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
    T best_fitness = gbest_.Fitness();
    for (const auto& p : particles_) {
      const auto f = p.Fitness();
      if (f < best_fitness) {
	best_fitness = f;
	gbest_ = p;
      }
    }
    
  }

  void FindMinParticleParallel() noexcept {
    T best_fitness = gbest_.Fitness();
#pragma omp parallel for
    for (std::size_t i = 0; i < particles_.size(); ++i) {
      const auto f = particles_[i].Fitness();
      if (f < best_fitness) {
#pragma omp critical
	if (f < best_fitness) {
	  best_fitness = f;
	  gbest_ = particles_[i];
	}
      }
    }
  }

  void UpdateParticlePositionsSequential() noexcept {
    for (auto& p : particles_) {
      p.UpdatePosition();
    }
  }

  void UpdateParticlePositionsParallel() noexcept {
#pragma omp parallel for
    for (std::size_t i = 0; i < particles_.size(); ++i) {
      particles_[i].UpdatePosition();
    }
  }

  void UpdateParticleVelocitiesSequential() noexcept {
    for (auto& p : particles_) {
      p.UpdateVelocity(gbest_);
    }
  }

  void UpdateParticleVelocitiesParallel() noexcept {
#pragma omp parallel for
    for (std::size_t i = 0; i < particles_.size(); ++i) {
      particles_[i].UpdateVelocity(gbest_);
    }
  }
  
  std::vector<Particle<T>> particles_;
  RunMode run_mode_;
  std::size_t n_particles_;
  std::size_t dim_;
  Particle<T> gbest_;
  
}; // class AoSExperiment

template <typename T>
class SoAExperiment {
 public:
  SoAExperiment(std::function<T(const std::vector<T>&)> f,
		ParticleUpdateMode p_mode,
		T x_min, T x_max, T v_min, T v_max,
		std::size_t n_particles, std::size_t dim,
		T omega, T c1, T c2)
    : particles_(Particles<T>(f, p_mode, n_particles, dim, x_min, x_max,
			      v_min, v_max, omega, c1, c2)),
      dim_(dim) {}

  void Run(std::size_t max_steps) {
    std::size_t n_steps = 0;

    std::cout << "Beginning simulation with "
	      << particles_.size()
	      << " particles\n";
    auto start = std::chrono::steady_clock::now();
    
    while (n_steps < max_steps) {
      // Update particle velocities
      particles_.UpdateVelocities();

      // Update particle positions
      particles_.UpdatePositions();

      ++n_steps;
    }
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    std::cout << "Simulation done after " << n_steps << " steps ";
    std::cout << "in " << elapsed.count() << " seconds\n";
    std::cout << "\tgbest: " << particles_.gbest() << '\n';

  }
   
 private:
  Particles<T> particles_;
  std::size_t dim_;
  
}; // class SoAExperiment

#endif // EXPERIMENT_H__
