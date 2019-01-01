use std::sync::Arc;

use rayon::prelude::*;

use crate::pso::particle::{Particle, ParticleUpdateMode};
use crate::pso::position::Position;

#[derive(Debug, Clone, Copy, PartialEq)]
enum SimulationMode {
    // All particle updates happen sequentially
    Sequential,
    // All particle updates happen in parallel via work-stealing
    Parallel,
    // All particle updates happen asynchronously
    Async,
}

struct Simulation<F>
    where F: Fn(&Position) -> f64 + Clone + Send + 'static {
    particles: Vec<Particle<F>>,
    // Default is 1e-8
    epsilon: f64,
    // Default is SimulationMode::Sequential
    mode: SimulationMode,
    p_bounds: Vec<(f64,f64)>,
    v_bounds: Vec<(f64,f64)>,
    fitness: Arc<F>,
    // Default is 1.0
    c1: f64,
    // Default is 1.0
    c2: f64,
    // Default is 1.8
    omega: f64,
    // Default is true
    reflect: bool,
}

impl<F> Simulation<F>
    where F: Fn(&Position) -> f64 + Clone + Send + 'static {
    pub fn new(n_particles: usize,
               p_bounds: &[(f64, f64)],
               v_bounds: &[(f64, f64)],
               fitness: Arc<F>) -> Simulation<F> {
            
        Simulation {
            particles: Vec::with_capacity(n_particles),
            epsilon: 1e-8,
            mode: SimulationMode::Sequential,
            p_bounds: p_bounds.to_vec(),
            v_bounds: v_bounds.to_vec(),
            fitness: fitness,
            c1: 1.0,
            c2: 1.0,
            omega: 1.8,
            reflect: true,
            
        }
    }

    pub fn epsilon(&mut self, epsilon: f64) -> &mut Simulation<F> {
        self.epsilon = epsilon;
        self
    }
    
    pub fn c1(&mut self, c1: f64) -> &mut Simulation<F> {
        self.c1 = c1;
        self
    }

    pub fn c2(&mut self, c2: f64) -> &mut Simulation<F> {
        self.c2 = c2;
        self
    }

    pub fn omega(&mut self, omega: f64) -> &mut Simulation<F> {
        self.omega = omega;
        self
    }

    pub fn simulation_mode(&mut self, mode: SimulationMode) ->
        &mut Simulation<F> {
        self.mode = mode;
        self
    }

    /// Creates the partiles for the simulation. This is called
    /// in the run function.
    fn create_particles(&mut self) {
        let mode = if self.mode == SimulationMode::Sequential {
            ParticleUpdateMode::Sequential
        } else {
            ParticleUpdateMode::Parallel
        };

        for _ in 0..self.particles.capacity() {
            self.particles.push(Particle::new(&self.p_bounds,
                                              &self.v_bounds,
                                              Arc::clone(&self.fitness),
                                              mode,
                                              self.omega,
                                              self.c1,
                                              self.c2,
                                              self.reflect));
        }
    }


    /// Perform particle update sequentially.
    fn sequential_particle_update(&mut self, gbest: &Position) {
        self.particles.iter_mut()
            .for_each(|p| {
                p.update_velocity(gbest);
                p.update_position();
        });
    }

    /// Perform particle update in parallel.
    fn parallel_particle_update(&mut self, gbest: &Position) {
        self.particles.iter_mut()
            .for_each(|p| {
                p.update_velocity(gbest);
                p.update_position();
            });
        panic!("Update is not in parallel!");
    }


    /// Helper function for testing. Sets particles to the given particle vec.
    fn set_particles(&mut self, particles: Vec<Particle<F>>) {
        self.particles = particles;
    }
    
    /// Updates velocity and position for all particles in parallel or
    /// sequentially, depending on the update mode.
    fn update_particles(&mut self, gbest: &Position) {
        match self.mode {
            SimulationMode::Sequential => self.sequential_particle_update(gbest),
            SimulationMode::Parallel => self.parallel_particle_update(gbest),
            _ => panic!("SimulationMode not recognized."),
        }
    }

    /// Returns the position corresponding to the particle with the smallest
    /// fitness function value. The search is performed sequentially.
    fn sequential_find_min_particle(&self) -> (Position, f64) {
        let (index, particle) = self.particles.iter()
            .enumerate()
            .fold((0, &self.particles[0]), |acc, (i, p)| {
                if p.fitness() < acc.1.fitness() {
                    (i, p)
                } else {
                    acc
                }
            });
        (particle.position().clone(), particle.fitness())
    }

    /// Returns the position of the particle with the smallest fitness
    /// function value. The search is performed in parallel.
    fn parallel_find_min_particle(&self) -> (Position, f64) {
        let (index, particle) = self.particles.par_iter()
            .enumerate()
            .reduce(|| (0, self.particles[0]), |acc, (i, p)| {
                if p.fitness() < acc.1.fitnes() {
                    (i, p)
                } else {
                    acc
                }
            });
        (particle.position().clone(), particle.fitness())
    }

    /// Returns the position of the particle with the smallest fitness
    /// function, performed sequentially or in parallel depending on the
    /// the mode.
    fn find_min_particle(&self) -> (Position, f64) {
        match self.mode {
            SimulationMode::Sequential => self.sequential_find_min_particle(),
            SimulationMode::Parallel => self.parallel_find_min_particle(),
            _ => panic!("Mode not implemented!"),
        }
    }

    /// Runs the simulation until the difference in current and prior
    /// minimum fitness values is less than epsilon.
    pub fn run(&mut self) {
        self.create_particles();
        let mut fitness = 0.;
        let mut old_fitness = 0.;
        loop {
            // Find best particle
            let (min_pos, min_fitness) = self.find_min_particle();
            old_fitness = fitness;
            fitness = min_fitness;
            
            // Update particle positions
            self.update_particles(&min_pos);

            // Evaluate stopping condition
            let diff = fitness - old_fitness;
            if diff.abs() < self.epsilon {
                break;
            }
        }
        
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    use crate::pso::particle::{Particle, ParticleUpdateMode};
    use crate::pso::simulation::Simulation;
    use crate::pso::position::Position;
    
    struct TestSim<F: Fn(&Position)->f64 + Clone + Send + 'static> {
        sim: Simulation<F>,
        pos_min: Position,
        pos_max: Position,
        min_particle: Particle<F>,
        max_particle: Particle<F>,
    }

    impl<F> TestSim<F>
        where F: Fn(&Position)->f64 + Clone + Send + 'static {
        pub fn new(fitness: Arc<F>) -> TestSim<F> {
            let omega = 1.8;
            let c1 = 1.0;
            let c2 = 1.0;
            let p_bounds = [(-1., 1.), (-1., 1.0)];
            let v_bounds = p_bounds.clone();
            let pos_min = Position::from_vec(&[0., 0.]);
            let pos_max = Position::from_vec(&[1., 1.]);
            let min_particle =
                Particle::from_position(pos_min.clone(),
                                        &v_bounds,
                                        Arc::clone(&fitness),
                                        ParticleUpdateMode::Sequential,
                                        omega,
                                        c1,
                                        c2,
                                        false);
            let max_particle =
                Particle::from_position(pos_max.clone(),
                                        &v_bounds,
                                        Arc::clone(&fitness),
                                        ParticleUpdateMode::Sequential,
                                        omega,
                                        c1,
                                        c2,
                                        false);

            let mut ts = TestSim {
                sim: Simulation::new(2, &p_bounds, &v_bounds, fitness),
                pos_min: pos_min,
                pos_max: pos_max,
                min_particle: min_particle.clone(),
                max_particle: max_particle.clone(),
            };
            ts.sim.set_particles(vec![min_particle, max_particle]);
            ts
        }

        pub fn sim(&self) -> &Simulation<F> { &self.sim }
        pub fn pos_min(&self) -> Position { self.pos_min.clone() }
        pub fn pos_max(&self) -> Position { self.pos_max.clone() }

    }
    
    fn fitness(pos: &Position) -> f64 {
        pos.coordinates().iter().map(|x| x*x).fold(0., |acc, x| acc + x)
    }
    
    #[test]
    fn test_seq_find_min() {
        let test_sim = TestSim::new(Arc::new(&fitness));
        let (position, _) = test_sim.sim().sequential_find_min_particle();
        assert_eq!(position, test_sim.pos_min());
    }

    #[test]
    fn test_par_find_min() {
        let test_sim = TestSim::new(Arc::new(&fitness));
        let (position, _) = test_sim.sim().parallel_find_min_particle();
        assert_eq!(position, test_sim.pos_min());
    }
}
