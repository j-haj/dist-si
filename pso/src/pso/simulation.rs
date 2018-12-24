use use pso::pso::particle::Particle;

enum SimulationMode {
    // All particle updates happen sequentially
    Sequential,
    // All particle updates happen in parallel via work-stealing
    Parallel,
    // All particle updates happen asynchronously
    Async,
}

struct Simulation {
    particles: Vec<Particle>,
    epsilon: f64,
    mode: SimulationMode,
}

impl Simulation {
    pub fn new(n_particles: u64,
               epsilon: f64,
               pos_bounds: &[(f64, f64)],
               v_bounds: &[(f64, f64)],
               f: F,
               omega: f64,
               c1: f64,
               c2: f64,
               simulation_mode: SimulationMode,
               reflect: bool) -> Simulation {
        let mode = if simulation_mode == SimulationMode::Sequential {
            ParticleUpdateMode::Parallel
        } else {
            ParticleUpdateMode::Sequential
        }
            
        // Create particles
        let mut particles = Vec::new();
        for _ in 0..n_particles {
            particles.push(Particle::new(pos_bounds,
                                         v_bounds,
                                         f,
                                         mode,
                                         omega,
                                         c1,
                                         c2,
                                         reflect));
        }
                                         
        Simulation {
            particles: particles, epsilon: epsilon, mode: mode
        }
    }


    fn sequential_particle_update(&mut self, gbest: &Position) {
        self.particles.iter_mut()
            .for_each(|p| {
                p.update_velocity();
                p.update_position();
        });
    }

    fn parallel_particle_update(&mut self, gbest: &Position) {
        self.particles.par_iter_mut()
            .for_each(|p| {
                p.update_velocity(gbest);
                p.update_position();
            });
    }

    fn update_particles(&mut self, gbest: &Position) {
        match self.mode {
            SimulationMode::Sequential => sequential_particle_update(gbest),
            SimulationMode::Parallel => parallel_particle_update(gbest),
            _ => panic!("SimulationMode not recognized."),
        }
    }

    fn sequential_find_min_particle(&self) -> Position {
        (_, index) = self.particles.iter()
            .enumerate()
            .fold(|acc, (p, i)| if p.fitness() < acc.0.fitness() {
                (p, i)
            } else {
                acc
            });
        self.particles[i].position().clone()
    }

    fn parallel_find_min_particle(&self) -> Position {
        panic!("Not implemented.");
    }

    fn find_min_particle(&self) -> Position {
        match self.mode {
            SimulationMode::Sequential => self.sequential_find_min_particle(),
            SimulationMode::Parallel => self.parallel_find_min_particle(),
        }
    }
    
    pub fn run(&self) {
        let mut fitness: f64;
        let mut old_fitness: f64 = std::f64::INFINITY;
        loop {
            // Find best particle
            let min_pos = self.find_min_particle();
            
            // Update particle positions
            let fitness = self.update_particles(&min_pos);

            // Evaluate stopping condition
            if (fitness - old_fitness).abs() < self.epsilon {
                break;
            }
        }
        
    }
}
