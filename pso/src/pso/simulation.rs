use pso::pso::particle::Particle;

enum SimulationMode {
    // All particle updates happen sequentially
    Sequential,
    // All particle updates happen in parallel via work-stealing
    Parallel,
    // All particle updates happen asynchronously
    Async,
}

struct Simulation<F>
    where F: Fn(&Position) -> f64 {
    particles: Vec<Particle>,
    // Default is 1e-8
    epsilon: f64,
    // Default is SimulationMode::Sequential
    mode: SimulationMode,
    p_bounds: Vec<[(f64,f64)]>,
    v_bounds: Vec<[(f64,f64)]>,
    fitness: F,
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
    where F: Fn(&Position) -> f64 {
    pub fn new(n_particles: u64,
               pos_bounds: &[(f64, f64)],
               v_bounds: &[(f64, f64)],
               fitness: F) -> Simulation<F> {
        let mode = if simulation_mode == SimulationMode::Sequential {
            ParticleUpdateMode::Sequential
        } else if simulation_mode == SimulationMode::Parallel {
            ParticleUpdateMode::Parallel
        } else {
            panic!("{} is not implemented!", simulation_mode);
        }
            
        Simulation {
            particles: Vec::with_capacity(n_particles),
            epsilon: 1e-8,
            mode: SimulationMode::Sequential,
            p_bounds: pos_bounds,
            v_bounds: v_bounds,
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

    /// Creates the partiles for the simulation. This is called
    /// in the run function.
    fn create_particles(&mut self) {
        for _ in 0..self.particles.capacity() {
            self.particles.push(Particle::new(self.pos_bounds,
                                              self.v_bounds,
                                              self.f,
                                              self.mode,
                                              self.c1,
                                              self.c2,
                                              self.reflect));
        }
    }


    /// Perform particle update sequentially.
    fn sequential_particle_update(&mut self, gbest: &Position) {
        self.particles.iter_mut()
            .for_each(|p| {
                p.update_velocity();
                p.update_position();
        });
    }

    /// Perform particle update in parallel.
    fn parallel_particle_update(&mut self, gbest: &Position) {
        self.particles.par_iter_mut()
            .for_each(|p| {
                p.update_velocity(gbest);
                p.update_position();
            });
    }


    /// Helper function for testing. Sets particles to the given particle vec.
    fn set_particles(&mut self, particles: Vec<Particle>) {
        self.particles = particles;
    }
    
    /// Updates velocity and position for all particles in parallel or
    /// sequentially, depending on the update mode.
    fn update_particles(&mut self, gbest: &Position) {
        match self.mode {
            SimulationMode::Sequential => sequential_particle_update(gbest),
            SimulationMode::Parallel => parallel_particle_update(gbest),
            _ => panic!("SimulationMode not recognized."),
        }
    }

    /// Returns the position corresponding to the particle with the smallest
    /// fitness function value. The search is performed sequentially.
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

    /// Returns the position of the particle with the smallest fitness
    /// function value. The search is performed in parallel.
    fn parallel_find_min_particle(&self) -> Position {
        panic!("Not implemented.");
    }

    /// Returns the position of the particle with the smallest fitness
    /// function, performed sequentially or in parallel depending on the
    /// the mode.
    fn find_min_particle(&self) -> Position {
        match self.mode {
            SimulationMode::Sequential => self.sequential_find_min_particle(),
            SimulationMode::Parallel => self.parallel_find_min_particle(),
        }
    }

    /// Runs the simulation until the difference in current and prior
    /// minimum fitness values is less than epsilon.
    pub fn run(&self) {
        self.create_particles();
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

#[cfg(test)]
mod tests {

    static omega: f64 = 1.8;
    static c1: f64 = 1.0;
    static c2: f64 = 1.0;
    static p_bounds: Vec<(f64,f64)> = vec![(-1.,1.), (-1.,1.)];
    static v_bounds: Vec<(f64,f64)> = p_bounds.clone();

    static pos_min: Position = Position::from_vec(vec![0., 0.]);
    static pos_max: Position = Position::from_vec(vec![1., 1.]);
    
    static p_min: Particle<Fn(&Position)->f64> =
        Particle::from_position(pos_min,
                                &v_bounds,
                                fitness,
                                ParticleUpdateMode::Sequential,
                                omega,
                                c1,
                                c2,
                                false);
    static p_max: Particle<Fn(&Position)->f64> =
        Particle::from_position(pos_max,
                                &v_bounds,
                                fitness,
                                ParticleUpdateMode::Sequential,
                                omega,
                                c1,
                                c2,
                                false);
    
    fn fitness(pos: &Position) -> f64 {
        p.coordinates().iter().map(|x| x*x).fold(0., |acc, x| acc + x)
    }
    
    #[test]
    fn test_seq_find_min() {
        let sim = Simulation::new(2, &p_bounds, &v_bounds, fitness);
        sim.set_particles(vec![p_min, p_max]);
        let actual = sim.sequential_find_min();
        assert_eq!(actual.position(), pos_min);
    }

    #[test]
    fn test_par_find_min() {
        let sim = Simulation::new(2, &p_pounds, &v_bounds, fitness);
        sim.set_particles(vec![p_min, p_max]);
        let actual = sim.parallel_find_min();
        assert_eq!(actual.position(), pos_min);
    }
}
