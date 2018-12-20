use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use rayon::prelude::*;

use crate::pso::position::Position;
use crate::pso::velocity::Velocity;

#[derive(Debug, Clone, Copy)]
pub enum ParticleUpdateMode {
    Sequential,
    Parallel,
}

#[derive(Debug, Clone)]
pub struct Particle<F>
    where F: Fn(&Position) -> f64 {
    position: Position,
    pbest_pos: Position,
    velocity: Velocity,
    fitness: F,
    mode: ParticleUpdateMode,
}


#[derive(Debug)]
pub struct ParticleUpdater {
    omega: f64,
    c1: f64,
    c2: f64,
    position_bounds: Vec<(f64,f64)>,
    velocity_bounds: Vec<(f64,f64)>,
    reflect: bool,
}

impl<F> Particle<F>
    where F: Fn(&Position) -> f64 {
    pub fn new(pos_bounds: &[(f64,f64)], v_bounds: &[(f64,f64)],
               f: F, mode: ParticleUpdateMode) -> Particle<F> {
        let initial_pos = Position::new(pos_bounds);
        Particle { position: initial_pos.clone(),
                   pbest_pos: initial_pos,
                   velocity: Velocity::new(v_bounds),
                   fitness: f,
                   mode: mode, }

    }

    pub fn update_position(&mut self, velocity: &Velocity) {
        match self.mode {
            ParticleUpdateMode::Sequential =>
                self.sequential_position_update(velocity),
            ParticleUpdateMode::Parallel =>
                self.data_parallel_position_update(velocity),
            _ => panic!("invalid mode - expected 'Sequential', \
                        or 'DataParallel'."),
        }
    }

    fn sequential_position_update(&mut self, velocity: &Velocity) {
        for (p, v) in self.position.coordinates_mut()
            .iter_mut()
            .zip(velocity.velocities()) {
            *p += *v;
        }
    }

    fn data_parallel_position_update(&mut self, velocity: &Velocity) {
        self.position.coordinates_mut()
            .par_iter_mut()
            .zip(velocity.velocities())
            .for_each(|(p, v)| *p += *v)
    }

    pub fn get_position(&self) -> &Position { &self.position }

    pub fn get_fitness(&self) -> f64 { (self.fitness)(&self.position) }
}


impl ParticleUpdater {
    pub fn new(omega: f64, c1: f64, c2: f64, position_bounds: Vec<(f64,f64)>,
               velocity_bounds: Vec<(f64,f64)>, reflect: bool) -> ParticleUpdater {
        ParticleUpdater { omega: omega, c1: c1, c2: c2,
                          position_bounds: position_bounds,
                          velocity_bounds: velocity_bounds,
                          reflect: reflect }
    }

    fn new_velocity<F>(&self, particle: &Particle<F>, gbest: &Particle<F>)
        -> Velocity
        where F: Fn(&Position) -> f64 {
        // Calculate r1 and r2
        let mut rng = thread_rng();
        let u_dist = Uniform::new(0.0, 1.0);
        let r1 = rng.sample(u_dist);
        let r2 = rng.sample(u_dist);

        // Get global best
        let gbest = self.global_best();

        // Get local best
        //let pbest = particle.pbest;

        // Crate new velocity
        let vs = Vec::new();
        Velocity::from_vec(vs)
    }

    fn global_best(&self) -> usize { 0 }
}

