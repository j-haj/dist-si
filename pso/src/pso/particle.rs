use ndarray;
use ndarray::{Array1, Zip};
use ndarray_parallel::prelude::*;
use rand::{Rng, thread_rng};
use rand::distributions::{Distribution, Uniform};
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
    pbest: f64,
    omega: f64,
    c1: f64,
    c2: f64,
    reflect: bool,
}


#[derive(Debug)]
pub struct ParticleUpdater {
    position_bounds: Vec<(f64,f64)>,
    velocity_bounds: Vec<(f64,f64)>,
}

impl<F> Particle<F>
    where F: Fn(&Position) -> f64 {
    pub fn new(pos_bounds: &[(f64,f64)],
               v_bounds: &[(f64,f64)],
               f: F,
               mode: ParticleUpdateMode,
               omega: f64,
               c1: f64,
               c2: f64,
               reflect: bool
    ) -> Particle<F> {
        let initial_pos = Position::new(pos_bounds);
        Particle { position: initial_pos.clone(),
                   pbest_pos: initial_pos,
                   velocity: Velocity::new(v_bounds),
                   fitness: f,
                   mode: mode,
                   pbest: std::f64::INFINITY,
                   omega: omega,
                   c1: c1,
                   c2: c2,
                   reflect: reflect,
        }

    }

    pub fn update_velocity(&mut self, gbest: &Position) {
        let mut rng = thread_rng();
        let u = Uniform::new(0., 1.);
        let r1 : f64 = u.sample(&mut rng);
        let r2 : f64 = u.sample(&mut rng);
        let omega = self.omega;
        self.velocity.components_mut().map_inplace(|v| *v = *v * omega);
        self.velocity.components_mut().scaled_add(self.c1 * r1,
                                                  self.pbest_pos.coordinates());
        self.velocity.components_mut().scaled_add(-1. * self.c1 * r1,
                                                  self.position.coordinates());
        self.velocity.components_mut().scaled_add(self.c2 * r2,
                                                  gbest.coordinates());
        self.velocity.components_mut().scaled_add(-1. * self.c2 * r2,
                                                  self.position.coordinates());
    }

    pub fn update_position(&mut self) {
        match self.mode {
            ParticleUpdateMode::Sequential =>
                self.sequential_position_update(),
            ParticleUpdateMode::Parallel =>
                self.data_parallel_position_update(),
            _ => panic!("invalid mode - expected 'Sequential', \
                        or 'DataParallel'."),
        }
        let f = self.fitness();
        if f < self.pbest {
            self.pbest = f;
            self.pbest_pos = self.position.clone();
        }
    }

    fn sequential_position_update(&mut self) {
        self.position.coordinates_mut()
            .zip_mut_with(self.velocity.components(),|p, &v| *p += v);
    }

    fn data_parallel_position_update(&mut self) {
        Zip::from(self.position.coordinates_mut())
            .and(self.velocity.components())
            .par_apply(|p, &v| *p += v)
    }

    pub fn position(&self) -> &Position { &self.position }

    pub fn fitness(&mut self) -> f64 { (self.fitness)(&self.position) }

    pub fn pbest(&self) -> f64 { self.pbest }
}


impl ParticleUpdater {
    pub fn new(omega: f64, c1: f64, c2: f64, position_bounds: Vec<(f64,f64)>,
               velocity_bounds: Vec<(f64,f64)>, reflect: bool) -> ParticleUpdater {
        ParticleUpdater { position_bounds: position_bounds,
                          velocity_bounds: velocity_bounds,
        }
    }

    fn global_best(&self) -> usize { 0 }
}

