use ndarray;
use ndarray::{Array1, Zip};
use ndarray_parallel::prelude::*;
use rand::{Rng, thread_rng};
use rand::distributions::{Distribution, Uniform};
use rayon::prelude::*;

use crate::pso::position::Position;
use crate::pso::velocity::Velocity;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParticleUpdateMode {
    Sequential,
    Parallel,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Particle<F>
    where F: Fn(&Position) -> f64 {
    position: Position,
    pbest_pos: Position,
    velocity: Velocity,
    max_speed: Array1<f64>,
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
        let max_speeds = Array1::from_vec(v_bounds.clone()
                                          .into_iter()
                                          .map(|&x| x.1)
                                          .collect());
        Particle { position: initial_pos.clone(),
                   pbest_pos: initial_pos,
                   velocity: Velocity::new(v_bounds.len()),
                   max_speed: max_speeds,
                   fitness: f,
                   mode: mode,
                   pbest: std::f64::INFINITY,
                   omega: omega,
                   c1: c1,
                   c2: c2,
                   reflect: reflect,
        }

    }

    pub fn from_position(position: Position,
                         v_bounds: &[(f64,f64)],
                         f: F,
                         mode: ParticleUpdateMode,
                         omega: f64,
                         c1: f64,
                         c2: f64,
                         reflect: bool) -> Particle<F> {
        let max_speeds = Array1::from_vec(v_bounds.clone()
                                          .into_iter()
                                          .map(|&x| x.1)
                                          .collect());
        Particle {
            position: position.clone(),
            pbest_pos: position,
            velocity: Velocity::new(v_bounds.len()),
            max_speed: max_speeds,
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
        let local_scale = self.c1 * r1;
        let global_scale = self.c2 * r2;
        let omega = self.omega;

        self.velocity.components_mut().map_inplace(|v| *v *= omega);
        self.velocity.components_mut().scaled_add(local_scale,
                                                  self.pbest_pos.coordinates());
        self.velocity.components_mut().scaled_add(-1. * local_scale,
                                                  self.position.coordinates());
        self.velocity.components_mut().scaled_add(global_scale,
                                                  gbest.coordinates());
        self.velocity.components_mut().scaled_add(-1. * global_scale,
                                                  self.position.coordinates());
        // Check speeds
        self.velocity.components_mut()
            .zip_mut_with(&self.max_speed, |v, &v_max| *v = v.min(v_max));
            
    }

    /// Updates particle position and returns the new fitness value.
    pub fn update_position(&mut self) -> f64 {
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
            println!("Updated pbest to {:?}", f);
            self.pbest = f;
            self.pbest_pos = self.position.clone();
        }
        f
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

    pub fn velocity(&self) -> &Velocity { &self.velocity }

    pub fn fitness(&mut self) -> f64 { (self.fitness)(&self.position) }

    pub fn pbest(&self) -> f64 { self.pbest }
}

