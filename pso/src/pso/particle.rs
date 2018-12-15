use std::ops::{Add, AddAssign};

use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;

use crate::pso::position::Position;
use crate::pso::velocity::Velocity;

#[derive(Debug, Clone)]
pub struct Particle<T, F>
    where T: AddAssign + Add<Output = T> + Clone + SampleUniform, F: Fn(&Position<T>) -> f64 {
    position: Position<T>,
    pbest_pos: Position<T>,
    velocity: Velocity<T>,
    fitness: F,
}


#[derive(Debug)]
pub struct ParticleUpdater<T>
    where T: AddAssign + Add<Output = T> + Clone + SampleUniform {
    omega: f64,
    c1: f64,
    c2: f64,
    position_bounds: Vec<(T,T)>,
    velocity_bounds: Vec<(T,T)>,
    reflect: bool,
}

impl<T, F> Particle<T, F>
    where T: AddAssign + Add<Output = T> + Clone + SampleUniform, F: Fn(&Position<T>) -> f64 {
    pub fn new(pos_bounds: &Vec<(T,T)>, v_bounds: &Vec<(T, T)>, f: F)
               -> Particle<T, F> {
        let initial_pos = Position::new(pos_bounds);
        Particle { position: initial_pos.clone(),
                   pbest_pos: initial_pos,
                   velocity: Velocity::new(v_bounds),
                   fitness: f, }

    }

    pub fn update_position(&mut self, v: &Velocity<T>) {
        for (i, p) in self.position.coordinates_mut().iter_mut().enumerate() {
            *p += v.get_velocities()[i].clone();
        }
    }

    pub fn get_position(&self) -> &Position<T> { &self.position }

    pub fn get_fitness(&self) -> f64 { (self.fitness)(&self.position) }
}


impl<T: AddAssign + Add<Output = T> + Clone + SampleUniform> ParticleUpdater<T> {
    pub fn new(omega: f64, c1: f64, c2: f64, position_bounds: Vec<(T,T)>,
               velocity_bounds: Vec<(T,T)>, reflect: bool)
               -> ParticleUpdater<T> {
        ParticleUpdater { omega: omega, c1: c1, c2: c2,
                          position_bounds: position_bounds,
                          velocity_bounds: velocity_bounds,
                          reflect: reflect }
    }

    fn new_velocity<F>(&self, particle: &Particle<T, F>, gbest: &Particle<T, F>)
        -> Velocity<T>
        where F: Fn(&Position<T>) -> f64 {
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

