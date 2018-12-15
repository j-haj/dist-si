use rand::{Rng,thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;

use std::ops::{Add,AddAssign};

///! Represents the position of a particle. A position does not necessarily
///! need to be a numeric value. A position object has a current coordinate,
///! represented by `coords`, a restriction on the min and max values allowed
///! for a respective coordinate, given by the respective element of `bounds`,
///! and a reflection property, given by `reflect`, which specifies whether a
///! particle should reflect when attempting to move beyond a boundary or
///! stay put. All particles must have a fitness function which is a mapping
///! f: Position -> f64
#[derive(Debug, Clone)]
pub struct Position<T: AddAssign + Add<Output = T> + Clone + SampleUniform> {
    coordinates: Vec<T>,
}

#[derive(Debug, Clone)]
pub struct Particle<T, F>
    where T: AddAssign + Add<Output = T> + Clone + SampleUniform, F: Fn(&Position<T>) -> f64 {
    position: Position<T>,
    pbest_pos: Position<T>,
    velocity: Velocity<T>,
    fitness: F,
}

///! Represents the velocity of a particle. Each coordinate component of the
///! particle has its own speed and direction -- collectively these make up
///! the particle's velocity. All velocities are clamped.
#[derive(Debug, Clone)]
pub struct Velocity<T> where T: AddAssign + Add<Output = T> + Clone + SampleUniform {
    velocities: Vec<T>,
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

impl<T: AddAssign + Add<Output = T> + Clone + SampleUniform> Position<T> {
    pub fn new(bounds: &Vec<(T,T)>) -> Position<T> {
        let mut rng = thread_rng();
        let mut initial_pos = Vec::new();
        for bound in bounds.iter() {
            let p_dist = Uniform::new(&bound.0, &bound.1);
            initial_pos.push(rng.sample(p_dist));
        }
        Position { coordinates: initial_pos, }
    }

    pub fn get_coordinates(&self) -> &Vec<T> { &self.coordinates }
}


impl<T> Velocity<T>
    where T: AddAssign + Add<Output = T> + Clone + SampleUniform {
    pub fn new(bounds: &Vec<(T,T)>) -> Velocity<T> {
        let mut rng = thread_rng();
        let mut initial_velocities = Vec::new();
        for bound in bounds.iter() {
            let v_dist = Uniform::new(&bound.0, &bound.1);
            initial_velocities.push(rng.sample(v_dist));
        }
        Velocity { velocities: initial_velocities, }
    }

    pub fn from_vec(velocities: Vec<T>) -> Velocity<T> {
        Velocity { velocities: velocities, }
    }

    pub fn get_velocities(&self) -> &Vec<T> { &self.velocities }
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
        for (i, p) in self.position.coordinates.iter_mut().enumerate() {
            *p += v.velocities[i].clone();
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
