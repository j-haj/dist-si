use std::ops::{Add, AddAssign};

use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;


///! Represents the velocity of a particle. Each coordinate component of the
///! particle has its own speed and direction -- collectively these make up
///! the particle's velocity. All velocities are clamped.
#[derive(Debug, Clone)]
pub struct Velocity<T> where T: AddAssign + Add<Output = T> + Clone + SampleUniform {
    velocities: Vec<T>,
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
