use rand::{Rng, thread_rng};
use rand::distributions::Uniform;

///! Represents the velocity of a particle. Each coordinate component of the
///! particle has its own speed and direction -- collectively these make up
///! the particle's velocity. All velocities are clamped.
#[derive(Debug, Clone)]
pub struct Velocity {
    velocities: Vec<f64>,
}


impl Velocity {
    pub fn new(bounds: &[(f64,f64)]) -> Velocity {
        let mut rng = thread_rng();
        let mut initial_velocities = Vec::new();
        for bound in bounds.iter() {
            let v_dist = Uniform::new(&bound.0, &bound.1);
            initial_velocities.push(rng.sample(v_dist));
        }
        Velocity { velocities: initial_velocities, }
    }

    pub fn from_vec(velocities: Vec<f64>) -> Velocity {
        Velocity { velocities: velocities.to_vec(), }
    }

    pub fn velocities(&self) -> &[f64] { &self.velocities }
}
