use ndarray;
use ndarray::Array1;

///! Represents the velocity of a particle. Each coordinate component of the
///! particle has its own speed and direction -- collectively these make up
///! the particle's velocity. All velocities are clamped.
#[derive(Debug, Clone, PartialEq)]
pub struct Velocity {
    components: Array1<f64>,
}


impl Velocity {
    pub fn new(dim: usize) -> Velocity {
        Velocity { components: Array1::zeros(dim) }
    }


    pub fn from_vec(velocities: Vec<f64>) -> Velocity {
        Velocity { components: Array1::from(velocities.to_vec()), }
    }

    pub fn components(&self) -> &Array1<f64> { &self.components }
    pub fn components_mut(&mut self) -> &mut Array1<f64> {
        &mut self.components
    }
}
