use ndarray::Array1;
use rand::{Rng, thread_rng};
use rand::distributions::Uniform;

///! Represents the position of a particle. A position does not necessarily
///! need to be a numeric value. A position object has a current coordinate,
///! represented by `coords`, a restriction on the min and max values allowed
///! for a respective coordinate, given by the respective element of `bounds`,
///! and a reflection property, given by `reflect`, which specifies whether a
///! particle should reflect when attempting to move beyond a boundary or
///! stay put. All particles must have a fitness function which is a mapping
///! f: Position -> f64
#[derive(Debug, Clone, PartialEq)]
pub struct Position {
    coordinates: Array1<f64>,
}

impl Position {
    pub fn new(bounds: &[(f64,f64)]) -> Position {
        let mut rng = thread_rng();
        let mut initial_pos = Vec::new();
        for bound in bounds.iter() {
            let p_dist = Uniform::new(&bound.0, &bound.1);
            initial_pos.push(rng.sample(p_dist));
        }
        Position { coordinates: Array1::from(initial_pos), }
    }

    pub fn from_vec(coordinates: &[f64]) -> Position {
        Position { coordinates: Array1::from(coordinates.to_vec()) }
    }

    pub fn coordinates(&self) -> &Array1<f64> { &self.coordinates }

    pub fn coordinates_mut(&mut self) -> &mut Array1<f64> { &mut self.coordinates }
}



