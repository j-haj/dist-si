use rand::{Rng,thread_rng};
use rand::distributions::Uniform;
use rand::distributions::uniform::SampleUniform;

///! Represents the position of a particle. A position does not necessarily
///! need to be a numeric value. A position object has a current coordinate,
///! represented by `coords`, a restriction on the min and max values allowed
///! for a respective coordinate, given by the respective element of `bounds`,
///! and a reflection property, given by `reflect`, which specifies whether a
///! particle should reflect when attempting to move beyond a boundary or
///! stay put.
#[derive(Debug, Clone)]
pub struct Position<T: Clone + SampleUniform> {
    coords: Vec<T>,
    bounds: Vec<(T,T)>,
    reflect: bool,
}

impl<T: Clone + SampleUniform> Position<T> {
    pub fn new(bounds: Vec<(T,T)>, reflect: bool) -> Position<T> {
        let mut rng = thread_rng();
        let mut initial_pos = Vec::new();
        for bound in bounds.iter() {
            let p_dist = Uniform::new(&bound.0, &bound.1);
            initial_pos.push(rng.sample(p_dist));
        }
        Position { coords: initial_pos, bounds: bounds, reflect: reflect }
    }
}

///! Represents the velocity of a particle. Each coordinate component of the
///! particle has its own speed and direction -- collectively these make up
///! the particle's velocity. All velocities are clamped.
#[derive(Debug, Clone)]
pub struct Velocity {
    velocities: Vec<f64>,
    bounds: Vec<(f64,f64)>,
}

impl Velocity {
    pub fn new(bounds: Vec<(f64,f64)>) -> Velocity{
        let mut rng = thread_rng();
        let mut initial_velocities = Vec::new();
        for bound in bounds.iter() {
            let v_dist = Uniform::new(&bound.0, &bound.1);
            initial_velocities.push(rng.sample(v_dist));
        }
        Velocity { velocities: initial_velocities, bounds: bounds }
    }
}

#[derive(Debug, Clone)]
pub struct Particle<T: Clone + SampleUniform> {
    position: Position<T>,
    pbest_pos: Position<T>,
    velocity: Velocity,
}

impl<T: Clone + SampleUniform> Particle<T> {
    pub fn new(pos_bounds: Vec<(T,T)>, v_bounds: Vec<(f64, f64)>,
               reflect: bool) -> Particle<T> {
        let initial_pos = Position::new(pos_bounds, reflect);
        Particle { position: initial_pos.clone(),
                   pbest_pos: initial_pos,
                   velocity: Velocity::new(v_bounds) }

    }
    
}

