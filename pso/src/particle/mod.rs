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
}

#[derive(Debug, Clone)]
pub struct Particle<T: Clone + SampleUniform> {
    position: Position<T>,
    pbest_pos: Position<T>,
    velocity: Velocity,
}

///! Represents the velocity of a particle. Each coordinate component of the
///! particle has its own speed and direction -- collectively these make up
///! the particle's velocity. All velocities are clamped.
#[derive(Debug, Clone)]
pub struct Velocity {
    velocities: Vec<f64>,
}

#[derive(Debug)]
pub struct ParticleUpdater {
    omega: f64,
    c1: f64,
    c2: f64,
    postion_bounds: Vec<(T,T)>,
    velocity_bounds: Vec<(f64,f64)>,
    reflect: bool,
}

impl<T: Clone + SampleUniform> Position<T> {
    pub fn new(bounds: &Vec<(T,T)>) -> Position<T> {
        let mut rng = thread_rng();
        let mut initial_pos = Vec::new();
        for bound in bounds.iter() {
            let p_dist = Uniform::new(&bound.0, &bound.1);
            initial_pos.push(rng.sample(p_dist));
        }
        Position { coords: initial_pos, }
    }
}


impl Velocity {
    pub fn new(bounds: &Vec<(f64,f64)>) -> Velocity{
        let mut rng = thread_rng();
        let mut initial_velocities = Vec::new();
        for bound in bounds.iter() {
            let v_dist = Uniform::new(&bound.0, &bound.1);
            initial_velocities.push(rng.sample(v_dist));
        }
        Velocity { velocities: initial_velocities, }
    }

    pub fn from_vec(velocities: Vec<f64>) -> Velocity {
        Velocity { velocities: velocities, }
    }
}


impl<T: Clone + SampleUniform> Particle<T> {
    pub fn new(pos_bounds: &Vec<(T,T)>, v_bounds: &Vec<(f64, f64)>)
               -> Particle<T> {
        let initial_pos = Position::new(pos_bounds);
        Particle { position: initial_pos.clone(),
                   pbest_pos: initial_pos,
                   velocity: Velocity::new(v_bounds) }

    }

    pub fn update_position(&self, v: &Velocity) {
        for (p, i) in &self.coords.iter().enumerate() {
            *p += *p + v[i];
        }
    }
}

impl<T: SampleUniform> ParticleUpdater {
    pub fn new(omega: f64, c1: f64, c2: f64, position_bounds: Vec<(T,T)>,
               velocity_bounds: Vec<(f64,f64)>, reflect: bool) -> ParticleUpdater {
        ParticleUpdater { omega: omega, c1: c1, c2: c2,
                          position_bounds: position_bounds,
                          velocity_bounds: velocity_bounds,
                          reflect: reflect }
    }

    pub fn update(&self, particle: mut &Particle<T>, gbest: &Particle<T>) {
        
    }


    fn new_velocity(&self, particle: &Particle<T>, gbest: &Particle<T>)
                    -> Velocity {
        // Calculate r1 and r2
        let mut rng = thread_rng();
        let u_dist = Uniform::new(0.0, 1.0);
        let r1 = rng.sample(u_dist);
        let r2 = rng.sample(u_dist);

        // Get global best
        let gbest = self.global_best();

        // Get local best
        let pbest = particle.pbest;

        // Crate new velocity
        Velocity::from_vec(vs)
    }
}
