use pso::pso::particle::{Particle, ParticleUpdateMode};
use pso::pso::position::Position;
use pso::pso::velocity::Velocity;

fn square(p: &Position) -> f64 {
    p.coordinates().iter().map(|x| x*x).fold(0., |acc, x| acc + x)
}

fn main() {
    let p_bounds : Vec<(f64,f64)> = vec![(-10., 10.), (-10., 10.)];
    let v_bounds : Vec<(f64,f64)> = vec![(-1., 1.), (-1., 1.)];
    let mode = ParticleUpdateMode::Parallel;
    let mut p = Particle::new(&p_bounds, &v_bounds, square, mode);
    let mut q = Particle::new(&p_bounds, &v_bounds, square, mode);
    let v = Velocity::new(&v_bounds);
    println!("Particle 1's fitness: {:?}", p.get_fitness());
    println!("Particle 2's fitness: {:?}", q.get_fitness());
    println!("Updating particles 1 and 2....");
    p.update_position(&v);
    q.update_position(&v);
    println!("Particle 1's new fitness: {:?}", p.get_fitness());
    println!("Particle 2's new fitness: {:?}", q.get_fitness());
}
