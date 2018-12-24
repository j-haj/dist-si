use ndarray;
use ndarray::Array;

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
    let mut p = Particle::new(&p_bounds, &v_bounds, square, mode,
                              1., 1., 1.8, true);
    let v = Velocity::new(v_bounds.len());
    println!("Particle 1's fitness: {:?}", p.fitness());
    println!("Updating particles 1 and 2....");
    let g_best = Position::from_vec(&vec![0., 0.]);
    for i in 0..10000 {
        p.update_velocity(&g_best);
        p.update_position();
        if i % 100 == 0 {
            println!("Round {:?}", i);
            println!("\tfitness = {:?}", p.fitness());
        }
    }
}
