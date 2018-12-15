use pso::particle::{Particle,Position};

fn square(p: &Position<f64>) -> f64 {
    p.get_coordinates().iter().map(|x| x*x).fold(0., |acc, x| acc + x)
}

fn main() {
    let p_bounds : Vec<(f64,f64)> = vec![(-10., 10.), (-10., 10.)];
    let v_bounds : Vec<(f64,f64)> = vec![(-1., 1.), (-1., 1.)];
    let p = Particle::new(&p_bounds, &v_bounds, square);
    let q = Particle::new(&p_bounds, &v_bounds, square);
    println!("Particle 1's fitness: {:?}", p.get_fitness());
    println!("Particle 2's fitness: {:?}", q.get_fitness());
}
