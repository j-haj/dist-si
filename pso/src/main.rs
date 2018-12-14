pub mod particle;

use particle::Particle;
fn main() {
    let p_bounds : Vec<(f64,f64)> = vec![(-10., 10.), (-10., 10.)];
    let v_bounds : Vec<(f64,f64)> = vec![(-1., 1.), (-1., 1.)];
    let p = Particle::new(p_bounds, v_bounds, true);
    println!("Particle: {:?}", p);
}
