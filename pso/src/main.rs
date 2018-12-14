use pso::particle::Particle;
fn main() {
    let p_bounds : Vec<(f64,f64)> = vec![(-10., 10.), (-10., 10.)];
    let v_bounds : Vec<(f64,f64)> = vec![(-1., 1.), (-1., 1.)];
    let p = Particle::new(&p_bounds, &v_bounds);
    let q = Particle::new(&p_bounds, &v_bounds);
    println!("Particle: {:?}", p);
    println!("Particle: {:?}", q);
}
