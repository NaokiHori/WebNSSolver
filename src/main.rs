#![deny(missing_docs)]

//! Entry point of a binary crate.

mod array;
mod dct;
mod simulator;
mod tdm;

fn main() {
    const RA: f64 = 8.5e7;
    const PR: f64 = 4.4e0;
    const TIME_MAX: f64 = 1e+1;
    const LOG_RATE: f64 = 1e+0;
    // initialise simulator and obtain the initial flow field
    let mut field: simulator::Field = simulator::Field::new(RA, PR);
    // schdulers
    let mut niter: u32 = 0;
    let mut time: f64 = 0.;
    let mut log_next: f64 = LOG_RATE;
    // main loop
    loop {
        // integrate (the first argument is the maximum time step size)
        let dt: f64 = simulator::evolve(LOG_RATE, &mut field);
        niter += 1;
        time += dt;
        if log_next < time {
            // compute statistics
            let div: f64 = simulator::check_div(&field);
            let ux: f64 = field.get_ux().get_ext();
            let uy: f64 = field.get_uy().get_ext();
            println!("step: {niter:8} time: {time:.1e}, dt: {dt:.1e}, div: {div:.1e}, max: ({ux:.1e}, {uy:.1e})");
            log_next += LOG_RATE;
        }
        if TIME_MAX < time {
            // completed
            break;
        }
    }
}
