use std::io::Write;
use web_ns_solver::{Array, Field};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    const RA: f64 = 8.5e7;
    const PR: f64 = 4.4e0;
    const LENGTHS: [f64; web_ns_solver::NDIMS] = [2f64, 1f64];
    const GRID_POINTS: [usize; web_ns_solver::NDIMS] = [128usize, 64usize];
    const TIME_MAX: f64 = 1e+1;
    const LOG_RATE: f64 = 1e+0;
    let mut field = Field::new(
        RA,
        PR,
        LENGTHS[0],
        LENGTHS[1],
        GRID_POINTS[0],
        GRID_POINTS[1],
    );
    let mut niter = 0u32;
    let mut time = 0f64;
    let mut log_next: f64 = LOG_RATE;
    // main loop
    loop {
        let dt: f64 = field.evolve();
        niter += 1;
        time += dt;
        if log_next < time {
            let div: f64 = field.check_div();
            println!("step: {niter:8} time: {time:.1e}, dt: {dt:.1e}, div: {div:.1e}");
            log_next += LOG_RATE;
        }
        if TIME_MAX < time {
            break;
        }
    }
    {
        let te: &Array = &field.get_te();
        let file_name = "sample.npy";
        let mut f = std::fs::File::create(file_name).unwrap();
        let shape: [usize; 2] = [GRID_POINTS[1] + 2usize, GRID_POINTS[0] + 2usize];
        let header = rust_npy_io::Header {
            descr: "'<f8'".to_string(),
            fortran_order: false,
            shape: shape.to_vec(),
        };
        match rust_npy_io::write_header(&mut f, &header) {
            Ok(_) => {}
            Err(e) => {
                println!("Failed to write header: {}", e);
                std::process::exit(1);
            }
        };
        for datum in te.get_data() {
            let bytes = datum.to_ne_bytes();
            f.write_all(&bytes)?;
        }
    }
    Ok(())
}
