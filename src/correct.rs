use super::dct;
use super::tdm;
use super::Array;
use super::Field;

pub fn correct(dt: f64, field: &mut Field) {
    let nx: usize = field.npoints[0];
    let ny: usize = field.npoints[1];
    let dxi: f64 = 1. / field.deltas[0];
    let dyi: f64 = 1. / field.deltas[1];
    let ux: &mut Array = &mut field.ux;
    let uy: &mut Array = &mut field.uy;
    let sc0: &mut Vec<f64> = field.sc0.get_data_mut();
    let sc1: &mut Vec<f64> = field.sc1.get_data_mut();
    let dct_plan: &mut dct::Plan = &mut field.dct_plan;
    let tdm_plan: &mut tdm::Plan = &mut field.tdm_plan;
    // compute local divergence
    // DCT normalisation (1 / 2N) is done here as well
    let coef: f64 = 1. / dt / (2 * nx) as f64;
    for j in 0..ny {
        for i in 0..nx {
            let duxdx: f64 = (-ux[j + 1][i] + ux[j + 1][i + 1]) * dxi;
            let duydy: f64 = (-uy[j][i + 1] + uy[j + 1][i + 1]) * dyi;
            let div: f64 = coef * (duxdx + duydy);
            // reorder input buffer for efficient DCT
            let ii: usize = dct_plan.map_index(i);
            sc0[j * nx + ii] = div;
        }
    }
    // solve Poisson equation
    // dct in x
    for j in 0..ny {
        let idx: usize = j * nx;
        dct_plan.exec_f(&sc0[idx..], &mut sc1[idx..]);
    }
    // compute left-hand side
    let dyi2: f64 = dyi * dyi;
    let waves: &Vec<f64> = &field.waves;
    for i in 0..nx {
        // prepare centre-diagonal components of the linear system
        let c: &mut Vec<f64> = &mut tdm_plan.c;
        let wave: f64 = waves[i];
        c[0] = wave - dyi2;
        for j in 1..ny - 1 {
            c[j] = wave - 2. * dyi2;
        }
        c[ny - 1] = wave - dyi2;
        // invoke tri-diagonal matrix solver
        tdm_plan.solve(&mut sc1[i..]);
    }
    // idct in x
    for j in 0..ny {
        let idx: usize = j * nx;
        dct_plan.exec_b(&sc1[idx..], &mut sc0[idx..]);
    }
    // velocity correction
    let coef: f64 = dt * dxi;
    for j in 0..ny {
        for i in 0..nx - 1 {
            // reorder output buffer for efficient DCT
            let im: usize = dct_plan.map_index(i);
            let ip: usize = dct_plan.map_index(i + 1);
            let pm: f64 = sc0[j * nx + im];
            let pp: f64 = sc0[j * nx + ip];
            ux[j + 1][i + 1] -= coef * (-pm + pp);
        }
    }
    let coef: f64 = dt * dyi;
    for j in 0..ny - 1 {
        for i in 0..nx {
            // reorder output buffer for efficient DCT
            let ii: usize = dct_plan.map_index(i);
            let pm: f64 = sc0[j * nx + ii];
            let pp: f64 = sc0[(j + 1) * nx + ii];
            uy[j + 1][i + 1] -= coef * (-pm + pp);
        }
    }
}
