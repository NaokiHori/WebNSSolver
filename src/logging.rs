use super::Array;
use super::Field;

pub fn check_div(field: &Field) -> f64 {
    let nx: usize = field.npoints[0];
    let ny: usize = field.npoints[1];
    let dxi: f64 = 1. / field.deltas[0];
    let dyi: f64 = 1. / field.deltas[1];
    let ux: &Array = &field.ux;
    let uy: &Array = &field.uy;
    let mut maxdiv: f64 = 0.;
    for j in 1..ny + 1 {
        for i in 1..nx + 1 {
            let duxdx: f64 = (ux[j][i] - ux[j][i - 1]) * dxi;
            let duydy: f64 = (uy[j][i] - uy[j - 1][i]) * dyi;
            let mut div: f64 = duxdx + duydy;
            div = if div < -1. * div { -1. * div } else { div };
            maxdiv = if div < maxdiv { maxdiv } else { div };
        }
    }
    maxdiv
}
