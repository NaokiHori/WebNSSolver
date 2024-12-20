use super::Array;
use super::{Field, NUMVIS};

pub fn predict(dt: f64, field: &mut Field) {
    let nx: usize = field.npoints[0];
    let ny: usize = field.npoints[1];
    let dx: f64 = field.deltas[0];
    let dy: f64 = field.deltas[1];
    let dxi: f64 = 1. / dx;
    let dyi: f64 = 1. / dy;
    let dx2: f64 = dx.powi(2);
    let dy2: f64 = dy.powi(2);
    let dxi2: f64 = dxi.powi(2);
    let dyi2: f64 = dyi.powi(2);
    let dxih: f64 = 0.5 * dxi;
    let dyih: f64 = 0.5 * dyi;
    let dxiq: f64 = 0.25 * dxi;
    let dyiq: f64 = 0.25 * dyi;
    let visu: f64 = field.visu;
    let vist: f64 = field.vist;
    let ux: &mut Array = &mut field.ux;
    let uy: &mut Array = &mut field.uy;
    let te: &mut Array = &mut field.te;
    let dux: &mut Array = &mut field.dux;
    let duy: &mut Array = &mut field.duy;
    let dte: &mut Array = &mut field.dte;
    // compute dux
    for j in 1..ny + 1 {
        for i in 1..nx {
            dux[j][i] = 0.;
            // x advection / diffusion
            //   +--^--+--^--+
            // i-1 j  i j  i+1 j
            //   +--^--+--^--+
            let l: f64 = (0. + dxiq) * (ux[j][i - 1] + ux[j][i]);
            let u: f64 = (0. - dxiq) * (ux[j][i] + ux[j][i + 1]);
            dux[j][i] += l * ux[j][i - 1] - (l + u) * ux[j][i] + u * ux[j][i + 1];
            let mag: f64 = NUMVIS[0] * dx * ux[j][i].abs();
            let l: f64 = (visu + mag) * dxi2;
            let u: f64 = (visu + mag) * dxi2;
            dux[j][i] += l * ux[j][i - 1] - (l + u) * ux[j][i] + u * ux[j][i + 1];
            // y advection / diffusion
            //     i j  i+1 j
            //   +--^--+--^--+
            //   |    i j    |
            //   +--^--+--^--+
            //   i j-1  i+1 j-1
            // NOTE: for advection, l/u corrections are ignored,
            //         assuming impermeable walls
            let l: f64 = (0. + dyiq) * (uy[j - 1][i] + uy[j - 1][i + 1]);
            let u: f64 = (0. - dyiq) * (uy[j][i] + uy[j][i + 1]);
            dux[j][i] += l * ux[j - 1][i] - (l + u) * ux[j][i] + u * ux[j + 1][i];
            let mag: f64 = NUMVIS[0] * dy2 * (l - u).abs();
            let l: f64 = (visu + mag) * dyi2 * if 1 == j { 2. } else { 1. };
            let u: f64 = (visu + mag) * dyi2 * if ny == j { 2. } else { 1. };
            dux[j][i] += l * ux[j - 1][i] - (l + u) * ux[j][i] + u * ux[j + 1][i];
        }
    }
    // compute duy
    for j in 1..ny {
        for i in 1..nx + 1 {
            duy[j][i] = 0.;
            // x advection / diffusion
            //     j i  j+1 i
            //   +--^--+--^--+
            //   |    j i    |
            //   +--^--+--^--+
            //   j i-1  j+1 i-1
            // NOTE: for advection, l/u corrections are ignored,
            //         assuming impermeable walls
            let l: f64 = (0. + dxiq) * (ux[j][i - 1] + ux[j + 1][i - 1]);
            let u: f64 = (0. - dxiq) * (ux[j][i] + ux[j + 1][i]);
            duy[j][i] += l * uy[j][i - 1] - (l + u) * uy[j][i] + u * uy[j][i + 1];
            let mag: f64 = NUMVIS[0] * dx2 * (l - u).abs();
            let l: f64 = (visu + mag) * dxi2 * if 1 == i { 2. } else { 1. };
            let u: f64 = (visu + mag) * dxi2 * if nx == i { 2. } else { 1. };
            duy[j][i] += l * uy[j][i - 1] - (l + u) * uy[j][i] + u * uy[j][i + 1];
            // y advection / diffusion
            //   +--^--+--^--+
            // j-1 i  j i  j+1 i
            //   +--^--+--^--+
            let l: f64 = (0. + dyiq) * (uy[j - 1][i] + uy[j][i]);
            let u: f64 = (0. - dyiq) * (uy[j][i] + uy[j + 1][i]);
            duy[j][i] += l * uy[j - 1][i] - (l + u) * uy[j][i] + u * uy[j + 1][i];
            let mag: f64 = NUMVIS[0] * dy * uy[j][i].abs();
            let l: f64 = (visu + mag) * dyi2;
            let u: f64 = (visu + mag) * dyi2;
            duy[j][i] += l * uy[j - 1][i] - (l + u) * uy[j][i] + u * uy[j + 1][i];
            // buoyancy
            duy[j][i] -= 0.5 * (te[j][i] + te[j + 1][i]);
        }
    }
    // compute dte
    for j in 1..ny + 1 {
        for i in 1..nx + 1 {
            dte[j][i] = 0.;
            // x advection / diffusion
            //   +--^--+
            // i-1 j  i j
            //   +--^--+
            let l: f64 = (0. + dxih) * ux[j][i - 1];
            let u: f64 = (0. - dxih) * ux[j][i];
            dte[j][i] += l * te[j][i - 1] - (l + u) * te[j][i] + u * te[j][i + 1];
            let mag: f64 = NUMVIS[1] * dx2 * (l - u).abs();
            let l: f64 = (vist + mag) * dxi2 * if 1 == i { 2. } else { 1. };
            let u: f64 = (vist + mag) * dxi2 * if nx == i { 2. } else { 1. };
            dte[j][i] += l * te[j][i - 1] - (l + u) * te[j][i] + u * te[j][i + 1];
            // y advection / diffusion
            //   +--^--+
            // j-1 i  j i
            //   +--^--+
            let l: f64 = (0. + dyih) * uy[j - 1][i];
            let u: f64 = (0. - dyih) * uy[j][i];
            dte[j][i] += l * te[j - 1][i] - (l + u) * te[j][i] + u * te[j + 1][i];
            let mag: f64 = NUMVIS[1] * dy2 * (l - u).abs();
            let l: f64 = (vist + mag) * dyi2 * if 1 == j { 2. } else { 1. };
            let u: f64 = (vist + mag) * dyi2 * if ny == j { 2. } else { 1. };
            dte[j][i] += l * te[j - 1][i] - (l + u) * te[j][i] + u * te[j + 1][i];
        }
    }
    // update ux
    for j in 1..ny + 1 {
        for i in 1..nx {
            ux[j][i] += dux[j][i] * dt;
        }
    }
    // update uy
    for j in 1..ny {
        for i in 1..nx + 1 {
            uy[j][i] += duy[j][i] * dt;
        }
    }
    // update te
    for j in 1..ny + 1 {
        for i in 1..nx + 1 {
            te[j][i] += dte[j][i] * dt;
        }
    }
}
