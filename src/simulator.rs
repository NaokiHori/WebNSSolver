#![deny(missing_docs)]

//! Core routines to handle the Navier-Stokes equation.
//!
//! # Overview
//! This module integrates the Navier-Stokes equation in time by means of the fractional-step method.
//!
//! # Caveat
//! The main objective of this crate is to yield something plausible; important aspects are the smoothness of the flow fields and the efficiency to achieve a reasonable FPS.
//! To this end I introduce numerical viscosities and assume equidistant grid sizes, which may deteriorate the quality of the solution.
//! As a result, *using this library for rigorous objectives is not recommended*.

use crate::array::Array;

/// Number of grid points, x and y.
pub const NPOINTS: [usize; 2] = [128, 128];

/// Domain lengths, x and y.
const LENGTHS: [f64; 2] = [1., 1.];

/// Grid sizes, x and y (assumed to be equidistant).
const DELTAS: [f64; 2] = [
    LENGTHS[0] / NPOINTS[0] as f64,
    LENGTHS[1] / NPOINTS[1] as f64,
];

/// Number of tracer particle history.
pub const NHISTORY: usize = 128;
const NTRACERS: [usize; 2] = [NPOINTS[0] / 16, NPOINTS[1] / 16];

/// Numerical viscosities, momentum / temperature.
///
/// Positive numbers are given to stabilise the integration by dumping out the non-linear oscillation.
/// Setting them to zero recovers the energy-conserving second-order central finite-difference schemes.
const NUMVIS: [f64; 2] = [0.10, 0.15];

/// Stores state of a flow field.
pub struct Field {
    /// X velocity
    ux: Array,
    /// Y velocity
    uy: Array,
    /// Temperature
    te: Array,
    /// X velocity increment (ux += dux * dt)
    dux: Array,
    /// Y velocity increment (uy += duy * dt)
    duy: Array,
    /// Temperature increment (te += dte * dt)
    dte: Array,
    /// Scalar potential for the correction step, buffer No. 0
    sc0: Array,
    /// Scalar potential for the correction step, buffer No. 1
    sc1: Array,
    /// Plan for the discrete cosine transforms
    dct_plan: crate::dct::Plan,
    /// Plan for the linear solver
    tdm_plan: crate::tdm::Plan,
    /// Modified wave numbers describing the discrete Poisson equation in x
    waves: Vec<f64>,
    /// Normalised momentum diffusivity
    visu: f64,
    /// Normalised temperature diffusivity
    vist: f64,
    /// Tracer positions
    pub tracers: Vec<[[f64; 2]; NHISTORY]>,
}

impl Field {
    /// Constructor of [`Field`] struct.
    ///
    /// # Arguments
    /// * `ra` : Rayleigh number.
    /// * `pr` : Prandtl number.
    ///
    /// # Returns
    /// A new initialised array.
    pub fn new(ra: f64, pr: f64) -> Field {
        // abbreviations
        let nx: usize = NPOINTS[0];
        let ny: usize = NPOINTS[1];
        let dx: f64 = DELTAS[0];
        let dy: f64 = DELTAS[1];
        // make scalar fields
        let mut ux = Array::new(nx + 1, ny + 2);
        let mut uy = Array::new(nx + 2, ny + 1);
        let mut te = Array::new(nx + 2, ny + 2);
        let dux = Array::new(nx + 1, ny + 2);
        let duy = Array::new(nx + 2, ny + 1);
        let dte = Array::new(nx + 2, ny + 2);
        // give BCs
        for j in 1..ny + 1 {
            ux[j][0] = 0.;
            ux[j][nx] = 0.;
        }
        for i in 0..nx + 1 {
            ux[0][i] = 0.;
            ux[ny + 1][i] = 0.;
        }
        for j in 0..ny + 1 {
            uy[j][0] = 0.;
            uy[j][nx + 1] = 0.;
        }
        for i in 1..nx + 1 {
            uy[0][i] = 0.;
            uy[ny][i] = 0.;
        }
        for j in 1..ny + 1 {
            te[j][0] = 1.;
            te[j][nx + 1] = 0.;
        }
        for i in 1..nx + 1 {
            te[0][i] = 0.;
            te[ny + 1][i] = 1.;
        }
        // prepare buffers of the scalar potential
        let sc0 = Array::new(nx + 0, ny + 0);
        let sc1 = Array::new(nx + 0, ny + 0);
        let dct_plan = crate::dct::Plan::new(nx);
        let mut tdm_plan = crate::tdm::Plan::new(ny, nx);
        let l: &mut Vec<f64> = &mut tdm_plan.l;
        let u: &mut Vec<f64> = &mut tdm_plan.u;
        for j in 0..ny {
            l[j] = 1. / dy / dy;
            u[j] = 1. / dy / dy;
        }
        let mut waves: Vec<f64> = vec![0.; nx];
        for i in 0..nx {
            // NOTE: for cosine transform, the signal size becomes 2N
            waves[i] =
                -(2. / dx * (std::f64::consts::PI * i as f64 / 2. / nx as f64).sin()).powi(2);
        }
        // diffusivities
        let visu: f64 = (pr / ra).sqrt();
        let vist: f64 = (1. / pr / ra).sqrt();
        // tracers
        let mut tracers = Vec::<[[f64; 2]; NHISTORY]>::new();
        for j in 0..NTRACERS[1] {
            let dy: f64 = LENGTHS[1] / NTRACERS[1] as f64;
            let y: f64 = 0.5 * (2 * j + 1) as f64 * dy;
            for i in 0..NTRACERS[0] {
                let dx: f64 = LENGTHS[0] / NTRACERS[0] as f64;
                let x: f64 = 0.5 * (2 * i + 1) as f64 * dx;
                let tracer: [[f64; 2]; NHISTORY] = [[x, y]; NHISTORY];
                tracers.push(tracer);
            }
        }
        Field {
            ux,
            uy,
            te,
            dux,
            duy,
            dte,
            sc0,
            sc1,
            dct_plan,
            tdm_plan,
            waves,
            visu,
            vist,
            tracers,
        }
    }

    /// Getter, [`Field::ux`].
    ///
    /// # Returns
    /// A reference to [`Field::ux`] member.
    #[allow(dead_code)]
    pub fn get_ux(&self) -> &Array {
        return &self.ux;
    }

    /// Getter, [`Field::uy`].
    ///
    /// # Returns
    /// A reference to [`Field::uy`] member.
    #[allow(dead_code)]
    pub fn get_uy(&self) -> &Array {
        return &self.uy;
    }

    /// Getter, [`Field::te`].
    ///
    /// # Returns
    /// A reference to [`Field::te`] member.
    #[allow(dead_code)]
    pub fn get_te(&self) -> &Array {
        return &self.te;
    }
}

/// Computes the time step size for the next iteration.
///
/// # Arguments
/// * `dt_max` - the biggest time step size.
///
/// # Returns
/// Time step size for the next iteration.
fn find_dt(dt_max: f64, field: &Field) -> f64 {
    // safety factors for the advective / diffusive components
    const FACTORS: [f64; 2] = [0.15, 0.125];
    // minimum grid size, just in case
    let delta: f64 = DELTAS[0].min(DELTAS[1]);
    // diffusivities
    let visui: f64 = 1. / field.visu;
    let visti: f64 = 1. / field.vist;
    // get extrema: max(|ux|) and max(|uy|)
    let ux: f64 = field.ux.get_ext();
    let uy: f64 = field.uy.get_ext();
    let ext: f64 = if ux > uy { ux } else { uy };
    // advective restriction
    let dt_adv: f64 = FACTORS[0] * delta / ext;
    // diffusive restriction
    let dt_dif: f64 = FACTORS[1] * visui.min(visti) * delta.powi(2);
    // find resulting dt
    let mut dt: f64 = dt_max;
    dt = if dt < dt_adv { dt } else { dt_adv };
    dt = if dt < dt_dif { dt } else { dt_dif };
    return dt;
}

/// The prediction step of the fractional-step method.
fn predict(dt: f64, field: &mut Field) -> () {
    let nx: usize = NPOINTS[0];
    let ny: usize = NPOINTS[1];
    let dx: f64 = DELTAS[0];
    let dy: f64 = DELTAS[1];
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
            let l: f64 = (0. + dxiq) * (ux[j + 0][i - 1] + ux[j + 0][i + 0]);
            let u: f64 = (0. - dxiq) * (ux[j + 0][i + 0] + ux[j + 0][i + 1]);
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
            let l: f64 = (0. + dyiq) * (uy[j - 1][i + 0] + uy[j - 1][i + 1]);
            let u: f64 = (0. - dyiq) * (uy[j + 0][i + 0] + uy[j + 0][i + 1]);
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
            let l: f64 = (0. + dxiq) * (ux[j + 0][i - 1] + ux[j + 1][i - 1]);
            let u: f64 = (0. - dxiq) * (ux[j + 0][i + 0] + ux[j + 1][i + 0]);
            duy[j][i] += l * uy[j][i - 1] - (l + u) * uy[j][i] + u * uy[j][i + 1];
            let mag: f64 = NUMVIS[0] * dx2 * (l - u).abs();
            let l: f64 = (visu + mag) * dxi2 * if 1 == i { 2. } else { 1. };
            let u: f64 = (visu + mag) * dxi2 * if nx == i { 2. } else { 1. };
            duy[j][i] += l * uy[j][i - 1] - (l + u) * uy[j][i] + u * uy[j][i + 1];
            // y advection / diffusion
            //   +--^--+--^--+
            // j-1 i  j i  j+1 i
            //   +--^--+--^--+
            let l: f64 = (0. + dyiq) * (uy[j - 1][i + 0] + uy[j + 0][i + 0]);
            let u: f64 = (0. - dyiq) * (uy[j + 0][i + 0] + uy[j + 1][i + 0]);
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
            let l: f64 = (0. + dxih) * ux[j + 0][i - 1];
            let u: f64 = (0. - dxih) * ux[j + 0][i + 0];
            dte[j][i] += l * te[j][i - 1] - (l + u) * te[j][i] + u * te[j][i + 1];
            let mag: f64 = NUMVIS[1] * dx2 * (l - u).abs();
            let l: f64 = (vist + mag) * dxi2 * if 1 == i { 2. } else { 1. };
            let u: f64 = (vist + mag) * dxi2 * if nx == i { 2. } else { 1. };
            dte[j][i] += l * te[j][i - 1] - (l + u) * te[j][i] + u * te[j][i + 1];
            // y advection / diffusion
            //   +--^--+
            // j-1 i  j i
            //   +--^--+
            let l: f64 = (0. + dyih) * uy[j - 1][i + 0];
            let u: f64 = (0. - dyih) * uy[j + 0][i + 0];
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

/// The correction step of the fractional-step method.
fn correct(dt: f64, field: &mut Field) -> () {
    let nx: usize = NPOINTS[0];
    let ny: usize = NPOINTS[1];
    let dxi: f64 = 1. / DELTAS[0];
    let dyi: f64 = 1. / DELTAS[1];
    let ux: &mut Array = &mut field.ux;
    let uy: &mut Array = &mut field.uy;
    let sc0: &mut Vec<f64> = &mut field.sc0.get_data_mut();
    let sc1: &mut Vec<f64> = &mut field.sc1.get_data_mut();
    let dct_plan: &mut crate::dct::Plan = &mut field.dct_plan;
    let tdm_plan: &mut crate::tdm::Plan = &mut field.tdm_plan;
    // compute local divergence
    // DCT normalisation (1 / 2N) is done here as well
    let coef: f64 = 1. / dt / (2 * nx) as f64;
    for j in 0..ny {
        for i in 0..nx {
            let duxdx: f64 = (-ux[j + 1][i + 0] + ux[j + 1][i + 1]) * dxi;
            let duydy: f64 = (-uy[j + 0][i + 1] + uy[j + 1][i + 1]) * dyi;
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
            let im: usize = dct_plan.map_index(i + 0);
            let ip: usize = dct_plan.map_index(i + 1);
            let pm: f64 = sc0[(j + 0) * nx + im];
            let pp: f64 = sc0[(j + 0) * nx + ip];
            ux[j + 1][i + 1] -= coef * (-pm + pp);
        }
    }
    let coef: f64 = dt * dyi;
    for j in 0..ny - 1 {
        for i in 0..nx {
            // reorder output buffer for efficient DCT
            let ii: usize = dct_plan.map_index(i);
            let pm: f64 = sc0[(j + 0) * nx + ii];
            let pp: f64 = sc0[(j + 1) * nx + ii];
            uy[j + 1][i + 1] -= coef * (-pm + pp);
        }
    }
}

/// Updating tracer particles.
fn update_tracers(dt: f64, field: &mut Field) -> () {
    let ux: &Array = &field.ux;
    let uy: &Array = &field.uy;
    let tracers: &mut Vec<[[f64; 2]; NHISTORY]> = &mut field.tracers;
    for tracer in tracers.iter_mut() {
        // take latest information
        let x: f64 = tracer[NHISTORY - 1][0];
        let y: f64 = tracer[NHISTORY - 1][1];
        let i: usize = (x / DELTAS[0]) as usize;
        let j: usize = (y / DELTAS[1]) as usize;
        let ux: f64 = 0.5 * ux[j + 1][i + 0] + 0.5 * ux[j + 1][i + 1];
        let uy: f64 = 0.5 * uy[j + 0][i + 1] + 0.5 * uy[j + 1][i + 1];
        let mut x: f64 = x + ux * dt;
        let mut y: f64 = y + uy * dt;
        if x < 0.5 * DELTAS[0] {
            x = DELTAS[0];
        }
        if LENGTHS[0] - DELTAS[0] < x {
            x = LENGTHS[0] - DELTAS[0];
        }
        if y < 0.5 * DELTAS[1] {
            y = DELTAS[1];
        }
        if LENGTHS[1] - DELTAS[1] < y {
            y = LENGTHS[1] - DELTAS[1];
        }
        // discard the head, add the latest position to the tail
        tracer.rotate_left(1);
        tracer[NHISTORY - 1] = [x, y];
    }
}

/// Entry point.
pub fn evolve(dt_max: f64, field: &mut Field) -> f64 {
    let dt: f64 = find_dt(dt_max, field);
    predict(dt, field);
    correct(dt, field);
    update_tracers(dt, field);
    return dt;
}

/// An auxiliary function to compute / monitor the maximum divergence of the flow field.
#[allow(dead_code)]
pub fn check_div(field: &Field) -> f64 {
    let nx: usize = NPOINTS[0];
    let ny: usize = NPOINTS[1];
    let dxi: f64 = 1. / DELTAS[0];
    let dyi: f64 = 1. / DELTAS[1];
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
    return maxdiv;
}
