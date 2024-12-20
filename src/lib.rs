mod array;
mod correct;
mod dct;
mod find_dt;
mod logging;
mod predict;
mod tdm;

use wasm_bindgen::prelude::*;

pub use array::Array;

pub const NDIMS: usize = 2usize;

pub const N_HALO: usize = 1usize;

const NUMVIS: [f64; 2] = [0.10, 0.15];

#[wasm_bindgen]
pub enum ScalarType {
    Temperature = 0,
    XVelocity = 1,
    YVelocity = 2,
}

#[cfg_attr(not(feature = "binary_crate"), wasm_bindgen)]
pub struct Field {
    npoints: [usize; NDIMS],
    deltas: [f64; NDIMS],
    ux: Array,
    uy: Array,
    te: Array,
    dux: Array,
    duy: Array,
    dte: Array,
    sc0: Array,
    sc1: Array,
    dct_plan: dct::Plan,
    tdm_plan: tdm::Plan,
    waves: Vec<f64>,
    visu: f64,
    vist: f64,
    buf: Vec<u8>,
    scalar_type: ScalarType,
}

#[cfg_attr(not(feature = "binary_crate"), wasm_bindgen)]
impl Field {
    #[cfg_attr(not(feature = "binary_crate"), wasm_bindgen(constructor))]
    pub fn new(ra: f64, pr: f64, lx: f64, ly: f64, nx: usize, ny: usize) -> Field {
        let dx: f64 = lx / nx as f64;
        let dy: f64 = ly / ny as f64;
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
        let sc0 = Array::new(nx, ny);
        let sc1 = Array::new(nx, ny);
        let dct_plan = dct::Plan::new(nx);
        let mut tdm_plan = tdm::Plan::new(ny, nx);
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
        Field {
            npoints: [nx, ny],
            deltas: [dx, dy],
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
            buf: vec![0u8; nx * ny],
            scalar_type: ScalarType::Temperature,
        }
    }

    pub fn evolve(&mut self) -> f64 {
        let dt: f64 = find_dt::find_dt(self);
        predict::predict(dt, self);
        correct::correct(dt, self);
        dt
    }

    pub fn check_div(&self) -> f64 {
        logging::check_div(self)
    }

    #[cfg(feature = "binary_crate")]
    pub fn get_te(&self) -> &Array {
        &self.te
    }

    #[cfg_attr(not(feature = "binary_crate"), wasm_bindgen)]
    pub fn set_scalar_type(&mut self, scalar_type: ScalarType) {
        self.scalar_type = scalar_type;
    }

    pub fn get_buf_ptr_u8(&mut self) -> *const u8 {
        let npoints: &[usize; NDIMS] = &self.npoints;
        let nx = npoints[0];
        let ny = npoints[1];
        let buf = &mut self.buf;
        match self.scalar_type {
            ScalarType::Temperature => {
                let te: &Array = &self.te;
                let extrema: [f64; 2] = [0f64, 1f64];
                for j in 0..ny {
                    for i in 0..nx {
                        let value: f64 = te[j + 1][i + 1];
                        let value: f64 = (value - extrema[0]) / (extrema[1] - extrema[0]);
                        let value: f64 = 1f64 - value.clamp(0f64, 1f64);
                        let value: f64 = 255f64 * value;
                        buf[j * nx + i] = value as u8;
                    }
                }
            }
            ScalarType::XVelocity => {
                let ux: &Array = &self.ux;
                let extrema: [f64; 2] = [-0.25f64, 0.25f64];
                for j in 0..ny {
                    for i in 0..nx {
                        let value: f64 = 0.5f64 * ux[j + 1][i] + 0.5f64 * ux[j + 1][i + 1];
                        let value: f64 = (value - extrema[0]) / (extrema[1] - extrema[0]);
                        let value: f64 = value.clamp(0f64, 1f64);
                        let value: f64 = 255f64 * value;
                        buf[j * nx + i] = value as u8;
                    }
                }
            }
            ScalarType::YVelocity => {
                let uy: &Array = &self.uy;
                let extrema: [f64; 2] = [-0.25f64, 0.25f64];
                for j in 0..ny {
                    for i in 0..nx {
                        let value: f64 = 0.5f64 * uy[j][i + 1] + 0.5f64 * uy[j + 1][i + 1];
                        let value: f64 = (value - extrema[0]) / (extrema[1] - extrema[0]);
                        let value: f64 = 1f64 - value.clamp(0f64, 1f64);
                        let value: f64 = 255f64 * value;
                        buf[j * nx + i] = value as u8;
                    }
                }
            }
        };
        buf.as_ptr()
    }
}

#[wasm_bindgen(start)]
pub fn init() {}
