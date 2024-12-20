use super::Field;

const DT_MAX: f64 = 1e-1f64;

// safety factors for the advective / diffusive components
const FACTORS: [f64; 2] = [0.25, 0.50];

pub fn find_dt(field: &Field) -> f64 {
    // minimum grid size, just in case
    let delta: f64 = field.deltas[0].min(field.deltas[1]);
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
    let mut dt: f64 = DT_MAX;
    dt = if dt < dt_adv { dt } else { dt_adv };
    dt = if dt < dt_dif { dt } else { dt_dif };
    dt
}
