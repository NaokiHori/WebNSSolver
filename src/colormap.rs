#![deny(missing_docs)]

//! Defines functions to map scalar values to RGB color codes.
//!
//! This module converts a [`f64`] value to a set of three [`u8`] values representing `RGB` value.

fn kernel(rgbcoefs: &[[f64; 3]], val: f64) -> [u8; 3] {
    // fit polynomial
    let mut rgb: [f64; 3] = [0., 0., 0.];
    for (n, rgbcoef) in rgbcoefs.iter().enumerate() {
        for m in 0..3 {
            rgb[m] += rgbcoef[m] * val.powi(n as i32);
        }
    }
    // truncate
    for m in 0..3 {
        rgb[m] = if rgb[m] < 0. { 0. } else { rgb[m] };
        rgb[m] = if 1. < rgb[m] { 1. } else { rgb[m] };
    }
    return [
        (255. * rgb[0]) as u8,
        (255. * rgb[1]) as u8,
        (255. * rgb[2]) as u8,
    ];
}

/// Seismic-like color map
/// # Arguments
/// * `val` : Normalised value to be converted.
///
/// # Returns
/// A set of `RGB` values based on the seismic color map.
///
/// # Reference
/// [matplotlib/lib/matplotlib/_cm.py](https://github.com/matplotlib/matplotlib).
pub fn seismic(val: f64) -> [u8; 3] {
    fn linear(val: f64, xs: [f64; 2], ys: [f64; 2]) -> f64 {
        return ys[0] + (ys[1] - ys[0]) / (xs[1] - xs[0]) * (val - xs[0]);
    }
    const XLEVELS: [f64; 5] = [0., 1. / 4., 2. / 4., 3. / 4., 4. / 4.];
    const YLEVELS: [f64; 3] = [1. / 4., 3. / 4., 4. / 4.];
    let r: f64 = if val < XLEVELS[1] {
        YLEVELS[0]
    } else if val < XLEVELS[2] {
        linear(val, [XLEVELS[1], XLEVELS[2]], [YLEVELS[0], YLEVELS[2]])
    } else if val < XLEVELS[3] {
        YLEVELS[2]
    } else if val < XLEVELS[4] {
        linear(val, [XLEVELS[3], XLEVELS[4]], [YLEVELS[2], YLEVELS[1]])
    } else {
        YLEVELS[1]
    };
    let g: f64 = if val < XLEVELS[1] {
        YLEVELS[0]
    } else if val < XLEVELS[2] {
        linear(val, [XLEVELS[1], XLEVELS[2]], [YLEVELS[0], YLEVELS[2]])
    } else if val < XLEVELS[3] {
        linear(val, [XLEVELS[2], XLEVELS[3]], [YLEVELS[2], YLEVELS[0]])
    } else {
        YLEVELS[0]
    };
    let b: f64 = if val < XLEVELS[0] {
        YLEVELS[1]
    } else if val < XLEVELS[1] {
        linear(val, [XLEVELS[0], XLEVELS[1]], [YLEVELS[1], YLEVELS[2]])
    } else if val < XLEVELS[2] {
        YLEVELS[2]
    } else if val < XLEVELS[3] {
        linear(val, [XLEVELS[2], XLEVELS[3]], [YLEVELS[2], YLEVELS[0]])
    } else {
        YLEVELS[0]
    };
    return [(255. * r) as u8, (255. * g) as u8, (255. * b) as u8];
}

/// Viridis color map
/// # Arguments
/// * `val` : Normalised value to be converted.
///
/// # Returns
/// A set of `RGB` values based on the viridis color map.
///
/// # Reference
/// [matplotlib/lib/matplotlib/_cm_listed.py](https://github.com/matplotlib/matplotlib).
pub fn viridis(val: f64) -> [u8; 3] {
    const RGBCOEFS: [[f64; 3]; 5] = [
        [
            2.672303238499781e-01,
            5.015408860973969e-03,
            3.290548382054911e-01,
        ],
        [
            8.867281107764821e-01,
            1.415434679048477e+00,
            6.427369217396137e-01,
        ],
        [
            -6.777660845884058e+00,
            -8.089902514371242e-01,
            2.998258532949060e+00,
        ],
        [
            1.102198635856048e+01,
            7.296293729490473e-01,
            -9.057970794130403e+00,
        ],
        [
            -4.404685706758277e+00,
            -4.355228476501643e-01,
            5.230151793650696e+00,
        ],
    ];
    return kernel(&RGBCOEFS, val);
}

/// Inferno color map
/// # Arguments
/// * `val` : Normalised value to be converted.
///
/// # Returns
/// A set of `RGB` values based on the inferno color map.
///
/// # Reference
/// [matplotlib/lib/matplotlib/_cm_listed.py](https://github.com/matplotlib/matplotlib).
pub fn inferno(val: f64) -> [u8; 3] {
    const RGBCOEFS: [[f64; 3]; 5] = [
        [
            1.449386144388309e-03,
            6.440597326930739e-04,
            1.422820030792134e-02,
        ],
        [
            -5.246766512961221e-02,
            6.032800351043351e-01,
            2.756869387343746e+00,
        ],
        [
            8.344091508648892e+00,
            -2.657917413225927e+00,
            -4.028847551173411e+00,
        ],
        [
            -1.324216881850942e+01,
            6.059021566144984e+00,
            -3.101644861737894e+00,
        ],
        [
            5.935434087778525e+00,
            -3.006240922365244e+00,
            5.006598743510135e+00,
        ],
    ];
    return kernel(&RGBCOEFS, val);
}

/// Cividis color map
/// # Arguments
/// * `val` : Normalised value to be converted.
///
/// # Returns
/// A set of `RGB` values based on the cividis color map.
///
/// # Reference
/// [matplotlib/lib/matplotlib/_cm_listed.py](https://github.com/matplotlib/matplotlib).
pub fn cividis(val: f64) -> [u8; 3] {
    const RGBCOEFS: [[f64; 3]; 4] = [
        [
            -1.892803366663002e-02,
            1.301247649044452e-01,
            3.504586514209761e-01,
        ],
        [
            9.485572691982757e-01,
            7.441875094643144e-01,
            5.745534820352121e-01,
        ],
        [
            2.819352539042831e-01,
            -1.746963935612102e-01,
            -7.409525976291664e-01,
        ],
        [
            -2.395763837375361e-01,
            2.067592501562131e-01,
            8.013828803842242e-02,
        ],
    ];
    return kernel(&RGBCOEFS, val);
}

/// Magma color map
/// # Arguments
/// * `val` : Normalised value to be converted.
///
/// # Returns
/// A set of `RGB` values based on the magma color map.
///
/// # Reference
/// [matplotlib/lib/matplotlib/_cm_listed.py](https://github.com/matplotlib/matplotlib).
pub fn magma(val: f64) -> [u8; 3] {
    const RGBCOEFS: [[f64; 3]; 4] = [
        [
            -4.926960662072657e-02,
            4.728491591684356e-02,
            -7.737655108244056e-03,
        ],
        [
            1.153630996787603e+00,
            -6.370045947380021e-02,
            3.223330701641625e+00,
        ],
        [
            1.706765451132571e+00,
            5.491051576171306e-01,
            -6.634189647774892e+00,
        ],
        [
            -1.868698117409254e+00,
            4.985349690285930e-01,
            4.142809899823387e+00,
        ],
    ];
    return kernel(&RGBCOEFS, val);
}

/// Turbo color map
/// # Arguments
/// * `val` : Normalised value to be converted.
///
/// # Returns
/// A set of `RGB` values based on the turbo color map.
///
/// # Reference
/// [matplotlib/lib/matplotlib/_cm_listed.py](https://github.com/matplotlib/matplotlib).
pub fn turbo(val: f64) -> [u8; 3] {
    const RGBCOEFS: [[f64; 3]; 5] = [
        [
            1.891320428997353e-01,
            7.171830990657042e-02,
            2.326283869688164e-01,
        ],
        [
            -3.227425752012548e-01,
            2.737858189663653e+00,
            7.920493523542985e+00,
        ],
        [
            1.319896573632580e+00,
            2.807408919718347e+00,
            -2.877955822159493e+01,
        ],
        [
            4.454585573767376e+00,
            -1.341530665932009e+01,
            3.352480296358407e+01,
        ],
        [
            -5.163443190685617e+00,
            7.815652904752115e+00,
            -1.288855294069439e+01,
        ],
    ];
    return kernel(&RGBCOEFS, val);
}
