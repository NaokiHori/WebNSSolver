#![deny(missing_docs)]

//! Discrete Cosine Transform of the type 2 and 3.
//!
//! # Overview
//! I need to solve a discrete Poisson equations to satisfy the incompressibility, where a discrete cosine transform is adopted in this project.
//!
//! # Caveats
//! The signal size `nitems` should be one of
//! * 4 x 2^N,
//! * 12 x 2^N,
//! where `N` is a non-negative integer.
//!
//! # Reference
//! * Cooley and Tukey, "An algorithm for the machine calculation of complex Fourier series", *Math. Comput.*, 1965
//! * Makhoul, "A Fast Cosine Transform in One and Two Dimensions", *IEEE T. Acoust. Speech*, 1980
//! * [Wikipedia - Cooley–Tukey FFT algorithm](https://en.wikipedia.org/wiki/Cooley–Tukey_FFT_algorithm)
pub struct Plan {
    // size of the input / output signals
    nitems: usize,
    // internal buffers in which
    //   real and imag numbers alternate,
    //   i.e. r[0] i[0] r[1] i[1] ...
    buf: Vec<f64>,
    // trigonometric tables
    //   "2 pi i / N"
    //   and
    //   "pi i / 2 / N",
    //   where i = 0, 1, ..., N/2-1
    // 0: cos, 1: sin
    mult1s: Vec<[f64; 2]>,
    mult2s: Vec<[f64; 2]>,
}

impl Plan {
    /// Constructor of [`Plan`] struct.
    ///
    /// # Arguments
    /// * `nitems` - size of a signal to be transformed.
    ///
    /// # Returns
    /// A new initialised array.
    pub fn new(nitems: usize) -> Plan {
        // sanitise
        let message = String::from(format!(
            "nitems {} should be 2^(N+2) or 3 x 2^(N+2)",
            nitems
        ));
        if nitems < 4 {
            panic!("{}", message);
        }
        if 0 == nitems % 3 {
            let n: usize = nitems / 3;
            if 0 != n & (n - 1) {
                panic!("{}", message);
            }
        } else if 0 != nitems & (nitems - 1) {
            panic!("{}", message);
        }
        let buf: Vec<f64> = vec![0f64; nitems + 2];
        let mut mult1s: Vec<[f64; 2]> = vec![[0f64; 2]; nitems / 2];
        let mut mult2s: Vec<[f64; 2]> = vec![[0f64; 2]; nitems / 2];
        for i in 0..nitems / 2 {
            let phase: f64 = 2. * std::f64::consts::PI * i as f64 / nitems as f64;
            mult1s[i][0] = phase.cos();
            mult1s[i][1] = phase.sin();
        }
        for i in 0..nitems / 2 {
            let phase: f64 = std::f64::consts::PI * i as f64 / 2. / nitems as f64;
            mult2s[i][0] = phase.cos();
            mult2s[i][1] = phase.sin();
        }
        return Plan {
            nitems,
            buf,
            mult1s,
            mult2s,
        };
    }

    /// Performs fast cosine transform by Makhoul 1980
    ///
    /// # Arguments
    /// * `xs` - input signal
    /// * `ys` - output signal
    ///
    /// # Caveats
    /// * `ys` should be properly allocated beforehand such that the signal size can be stored.
    /// * The input signal should be reordered beforehand (see [`Plan::map_index`]), which is to avoid the unnecessary buffer copies.
    pub fn exec_f(&mut self, xs: &[f64], ys: &mut [f64]) -> () {
        const IS_FORWARD: bool = true;
        let nitems: usize = self.nitems;
        let buf: &mut Vec<f64> = &mut self.buf;
        let mult1s: &Vec<[f64; 2]> = &self.mult1s;
        let mult2s: &Vec<[f64; 2]> = &self.mult2s;
        // step 1: reorder input, should be done by user
        // step 2: fft, followed by fig. 5a
        fft(IS_FORWARD, nitems / 2, 1, mult1s, xs, buf);
        ys[0] = 2. * (buf[0] + buf[1]);
        ys[nitems / 2] = std::f64::consts::SQRT_2 * (buf[0] - buf[1]);
        buf[nitems / 2 + 0] *= 0. + 2.;
        buf[nitems / 2 + 1] *= 0. - 2.;
        for i in 1..nitems / 4 {
            let j: usize = nitems / 2 - i;
            exchange(IS_FORWARD, &mult1s[i], [2 * i, 2 * j], buf);
        }
        // step 3
        for i in 1..nitems / 2 {
            let j: usize = nitems - i;
            let mult_c: f64 = mult2s[i][0];
            let mult_s: f64 = mult2s[i][1];
            ys[i] = 0. + buf[i * 2 + 0] * mult_c + buf[i * 2 + 1] * mult_s;
            ys[j] = 0. - buf[i * 2 + 1] * mult_c + buf[i * 2 + 0] * mult_s;
        }
    }

    /// Performs inverse fast cosine transform by Makhoul 1980
    ///
    /// # Arguments
    /// * `xs` - input signal
    /// * `ys` - output signal
    ///
    /// # Caveats
    /// * `ys` should be properly allocated beforehand such that the signal size can be stored.
    /// * The results are doubled so that it is equivalent to FFTW's DCT 3.
    /// * The output signal should be reordered afterwards (see [`Plan::map_index`]), which is to avoid the unnecessary buffer copies.
    pub fn exec_b(&mut self, xs: &[f64], ys: &mut [f64]) -> () {
        const IS_FORWARD: bool = false;
        let nitems: usize = self.nitems;
        let buf: &mut Vec<f64> = &mut self.buf;
        let mult1s: &Vec<[f64; 2]> = &self.mult1s;
        let mult2s: &Vec<[f64; 2]> = &self.mult2s;
        // step 1: eq. 28
        for i in 1..nitems / 2 {
            let j: usize = nitems - i;
            let mult_c: f64 = mult2s[i][0];
            let mult_s: f64 = mult2s[i][1];
            buf[i * 2 + 0] = xs[i] * mult_c + xs[j] * mult_s;
            buf[i * 2 + 1] = xs[i] * mult_s - xs[j] * mult_c;
        }
        // step 2: fig. 5(b), followed by iFFT
        buf[1] = xs[0] - std::f64::consts::SQRT_2 * xs[nitems / 2];
        buf[0] = xs[0] + std::f64::consts::SQRT_2 * xs[nitems / 2];
        buf[nitems / 2 + 0] *= 0. + 2.;
        buf[nitems / 2 + 1] *= 0. - 2.;
        for i in 1..nitems / 4 {
            let j: usize = nitems / 2 - i;
            exchange(IS_FORWARD, &mult1s[i], [2 * i, 2 * j], buf);
        }
        fft(IS_FORWARD, nitems / 2, 1, mult1s, buf, ys);
        // step 3: reorder output, should be done by user
    }

    /// A helper function to reorder input / output signals to use fast cosine transforms by Makhoul 1980.
    ///
    /// # Arguments
    /// * `index` - index in the user space.
    ///
    /// # Returns
    /// * Index in the DCT space.
    ///
    /// # Reference
    /// Makhoul 1980, Equations 20 and A-1.
    pub fn map_index(&self, index: usize) -> usize {
        return map_index(self.nitems, index);
    }
}

/// Computes discrete Fourier transform, using Cooley-Tukey algorithm.
///
/// * `is_forward` - a flag to specify forward / backward transforms.
/// * `nitems`     - the size of the signal.
/// * `stride`     - stride between two neighbouring elements of the input signal.
/// * `mults`      - trigonometic tables.
/// * `xs`         - input signal.
/// * `ys`         - output signal.
fn fft(
    is_forward: bool,
    nitems: usize,
    stride: usize,
    mults: &Vec<[f64; 2]>,
    xs: &[f64],
    ys: &mut [f64],
) -> () {
    // small sizes, use analytical solutions
    if 1 == nitems {
        ys[0] = xs[0];
        ys[1] = xs[1];
        return;
    }
    if 2 == nitems {
        let x00: f64 = xs[0 * stride * 2 + 0];
        let x01: f64 = xs[0 * stride * 2 + 1];
        let x10: f64 = xs[1 * stride * 2 + 0];
        let x11: f64 = xs[1 * stride * 2 + 1];
        ys[0] = x00 + x10;
        ys[1] = x01 + x11;
        ys[2] = x00 - x10;
        ys[3] = x01 - x11;
        return;
    }
    let sign: f64 = if is_forward { -1. } else { 1. };
    if 3 == nitems {
        let x00: f64 = xs[0 * stride * 2 + 0];
        let x01: f64 = xs[0 * stride * 2 + 1];
        let x10: f64 = xs[1 * stride * 2 + 0];
        let x11: f64 = xs[1 * stride * 2 + 1];
        let x20: f64 = xs[2 * stride * 2 + 0];
        let x21: f64 = xs[2 * stride * 2 + 1];
        // sqrt(3) / 2
        const SQRT3H: f64 = 0.5 * 1.732050807568877293527446341505872367_f64;
        ys[0] = x00 + x10 + x20;
        ys[1] = x01 + x11 + x21;
        ys[2] = (x00 - 0.5 * x10 - 0.5 * x20) - sign * SQRT3H * (x11 - x21);
        ys[3] = (x01 - 0.5 * x11 - 0.5 * x21) + sign * SQRT3H * (x10 - x20);
        ys[4] = (x00 - 0.5 * x10 - 0.5 * x20) + sign * SQRT3H * (x11 - x21);
        ys[5] = (x01 - 0.5 * x11 - 0.5 * x21) - sign * SQRT3H * (x10 - x20);
        return;
    }
    if 4 == nitems {
        let x00: f64 = xs[0 * stride * 2 + 0];
        let x01: f64 = xs[0 * stride * 2 + 1];
        let x10: f64 = xs[1 * stride * 2 + 0];
        let x11: f64 = xs[1 * stride * 2 + 1];
        let x20: f64 = xs[2 * stride * 2 + 0];
        let x21: f64 = xs[2 * stride * 2 + 1];
        let x30: f64 = xs[3 * stride * 2 + 0];
        let x31: f64 = xs[3 * stride * 2 + 1];
        ys[0] = (x00 + x20) + (x10 + x30);
        ys[1] = (x01 + x21) + (x11 + x31);
        ys[2] = (x00 - x20) - sign * (x11 - x31);
        ys[3] = (x01 - x21) + sign * (x10 - x30);
        ys[4] = (x00 + x20) - (x10 + x30);
        ys[5] = (x01 + x21) - (x11 + x31);
        ys[6] = (x00 - x20) + sign * (x11 - x31);
        ys[7] = (x01 - x21) - sign * (x10 - x30);
        return;
    }
    // general sizes, use recurrent relations
    // an FFT is divided into two sub-FFTs: even-indexed and odd-indexed
    // e.g., N = 8
    //   even-indexed
    //   x[0], x[2], x[4], x[6], ... -> y[0], y[1], y[2], y[3], ...
    //   odd -indexed
    //   x[1], x[3], x[5], x[7], ... -> y[4], y[5], y[6], y[7], ...
    fft(
        is_forward,
        nitems / 2,
        stride * 2,
        mults,
        &xs[0..],
        &mut ys[0..],
    );
    fft(
        is_forward,
        nitems / 2,
        stride * 2,
        mults,
        &xs[stride * 2..],
        &mut ys[nitems..],
    );
    // combine even and odd results
    for i in 0..nitems / 2 {
        let j: usize = i + nitems / 2;
        // multiplier, equivalent to
        //   mult_c = cos( 2 * pi * i / nitems )
        //   mult_s = sin( 2 * pi * i / nitems )
        let mult_c: f64 = mults[stride * 2 * i][0];
        let mult_s: f64 = mults[stride * 2 * i][1];
        // exchange between i and j
        let y00: f64 = ys[i * 2 + 0];
        let y01: f64 = ys[i * 2 + 1];
        let y10: f64 = ys[j * 2 + 0];
        let y11: f64 = ys[j * 2 + 1];
        ys[i * 2 + 0] = y00 + y10 * mult_c - sign * y11 * mult_s;
        ys[i * 2 + 1] = y01 + y11 * mult_c + sign * y10 * mult_s;
        ys[j * 2 + 0] = y00 - y10 * mult_c + sign * y11 * mult_s;
        ys[j * 2 + 1] = y01 - y11 * mult_c - sign * y10 * mult_s;
    }
}

/// Exchanges information between i-th and j-th elemenst
///
/// # Reference
/// Makhoul 1980, Figure 5.
///
/// * `is_forward` - a flag to specify forward / backward transforms.
/// * `mult`       - an element of the trigonometric table.
/// * `indices`    - `i` and `j`.
/// * `buf`        - a buffer whose i-th and j-th elements are swapped.
fn exchange(is_forward: bool, mult: &[f64; 2], indices: [usize; 2], buf: &mut [f64]) -> () {
    let sign: f64 = if is_forward { -1. } else { 1. };
    let x0: f64 = buf[indices[0] + 0];
    let x1: f64 = buf[indices[0] + 1];
    let y0: f64 = buf[indices[1] + 0];
    let y1: f64 = buf[indices[1] + 1];
    let a0: f64 = x0 + y0;
    let s0: f64 = x0 - y0;
    let a1: f64 = x1 + y1;
    let s1: f64 = x1 - y1;
    buf[indices[0] + 0] = 0. + a0 - sign * mult[0] * a1 - mult[1] * s0;
    buf[indices[0] + 1] = 0. + s1 + sign * mult[0] * s0 - mult[1] * a1;
    buf[indices[1] + 0] = 0. + a0 + sign * mult[0] * a1 + mult[1] * s0;
    buf[indices[1] + 1] = 0. - s1 + sign * mult[0] * s0 - mult[1] * a1;
}

/// A helper function to reorder input / output signals to use fast cosine transforms by Makhoul 1980.
///
/// # Arguments
/// * `nitems` - signal length.
/// * `i`      - index in the user space.
///
/// # Returns
/// * Index in the DCT space.
///
/// # Reference
/// Makhoul 1980, Equations 20 and A-1.
fn map_index(nitems: usize, i: usize) -> usize {
    if 0 == i % 2 {
        return i / 2;
    } else {
        return nitems - (i - 1) / 2 - 1;
    }
}

#[cfg(test)]
mod test_map_index {
    #[test]
    fn case1() -> () {
        let nitems: usize = 8;
        let result: Vec<usize> = (0..nitems).map(|n| super::map_index(nitems, n)).collect();
        let answer = [0, 7, 1, 6, 2, 5, 3, 4];
        for n in 0..nitems {
            assert_eq!(result[n], answer[n]);
        }
    }
    #[test]
    fn case2() -> () {
        let nitems: usize = 12;
        let result: Vec<usize> = (0..nitems).map(|n| super::map_index(nitems, n)).collect();
        let answer = [0, 11, 1, 10, 2, 9, 3, 8, 4, 7, 5, 6];
        for n in 0..nitems {
            assert_eq!(result[n], answer[n]);
        }
    }
}

#[cfg(test)]
mod test_dct {
    /// A naive DCT-2 implementation.
    fn dct_f(nitems: usize, xs: &[f64], ys: &mut [f64]) -> () {
        for j in 0..nitems {
            for i in 0..nitems {
                let phase: f64 =
                    std::f64::consts::PI * (2 * i + 1) as f64 * j as f64 / (2. * nitems as f64);
                ys[j] += 2. * xs[i] * phase.cos();
            }
        }
    }
    /// A naive DCT-3 implementation.
    fn dct_b(nitems: usize, xs: &[f64], ys: &mut [f64]) -> () {
        for j in 0..nitems {
            ys[j] = xs[0];
            for i in 1..nitems {
                let phase: f64 =
                    std::f64::consts::PI * (2 * j + 1) as f64 * i as f64 / (2. * nitems as f64);
                ys[j] += 2. * xs[i] * phase.cos();
            }
        }
    }
    #[test]
    /// Performs naive DCT-2 and FCT-2 and compares results.
    fn f() -> () {
        let nitemss = [4, 8, 12, 16, 32, 48, 64, 96, 128, 192, 256];
        for nitems in nitemss {
            let input0: Vec<f64> = (0..nitems).map(|x| x as f64 / nitems as f64).collect();
            let mut input1: Vec<f64> = vec![0.; nitems];
            let mut output0: Vec<f64> = vec![0.; nitems];
            let mut output1: Vec<f64> = vec![0.; nitems];
            let mut plan = super::Plan::new(nitems);
            // reorder input buffer
            for i in 0..nitems {
                let ii: usize = plan.map_index(i);
                input1[ii] = input0[i];
            }
            dct_f(nitems, &input0, &mut output0);
            plan.exec_f(&input1, &mut output1);
            let mut dif: f64 = 0.;
            for n in 0..nitems {
                dif += (output0[n] - output1[n]).abs() / nitems as f64;
            }
            assert!(dif < 1e-14 * nitems as f64);
        }
    }
    #[test]
    /// Performs naive DCT-3 and FCT-3 and compares results.
    fn b() -> () {
        let nitemss = [4, 8, 12, 16, 32, 48, 64, 96, 128, 192, 256];
        for nitems in nitemss {
            let input: Vec<f64> = (0..nitems).map(|x| x as f64 / nitems as f64).collect();
            let mut output0: Vec<f64> = vec![0.; nitems];
            let mut output1: Vec<f64> = vec![0.; nitems];
            let mut output2: Vec<f64> = vec![0.; nitems];
            let mut plan = super::Plan::new(nitems);
            dct_b(nitems, &input, &mut output0);
            plan.exec_b(&input, &mut output1);
            // reorder output buffer
            for i in 0..nitems {
                let ii: usize = plan.map_index(i);
                output2[i] = output1[ii];
            }
            let mut dif: f64 = 0.;
            for n in 0..nitems {
                dif += (output0[n] - output2[n]).abs() / nitems as f64;
            }
            assert!(dif < 1e-14 * nitems as f64);
        }
    }
}
