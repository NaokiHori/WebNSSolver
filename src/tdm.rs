#![deny(missing_docs)]

//! Solves tri-diagonal linear systems by means of the Thomas algorithm.

/// Stores state and the internal buffers.
pub struct Plan {
    /// size of the system.
    nitems: usize,
    /// stride between two neighbouring input / output buffer.
    stride: usize,
    /// lower-diagonal part of the system.
    pub l: Vec<f64>,
    /// centre-diagonal part of the system.
    pub c: Vec<f64>,
    /// upper-diagonal part of the system.
    pub u: Vec<f64>,
    /// buffer for internal use.
    v: Vec<f64>,
}

impl Plan {
    /// Constructor of [`Plan`] struct.
    ///
    /// # Arguments
    /// * `nitems` - size of the system.
    /// * `stride` - stride between two neighbouring input / output buffer.
    ///
    /// # Returns
    /// A new initialised array.
    pub fn new(nitems: usize, stride: usize) -> Plan {
        let l: Vec<f64> = vec![0.; nitems];
        let c: Vec<f64> = vec![0.; nitems];
        let u: Vec<f64> = vec![0.; nitems];
        let v: Vec<f64> = vec![0.; nitems];
        return Plan {
            nitems,
            stride,
            l,
            c,
            u,
            v,
        };
    }

    /// Solves tri-diagonal linear systems by means of the Thomas algorithm.
    ///
    /// Note that the rank of the linear system can be `nitems-1`.
    pub fn solve(&mut self, q: &mut [f64]) -> () {
        // Thomas algorithm for a strided input / output (q)
        let nitems: usize = self.nitems;
        let stride: usize = self.stride;
        let l: &Vec<f64> = &self.l;
        let c: &Vec<f64> = &self.c;
        let u: &Vec<f64> = &self.u;
        let v: &mut Vec<f64> = &mut self.v;
        // forward
        v[0] = u[0] / c[0];
        q[0] = q[0] / c[0];
        for i in 1..nitems - 1 {
            let val: f64 = 1. / (c[i] - l[i] * v[i - 1]);
            v[i] = val * u[i];
            q[i * stride] = val * (q[i * stride] - l[i] * q[(i - 1) * stride]);
        }
        let val: f64 = c[nitems - 1] - l[nitems - 1] * v[nitems - 2];
        // the last row
        if val.abs() > f64::EPSILON {
            q[(nitems - 1) * stride] =
                1. / val * (q[(nitems - 1) * stride] - l[nitems - 1] * q[(nitems - 2) * stride]);
        } else {
            // eigen-value degenerates, caused by Neumann BCs.
            q[(nitems - 1) * stride] = 0.;
        }
        // backward
        for i in (0..nitems - 1).rev() {
            q[i * stride] -= v[i] * q[(i + 1) * stride];
        }
    }
}

#[cfg(test)]
mod test {
    const EPS: f64 = 1e-8;
    #[test]
    fn case1() -> () {
        let nitems: usize = 3;
        let stride: usize = 1;
        let mut plan = super::Plan::new(nitems, stride);
        plan.l[1] = 1.;
        plan.l[2] = 1.;
        plan.c[0] = -2.;
        plan.c[1] = -2.;
        plan.c[2] = -2.;
        plan.u[0] = 1.;
        plan.u[1] = 1.;
        let mut q: Vec<f64> = vec![1., 2., 3.];
        let ans: Vec<f64> = vec![-2.5, -4., -3.5];
        plan.solve(&mut q[stride - 1..]);
        for n in 0..nitems {
            assert!(ans[n] - EPS < q[n]);
            assert!(q[n] < ans[n] + EPS);
        }
    }
    #[test]
    fn case2() -> () {
        let nitems: usize = 3;
        let stride: usize = 2;
        let mut plan = super::Plan::new(nitems, stride);
        plan.l[1] = 1.;
        plan.l[2] = 1.;
        plan.c[0] = -2.;
        plan.c[1] = -2.;
        plan.c[2] = -2.;
        plan.u[0] = 1.;
        plan.u[1] = 1.;
        let mut q: Vec<f64> = vec![0., 1., 0., 2., 0., 3.];
        let ans: Vec<f64> = vec![-2.5, -4., -3.5];
        plan.solve(&mut q[stride - 1..]);
        for n in 0..nitems {
            assert!(ans[n] - EPS < q[n * stride + stride - 1]);
            assert!(q[n * stride + stride - 1] < ans[n] + EPS);
        }
    }
}
