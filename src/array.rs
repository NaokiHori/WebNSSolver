#![deny(missing_docs)]

//! Represents two-dimensional arrays.
//!
//! # Overview
//! Internally I use one-dimensional vectors [`Vec<f64>`] to store two-dimensional datasets.
//! Describing finite-difference schemes using such vectors is, however, tedious since the size of the arrays (i.e. strides) are always to be written explicitly, e.g.:
//! ```rust
//! data[(j + 1) * (nx + 1) + (i - 1)]
//! ```
//! This module defines two-dimensional index access such that each element can be accessed as
//! ```rust
//! data[j + 1][i - 1]
//! ```
//! by defining traits.
//! Note that two traits are to be defined to allow immutable (right-hand side, via [`std::ops::Index`]) and mutable (left-hand side, via [`std::ops::IndexMut`]) accesses.

/// Stores the size of an array (whose type is [`f64`]) and its dataset.
pub struct Array {
    /// Number of grid points in the x direction.
    nx: usize,
    /// Number of grid points in the y direction.
    ny: usize,
    /// Dataset.
    data: Vec<f64>,
}

impl Array {
    /// Constructor of [`Array`] struct.
    ///
    /// # Arguments
    /// * `nx` : number of grid points in the x direction.
    /// * `ny` : number of grid points in the y direction.
    ///
    /// # Returns
    /// A new initialised array.
    pub fn new(nx: usize, ny: usize) -> Array {
        let data: Vec<f64> = vec![0.; nx * ny];
        return Array { nx, ny, data };
    }

    /// Getter of the array sizes: [nx](Array::nx) and [ny](Array::ny).
    ///
    /// # Returns
    /// Array sizes, [`nx`, `ny`].
    #[allow(dead_code)]
    pub fn get_size(&self) -> [usize; 2] {
        return [self.nx, self.ny];
    }

    /// Getter of the dataset: [data](Array::data), immutably.
    ///
    /// # Returns
    /// Immutable reference to the dataset.
    #[allow(dead_code)]
    pub fn get_data(&self) -> &Vec<f64> {
        return &self.data;
    }

    /// Getter of the dataset: [data](Array::data), mutably.
    ///
    /// # Returns
    /// Mutable reference to the dataset.
    #[allow(dead_code)]
    pub fn get_data_mut(&mut self) -> &mut Vec<f64> {
        return &mut self.data;
    }

    /// Computes the minimum and the maximum values of the dataset.
    ///
    /// # Returns
    /// The minimum and the maximum values: [`min`, `max`].
    #[allow(dead_code)]
    pub fn get_minmax(&self) -> [f64; 2] {
        let data: &Vec<f64> = &self.data;
        let mut min: f64 = f64::MAX;
        let mut max: f64 = f64::MIN;
        for &val in data.iter() {
            let tmp = val;
            min = if tmp < min { tmp } else { min };
            max = if max < tmp { tmp } else { max };
        }
        return [min, max];
    }

    /// Computes the extreme value `max(|data|)` of the dataset.
    ///
    /// # Returns
    /// The extreme value.
    #[allow(dead_code)]
    pub fn get_ext(&self) -> f64 {
        let minmax: [f64; 2] = self.get_minmax();
        let absmin = minmax[0].abs();
        let absmax = minmax[1].abs();
        return absmin.max(absmax);
    }
}

/// Defines immutable index access.
impl std::ops::Index<usize> for Array {
    type Output = [f64];
    #[inline(always)]
    fn index(&self, j: usize) -> &Self::Output {
        return &self.data[j * self.nx..(j + 1) * self.nx];
    }
}

/// Defines mutable index access.
impl std::ops::IndexMut<usize> for Array {
    #[inline(always)]
    fn index_mut(&mut self, j: usize) -> &mut Self::Output {
        return &mut self.data[j * self.nx..(j + 1) * self.nx];
    }
}
