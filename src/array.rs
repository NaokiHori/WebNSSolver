pub struct Array {
    nx: usize,
    ny: usize,
    data: Vec<f64>,
}

impl Array {
    pub fn new(nx: usize, ny: usize) -> Array {
        let data: Vec<f64> = vec![0.; nx * ny];
        Array { nx, ny, data }
    }

    pub fn get_size(&self) -> [usize; 2] {
        [self.nx, self.ny]
    }

    pub fn get_data(&self) -> &Vec<f64> {
        &self.data
    }

    pub fn get_data_mut(&mut self) -> &mut Vec<f64> {
        &mut self.data
    }

    pub fn get_minmax(&self) -> [f64; 2] {
        let data: &Vec<f64> = &self.data;
        let mut min: f64 = f64::MAX;
        let mut max: f64 = f64::MIN;
        for &val in data.iter() {
            let tmp = val;
            min = if tmp < min { tmp } else { min };
            max = if max < tmp { tmp } else { max };
        }
        [min, max]
    }

    pub fn get_ext(&self) -> f64 {
        let minmax: [f64; 2] = self.get_minmax();
        let absmin = minmax[0].abs();
        let absmax = minmax[1].abs();
        absmin.max(absmax)
    }
}

impl std::ops::Index<usize> for Array {
    type Output = [f64];
    fn index(&self, j: usize) -> &Self::Output {
        &self.data[j * self.nx..(j + 1) * self.nx]
    }
}

impl std::ops::IndexMut<usize> for Array {
    fn index_mut(&mut self, j: usize) -> &mut Self::Output {
        &mut self.data[j * self.nx..(j + 1) * self.nx]
    }
}
