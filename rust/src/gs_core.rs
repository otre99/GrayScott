use std::fs::File;
use std::io::Write;
use std::mem::swap;
use std::slice::from_raw_parts;

pub(crate) struct GrayScottSolver {
    dim: usize,
    mu: f64,
    mv: f64,
    f: f64,
    k: f64,
    aux_u: Vec<f64>,
    aux_v: Vec<f64>,
    data_u: Vec<f64>,
    data_v: Vec<f64>,
    iter_count: usize,
}

impl GrayScottSolver {
    pub(crate) fn new(dim: usize, mu: f64, mv: f64, f: f64, k: f64) -> Self {
        let n = dim * dim;
        let mut gs_instance = Self {
            dim,
            mu,
            mv,
            f,
            k,
            aux_u: vec![0.0; n],
            aux_v: vec![0.0; n],
            data_u: vec![0.0; n],
            data_v: vec![0.0; n],
            iter_count: 0,
        };

        gs_instance.init_ic();
        gs_instance
    }

    fn init_ic(&mut self) {
        let r = self.dim / (self.dim / 16);
        let r1 = self.dim / 2 - r;
        let r2 = self.dim / 2 + r;
        for row in 0..self.dim {
            let i = row * self.dim;

            let cond0 = row >= r1 && row < r2;
            for col in 0..self.dim {
                if cond0 && col >= r1 && col < r2 {
                    self.data_u[i + col] = 0.5;
                    self.data_v[i + col] = 0.25;
                } else {
                    self.data_u[i + col] = 1.0;
                    self.data_v[i + col] = 0.0;
                }
            }
        }
    }

    pub(crate) fn step_seq(&mut self) {
        for i in 0..self.dim {
            self.solve_one_row(i)
        }

        swap(&mut self.aux_u, &mut self.data_u);
        swap(&mut self.aux_v, &mut self.data_v);

        self.iter_count += 1;
    }

    fn solve_one_row(&mut self, row: usize) {
        let (tu, bu): (&[f64], &[f64]);
        let (tv, bv): (&[f64], &[f64]);
        let l2 = row * self.dim;
        let l3 = l2 + self.dim;

        let vv = &mut self.aux_v[l2..l3];
        let uu = &mut self.aux_u[l2..l3];
        let v = &self.data_v[l2..l3];
        let u = &self.data_u[l2..l3];

        if row > 0 && row < self.dim - 1 {
            tu = &self.data_u[(l2 - self.dim)..l2];
            bu = &self.data_u[l3..(l3 + self.dim)];
            tv = &self.data_v[(l2 - self.dim)..l2];
            bv = &self.data_v[l3..(l3 + self.dim)];
        } else if row == 0 {
            let n = self.dim * self.dim;
            tu = &self.data_u[(n - self.dim)..n];
            bu = &self.data_u[l3..(l3 + self.dim)];
            tv = &self.data_v[(n - self.dim)..n];
            bv = &self.data_v[l3..(l3 + self.dim)];
        } else {
            tu = &self.data_u[(l2 - self.dim)..l2];
            bu = &self.data_u[0..self.dim];
            tv = &self.data_v[(l2 - self.dim)..l2];
            bv = &self.data_v[0..self.dim];
        }

        let mut uv2 = u[0] * v[0] * v[0];
        uu[0] = u[0] + self.mu * (u[self.dim - 1] + u[1] + tu[0] + bu[0] - 4.0 * u[0]) - uv2
            + self.f * (1.0 - u[0]);
        vv[0] = v[0] + self.mv * (v[self.dim - 1] + v[1] + tv[0] + bv[0] - 4.0 * v[0]) + uv2
            - (self.f + self.k) * v[0];
        for i in 1..(self.dim - 1) {
            uv2 = u[i] * v[i] * v[i];
            uu[i] = u[i] + self.mu * (u[i - 1] + u[i + 1] + tu[i] + bu[i] - 4.0 * u[i]) - uv2
                + self.f * (1.0 - u[i]);
            vv[i] = v[i] + self.mv * (v[i - 1] + v[i + 1] + tv[i] + bv[i] - 4.0 * v[i]) + uv2
                - (self.f + self.k) * v[i];
        }
        let i = self.dim - 1;
        uv2 = u[i] * v[i] * v[i];
        uu[i] = u[i] + self.mu * (u[i - 1] + u[0] + tu[i] + bu[i] - 4.0 * u[i]) - uv2
            + self.f * (1.0 - u[i]);
        vv[i] = v[i] + self.mv * (v[i - 1] + v[0] + tv[i] + bv[i] - 4.0 * v[i]) + uv2
            - (self.f + self.k) * v[i];
    }

    pub fn save_solution(&mut self, ofolder: &str) {
        let mut fname = String::new();
        std::fmt::write(
            &mut fname,
            format_args!("{}/udata_iter{:08}", ofolder, self.iter_count),
        )
        .unwrap();
        self.save_to_file(&self.data_u, fname.as_str());

        fname.clear();
        std::fmt::write(
            &mut fname,
            format_args!("{}/vdata_iter{:08}", ofolder, self.iter_count),
        )
        .unwrap();
        self.save_to_file(&self.data_v, fname.as_str());
    }

    fn save_to_file(&self, data: &Vec<f64>, fname: &str) {
        let dsize = data.len() * std::mem::size_of::<f64>();
        let buffer: &[u8] = unsafe { from_raw_parts(data.as_ptr() as *const u8, dsize) };
        let mut file = match File::create(fname) {
            Ok(file) => file,
            Err(msg) => {
                println!("Error opening file {} -> {}", fname, msg);
                return;
            }
        };
        file.write(buffer).unwrap();
    }
}
