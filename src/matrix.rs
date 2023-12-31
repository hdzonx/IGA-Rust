use std::fmt::{Debug, Formatter, Result};

//#[derive(Debug)]
#[derive(PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    value: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn n_rows(&self) -> usize {
        self.rows
    }
    pub fn n_cols(&self) -> usize {
        self.cols
    }

    pub fn set_value(&mut self, row: usize, col: usize, data: f64) {
        self.value[row][col] = data;
    }

    pub fn get_value(&self, row: usize, col: usize) -> f64 {
        self.value[row][col]
    }

    pub fn new(rows: usize, cols: usize) -> Matrix {
        Matrix {
            rows,
            cols,
            value: vec![vec![0.0; cols]; rows],
        }
    }

    pub fn identity(order: usize) -> Matrix {
        let mut mat = Matrix::new(order, order);

        for i in 0..order {
            for j in 0..order {
                if i == j {
                    mat.value[i][j] = 1.0;
                }
            }
        }
        mat
    }

    pub fn add_matrices(&self, other: &Matrix) -> Matrix {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("add matrix of different dimension is not possible")
        }

        let mut mat = Matrix::new(self.rows, self.cols);

        for i in 0..self.rows {
            for j in 0..self.cols {
                mat.value[i][j] = self.value[i][j] + other.value[i][j];
            }
        }

        mat
    }

    pub fn subtraction_matrices(&self, other: &Matrix) -> Matrix {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("subtract matrix of different dimension is not possible")
        }

        let mut mat = Matrix::new(self.rows, self.cols);

        for i in 0..self.rows {
            for j in 0..self.cols {
                mat.value[i][j] = self.value[i][j] - other.value[i][j];
            }
        }

        mat
    }

    pub fn transpose(&self) -> Matrix {
        let mut mat = Matrix::new(self.rows, self.cols);

        for i in 0..self.rows {
            for j in 0..self.cols {
                mat.value[j][i] = self.value[i][j]
            }
        }
        mat
    }

    pub fn multiply(&self, other: &Matrix) -> Matrix {
        if other.rows != self.cols {
            panic!("multiply matrices of incompatible dimension is not possible");
        }

        let mut mat = Matrix::new(self.rows, other.cols);

        for i in 0..self.rows {
            for j in 0..self.cols {
                let mut sum = 0.0;
                for n in 0..self.cols {
                    sum += self.value[i][n] * other.value[n][j]
                }
                mat.value[i][j] = sum;
            }
        }
        mat
    }

    pub fn dot_mult(&self, other: &Matrix) -> Matrix {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("dot multiplication of imcopatible dimension is not possible");
        }

        let mut mat = Matrix::new(self.rows, self.cols);

        for i in 0..self.rows {
            for j in 0..self.cols {
                mat.value[i][j] = self.value[i][j] * other.value[i][j];
            }
        }
        mat
    }
}

impl Debug for Matrix {
    fn fmt(&self, f: &mut Formatter) -> Result {
        write!(
            f,
            "Matrix {{\n{}\n}}",
            (&self.value)
                .into_iter()
                .map(|row| "  ".to_string()
                    + &row
                        .into_iter()
                        .map(|value| value.to_string())
                        .collect::<Vec<String>>()
                        .join(" "))
                .collect::<Vec<String>>()
                .join("\n")
        )
    }
}
