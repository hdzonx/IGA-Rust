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
        let mut mat = Matrix::new(self.cols, self.rows);

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
            for j in 0..other.cols {
                let mut sum = 0.0;
                for n in 0..self.cols {
                    let a = self.value[i][n];
                    let b = other.value[n][j];
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

    pub fn mult_scalar(&self, scalar: f64) -> Matrix {
        let mut mat: Matrix = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                mat.value[i][j] = scalar * self.value[i][j];
            }
        }
        mat
    }
}

impl Debug for Matrix {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        let result = write!(
            f,
            "Matrix {{\n{}\n}}",
            self.value
                .iter()
                .map(|row| {
                    format!(
                        "  {}",
                        row.iter()
                            .map(|value| value.to_string())
                            .collect::<Vec<_>>()
                            .join(" ")
                    )
                })
                .collect::<Vec<_>>()
                .join("\n")
        );
        result
    }
}

#[cfg(test)]
mod tests {

    use crate::matrix::Matrix;

    #[test]
    fn multiply_test() {
        let mut b = Matrix::new(2, 5);
        b.set_value(0, 0, 1.0);
        b.set_value(0, 1, 2.0);
        b.set_value(0, 2, 3.0);
        b.set_value(0, 3, 4.0);
        b.set_value(0, 4, 5.0);

        b.set_value(1, 0, 6.0);
        b.set_value(1, 1, 7.0);
        b.set_value(1, 2, 8.0);
        b.set_value(1, 3, 9.0);
        b.set_value(1, 4, 10.0);

        let mut d = Matrix::new(2, 2);
        d.set_value(0, 0, 2.0);
        d.set_value(0, 1, 6.0);

        d.set_value(1, 0, 9.0);
        d.set_value(1, 1, 7.0);

        let db = d.multiply(&b);
        println!("{:?}", db);

        let tranpose_b = b.transpose();

        let k = tranpose_b.multiply(&db);
        println!("{:?}", k);
    }

    #[test]
    fn multiply_01() {
        let mut b_matrix = Matrix::new(3, 6);
        //First line
        b_matrix.set_value(0, 0, -1.5);
        b_matrix.set_value(0, 1, 0.0);
        b_matrix.set_value(0, 2, -3.5);
        b_matrix.set_value(0, 3, 0.0);
        b_matrix.set_value(0, 4, 5.0);
        b_matrix.set_value(0, 5, 0.0);
        //Second line
        b_matrix.set_value(1, 0, 0.0);
        b_matrix.set_value(1, 1, 5.5);
        b_matrix.set_value(1, 2, 0.0);
        b_matrix.set_value(1, 3, -3.0);
        b_matrix.set_value(1, 4, 0.0);
        b_matrix.set_value(1, 5, -2.5);
        //Third line
        b_matrix.set_value(2, 0, 5.5);
        b_matrix.set_value(2, 1, -1.5);
        b_matrix.set_value(2, 2, -3.0);
        b_matrix.set_value(2, 3, 3.5);
        b_matrix.set_value(2, 4, -2.5);
        b_matrix.set_value(2, 5, 5.0);

        let det_jacobian = 0.042105263;

        let new_b = b_matrix.mult_scalar(det_jacobian);

        // println!("{:?}", new_b);

        let new_b_transpose = new_b.transpose();

        // println!("{:?}", new_b_transpose);

        let mut constitutive_matrix = Matrix::new(3, 3);
        constitutive_matrix.set_value(0, 0, 220800.0);
        constitutive_matrix.set_value(0, 1, 55200.0);
        constitutive_matrix.set_value(0, 2, 0.0);
        constitutive_matrix.set_value(1, 0, 55200.0);
        constitutive_matrix.set_value(1, 1, 220800.0);
        constitutive_matrix.set_value(1, 1, 0.0);
        constitutive_matrix.set_value(2, 0, 0.0);
        constitutive_matrix.set_value(2, 1, 0.0);
        constitutive_matrix.set_value(2, 2, 82800.0);

        let db_matrix = constitutive_matrix.multiply(&new_b);
        // println!("{:?}", db_matrix);

        let k_matrix = new_b_transpose.multiply(&db_matrix);
        println!("{:?}", k_matrix);
    }
}
