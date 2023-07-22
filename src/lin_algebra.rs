use crate::matrix::Matrix;
use crate::vector::Vector;

pub struct LinAlgebra {}

impl LinAlgebra {
    pub fn mult_matrix_by_vector(matrix: &Matrix, vector: &Vector) -> Vector {
        if vector.rows != matrix.cols {
            panic!("Invalid multiplication due to mismatched dimensions.")
        }
        let mut new_vector: Vector = Vector::zeros(matrix.rows);
        for i in 0..new_vector.rows {
            for j in 0..matrix.cols {
                new_vector.value[i] += vector.value[j] * matrix.value[i][j];
            }
        }
        new_vector
    }

    pub fn crout_method(matrix: &Matrix, vector: &Vector) {
        //let mut crout_vec: Vector = Vector::zeros(matrix.rows);
        //let mut crout_mat: Matrix = Matrix::zeros(matrix.rows, matrix.cols);
        let permut_vec: Vector = Vector::zeros(vector.rows);
    }

    fn lu_decompostion(matrix: &Matrix) {
        if matrix.rows != matrix.cols {
            panic!("matrix for LU  decomposition must be square")
        }

        let mut lu_matrix: Matrix = Matrix::zeros(matrix.rows, matrix.cols);
        let mut permut_vector: Vector = Vector::zeros(matrix.rows);

        for i in 0..matrix.rows {
            for j in 0..matrix.cols {
                lu_matrix.set_value(i, j, matrix.get_value(i, j));
            }
        }

        let n = lu_matrix.rows - 1;
        let mut t: f64 = 0.0;

        for i in 0..lu_matrix.rows {
            permut_vector.set_value(i, t);
            t += 1.0;
        }

        for i in 0..n {
            let mut z: usize = 0;
            z = i;

            let mut max_val: f64 = lu_matrix.get_value(i, i).abs();

            for j in i..n {
                if lu_matrix.get_value(j, i) > max_val {
                    max_val = lu_matrix.get_value(j, i);
                    z = j;
                }
            }

            if z != i {





            }
        }
    }
}
