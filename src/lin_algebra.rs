use crate::matrix::Matrix;
use crate::vector::Vector;

pub struct LinAlgebra {}

impl LinAlgebra {
    pub fn mult_matrix_by_vector(matrix: &Matrix, old_vector: &Vector) -> Vector {
        if old_vector.n_rows() != matrix.cols {
            panic!("Invalid multiplication due to mismatched dimensions.")
        }
        let mut new_vector: Vector = Vector::zeros(matrix.rows);
        for i in 0..new_vector.n_rows() {
            let mut val = 0.0;
            for j in 0..matrix.cols {
                val += old_vector.get_value(j) * matrix.value[i][j];
                new_vector.set_value(i, val);
                //  new_vector.value[i] += vector.value[j] * matrix.value[i][j];
            }
        }
        new_vector
    }

    //Solve linear equation sistem Ax = b using crout method
    // The LU Decomposition always work correctly when the matrices is
    // symetric positive-definite.
    pub fn crout_method_solve(a_matrix: &Matrix, b_vector: &Vector) {
        let n = a_matrix.rows;
        let mut l_matrix: Matrix = Matrix::zeros(n, n);
        let mut u_matrix: Matrix = Matrix::identity(n);

        //get low matrix
        for j in 0..n {
            for i in j..n {
                let mut sum = 0.0;
                for k in 0..j {
                    sum += l_matrix.value[i][k] * u_matrix.value[k][j];
                }
                l_matrix.value[i][j] = a_matrix.value[i][j] - sum;
            }
            //get upper matrix
            for i in j..n {
                let mut sum = 0.0;
                for k in 0..j {
                    sum += l_matrix.value[j][k] * u_matrix.value[k][i];
                }

                if l_matrix.value[j][j] == 0.0 {
                    panic!("Matrix is singular, LU decomposition failed.");
                }

                u_matrix.value[j][i] = (a_matrix.value[j][i] - sum) / l_matrix.value[j][j];
            }
        }

        println!("L in lin_algerbra = {:?}", l_matrix);
        println!("U in lin_algerbra = {:?}", u_matrix);

        let y_vector: Vector = LinAlgebra::foward_substitution(&l_matrix, &b_vector);

        let x_vector: Vector = LinAlgebra::back_substitution(&u_matrix, &y_vector);
    }

    fn foward_substitution(l_matrix: &Matrix, b_vector: &Vector) -> Vector {
        let n = l_matrix.rows;
        let mut y_vector: Vector = Vector::zeros(n);
        y_vector.set_value(0, b_vector.get_value(0) / l_matrix.value[0][0]);
        //y_vector.value[0] = b_vector.value[0] / l_matrix.value[0][0];

        for i in 1..(n) {
            let mut sum = 0.0;
            for j in 0..(n - 1) {
                sum += y_vector.get_value(j) * l_matrix.value[i][j];
                //sum += y_vector.value[j] * l_matrix.value[i][j];
            }
            y_vector.set_value(i, (b_vector.get_value(i) - sum) / l_matrix.value[i][i]);
            //y_vector.value[i] = (b_vector.value[i] - sum) / l_matrix.value[i][i];
        }
        println!(" y_vector in lin_algebra = {:?}", y_vector);
        y_vector
    }

    fn back_substitution(u_matrix: &Matrix, y_vector: &Vector) -> Vector {
        let n = u_matrix.rows;
        let mut x_vector: Vector = Vector::zeros(n);

        x_vector.set_value(
            n - 1,
            y_vector.get_value(n - 1) / u_matrix.value[n - 1][n - 1],
        );
        //x_vector.value[n-1] = y_vector.value[n-1] / u_matrix.value[n-1][n-1];
        //reverse loop "for" counted step by step
        for i in (0..=n - 1).rev().step_by(1) {
            let mut sum = 0.0;
            for j in i + 1..n {
                sum += x_vector.get_value(j) * u_matrix.value[i][j];
                //sum += x_vector.value[j] * u_matrix.value[i][j];
            }
            x_vector.set_value(i, (y_vector.get_value(i) - sum) / u_matrix.value[i][i]);
            // x_vector.value[i] = (y_vector.value[i] - sum) / u_matrix.value[i][i];
        }
        println!(" x_vector in lin_algebra = {:?}", x_vector);
        x_vector
    }
}
