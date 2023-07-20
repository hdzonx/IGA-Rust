use crate::matrix::Matrix;
use crate::vector::Vector;

pub struct LinAlgebra {

}

impl LinAlgebra {


    pub fn mult_matrix_by_vector(matrix: &Matrix,vector: &Vector) -> Vector {
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


    pub fn crout_method(matrix: &Matrix, vector: &Vector){
        let mut crout_vec: Vector = Vector::zeros(matrix.rows);
        let mut crout_mat:Matrix = Matrix::zeros(matrix.rows, matrix.cols);
        let permut_vec:Vector = Vector::zeros(vector.rows);

    }
}