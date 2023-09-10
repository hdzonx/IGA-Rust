use crate::matrix::Matrix;
use crate::vector::Vector;

pub struct LinAlgebra {}

impl LinAlgebra {

    pub fn mult_matrix_by_vector(matrix: &Matrix, old_vector: &Vector) -> Vector {
        if old_vector.n_rows() != matrix.n_cols() {
            panic!("Invalid multiplication due to mismatched dimensions.")
        }
        let n = matrix.n_cols();
        let mut new_vector: Vector = Vector::zeros(n);
        for i in 0..new_vector.n_rows() {
            let mut val = 0.0;
            for j in 0..n {
                val += old_vector.get_value(j) * matrix.get_value(i, j);
                new_vector.set_value(i, val);
            }
        }
        new_vector
    }

    //Solve linear equation sistem Ax = b using crout method
    // The LU Decomposition always work correctly when the matrices is
    // symetric positive-definite.
    pub fn crout_method_solve(a_matrix: &Matrix, b_vector: &Vector) ->  Vector{
        let n = a_matrix.n_cols();
        let mut l_matrix: Matrix = Matrix::zeros(n, n);
        let mut u_matrix: Matrix = Matrix::identity(n);

        //get low matrix
        for j in 0..n {
            for i in j..n {
                let mut sum = 0.0;
                for k in 0..j {
                    sum += l_matrix.get_value(i,k) * u_matrix.get_value(k,j);
                }
                l_matrix.set_value(i,j,a_matrix.get_value(i,j) - sum);
            }
            //get upper matrix
            for i in j..n {
                let mut sum = 0.0;
                for k in 0..j {
                    sum += l_matrix.get_value(j,k) * u_matrix.get_value(k,i);
                }

                if l_matrix.get_value(j,j) == 0.0 {
                    panic!("Matrix is singular, LU decomposition failed.");
                }
                u_matrix.set_value(j,i, (a_matrix.get_value(j,i) - sum) / l_matrix.get_value(j,j));
            }
        }
        println!("L in lin_algerbra = {:?}", l_matrix);
        println!("U in lin_algerbra = {:?}", u_matrix);

        let y_vector: Vector = LinAlgebra::foward_substitution(&l_matrix, &b_vector);
        let x_vector: Vector = LinAlgebra::back_substitution(&u_matrix, &y_vector);
        x_vector
    }

    fn foward_substitution(l_matrix: &Matrix, b_vector: &Vector) -> Vector {
        let n = l_matrix.n_cols();
        let mut y_vector: Vector = Vector::zeros(n);
        y_vector.set_value(0, b_vector.get_value(0) / l_matrix.get_value(0,0));

        for i in 1..(n) {
            let mut sum = 0.0;
            for j in 0..(n - 1) {
                sum += y_vector.get_value(j) * l_matrix.get_value(i,j);
            }
            y_vector.set_value(i, (b_vector.get_value(i) - sum) / l_matrix.get_value(i,i));
        }
        println!(" y_vector in lin_algebra = {:?}", y_vector);
        y_vector
    }

    fn back_substitution(u_matrix: &Matrix, y_vector: &Vector) -> Vector {
        let n = u_matrix.n_cols();
        let mut x_vector: Vector = Vector::zeros(n);

        x_vector.set_value(
            n - 1,
            y_vector.get_value(n - 1) / u_matrix.get_value(n - 1,n - 1),
        );
        //reverse loop "for" counted step by step
        for i in (0..=n - 1).rev().step_by(1) {
            let mut sum = 0.0;
            for j in i + 1..n {
                sum += x_vector.get_value(j) * u_matrix.get_value(i,j);
            }
            x_vector.set_value(i, (y_vector.get_value(i) - sum) / u_matrix.get_value(i,i));
        }
        println!(" x_vector in lin_algebra = {:?}", x_vector);
        x_vector
    }
}




//////////////////////////////
//Tests
/////////////////////////////


#[cfg(test)]
mod test{
    use crate::lin_algebra;
    use crate::matrix;
    use crate::vector;
    #[test]
    fn test_linear_eq_system_0(){
        let mut matrix_a = matrix::Matrix::zeros(3, 3);

        matrix_a.set_value(0, 0, 2.0);
        matrix_a.set_value(0, 1, -3.0);
        matrix_a.set_value(0, 2, 1.0);
        matrix_a.set_value(1, 0, 1.0);
        matrix_a.set_value(1, 1, 2.0);
        matrix_a.set_value(1, 2, -3.0);
        matrix_a.set_value(2, 0, 4.0);
        matrix_a.set_value(2, 1, -1.0);
        matrix_a.set_value(2, 2, -2.0);
    
        let mut vector_a = vector::Vector::zeros(3);
        vector_a.set_value(0, 1.0);
        vector_a.set_value(1, 4.0);
        vector_a.set_value(2, 8.0);
    
        let vector_b = lin_algebra::LinAlgebra::crout_method_solve(&matrix_a, &vector_a);

        assert_eq!(vec![3.0, 2.0, 1.0], vector_b.vector_values());
    }

    #[test]
    fn test_linear_eq_system_1(){
        let mut matrix_a = matrix::Matrix::zeros(3, 3);

        matrix_a.set_value(0, 0, 1.0);
        matrix_a.set_value(0, 1, 3.0);
        matrix_a.set_value(0, 2, -1.0);
        matrix_a.set_value(1, 0, 2.0);
        matrix_a.set_value(1, 1, 1.0);
        matrix_a.set_value(1, 2, 1.0);
        matrix_a.set_value(2, 0, 3.0);
        matrix_a.set_value(2, 1, -1.0);
        matrix_a.set_value(2, 2, 1.0);
    
        let mut vector_a = vector::Vector::zeros(3);
        vector_a.set_value(0, 0.0);
        vector_a.set_value(1, 1.0);
        vector_a.set_value(2, 3.0);
    
        let vector_b = lin_algebra::LinAlgebra::crout_method_solve(&matrix_a, &vector_a);

        assert_eq!(vec![1.0, -0.5, -0.5], vector_b.vector_values());
    }

    #[test]
    fn test_linear_eq_system_2(){
        let mut matrix_a = matrix::Matrix::zeros(3, 3);

        matrix_a.set_value(0, 0, 1.0);
        matrix_a.set_value(0, 1, 1.0);
        matrix_a.set_value(0, 2, 1.0);
        matrix_a.set_value(1, 0, 1.0);
        matrix_a.set_value(1, 1, 2.0);
        matrix_a.set_value(1, 2, 2.0);
        matrix_a.set_value(2, 0, 1.0);
        matrix_a.set_value(2, 1, 4.0);
        matrix_a.set_value(2, 2, 5.0);
    
        let mut vector_a = vector::Vector::zeros(3);
        vector_a.set_value(0, 1.0);
        vector_a.set_value(1, 2.0);
        vector_a.set_value(2, 4.0);
    
        let vector_b = lin_algebra::LinAlgebra::crout_method_solve(&matrix_a, &vector_a);

        assert_eq!(vec![0.0, 1.0, 0.0], vector_b.vector_values());
    }

    #[test]
    fn test_linear_eq_system_3(){
        let mut matrix_a = matrix::Matrix::zeros(2, 2);

        matrix_a.set_value(0, 0, 1.0);
        matrix_a.set_value(0, 1, -1.0);
        matrix_a.set_value(1, 0, 2.0);
        matrix_a.set_value(1, 1, 1.0);

    
        let mut vector_a = vector::Vector::zeros(3);
        vector_a.set_value(0, 4.0);
        vector_a.set_value(1, 8.0);
    
        let vector_b = lin_algebra::LinAlgebra::crout_method_solve(&matrix_a, &vector_a);

        assert_eq!(vec![4.0, 0.0], vector_b.vector_values());
    }
}


