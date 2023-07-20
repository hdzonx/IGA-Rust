pub mod matrix;
pub mod vector;
pub mod lin_algebra;


fn main() {
    let mut matrix_a = matrix::Matrix::zeros(2,3);
    let mut matrix_b = matrix::Matrix::zeros(2,3);

    matrix_a.set_value(0, 0, 2.0);
    matrix_a.set_value(0, 1, 1.0);
    matrix_a.set_value(0, 2, 3.0);
    matrix_b.set_value(1, 0, 4.0);
    matrix_b.set_value(1, 1, -3.0);
    matrix_b.set_value(1, 2, 5.0);

    let matrix_c = matrix_a.add_matrices(&matrix_b);

    //let matrx_d = matrix_c.transpose();

    //let value_in_pos = matrix_c.get_value(0, 1);

    //let matrix_identity = matrix::Matrix::identity(6);

    let mut vector_a = vector::Vector::zeros(3);
    vector_a.set_value(0, 2.0);
    vector_a.set_value(1, 3.0);
    vector_a.set_value(2, 5.0);

    let vector_b = lin_algebra::LinAlgebra::mult_matrix_by_vector(&matrix_c , &vector_a);

    println!("{:?}", matrix_c);
    //println!("{:?}", matrix_identity);
    //println!(" value is {}", value_in_pos);
    println!("vector is {:?}", vector_b);
}
