pub mod matrix;
pub mod vector;


fn main() {
    let mut matrix_a = matrix::Matrix::zeros(2,2);
    let mut matrix_b = matrix::Matrix::zeros(2,2);

    matrix_a.set_value(0, 1, 6.5);
    matrix_a.set_value(1, 1, 2.8);
    matrix_b.set_value(0, 0, 1.4);

    let matrix_c = matrix_a.add_matrices(&matrix_b);

    let matrx_d = matrix_c.transpose();

    let value_in_pos = matrix_c.get_value(0, 1);

    let matrix_identity = matrix::Matrix::identity(6);

    let vector_a = vector::Vector::zeros(3);

    println!("{:?}", matrix_c);
    println!("{:?}", matrix_identity);
    println!(" value is {}", value_in_pos);
    println!("vector is {:?}", vector_a);
}
