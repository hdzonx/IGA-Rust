use atempts::atempt_lu_decomposition;

pub mod atempts;
pub mod lin_algebra;
pub mod matrix;
pub mod numerical_integration;
pub mod vector;

fn lu_decomposition_crout(
    matrix: &Vec<Vec<f64>>,
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), &'static str> {
    let n = matrix.len();
    let mut l_matrix = vec![vec![0.0; n]; n];
    let mut u_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        u_matrix[i][i] = 1.0;
    }

    for j in 0..n {
        for i in j..n {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l_matrix[i][k] * u_matrix[k][j];
            }
            l_matrix[i][j] = matrix[i][j] - sum;
        }

        for i in j..n {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l_matrix[j][k] * u_matrix[k][i];
            }

            if l_matrix[j][j] == 0.0 {
                return Err("Matrix is singular, LU decomposition failed.");
            }

            u_matrix[j][i] = (matrix[j][i] - sum) / l_matrix[j][j];
        }
    }
    Ok((l_matrix, u_matrix))
}

fn foward_substitution(l_matrix: &Vec<Vec<f64>>, b_vector: &Vec<f64>) -> Vec<f64> {
    let n = l_matrix.len();
    let mut y_vector: Vec<f64> = vec![0.0; n];

    y_vector[0] = b_vector[0] / l_matrix[0][0];

    for i in 1..(n) {
        let mut sum = 0.0;
        for j in 0..(n - 1) {
            sum += y_vector[j] * l_matrix[i][j];
        }
        y_vector[i] = (b_vector[i] - sum) / l_matrix[i][i];
    }

    println!(" y_vector in forward substitution");
    for row in &y_vector {
        println!("{:?}", row);
    }

    y_vector
}

fn back_substitution(u_matrix: &Vec<Vec<f64>>, c_vector: &Vec<f64>) -> Vec<f64> {
    let n = u_matrix.len();
    let mut x_vector: Vec<f64> = vec![0.0; n];

    x_vector[n - 1] = c_vector[n - 1] / u_matrix[n - 1][n - 1];
    //reverse loop "for" counted step by step
    for i in (0..=n - 1).rev().step_by(1) {
        let mut sum = 0.0;
        for j in i + 1..n {
            sum += x_vector[j] * u_matrix[i][j];
        }
        x_vector[i] = (c_vector[i] - sum) / u_matrix[i][i];
    }
    println!(" x_vector in back substitution");
    for row in &x_vector {
        println!("{:?}", row);
    }

    x_vector
}

fn main() {
    let mut matrix_a = matrix::Matrix::zeros(3, 3);
    //let mut matrix_b = matrix::Matrix::zeros(2,3);

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

    // println!("{:?}", matrix_c);
    //println!("{:?}", matrix_identity);
    //println!(" value is {}", value_in_pos);
    // println!("vector is {:?}", vector_b);

    //for atempst.rs
    #[cfg(test)]
    mod tests {
        #[test]
        fn exploration() {
            assert_eq!(2 + 2, 4);
        }

        #[test]
        fn another() {
            panic!("Make this test fail");
        }
    }

    //atempt_lu_decomposition(&mut matrix_a);

    //others implementations
    let matrix = vec![
        vec![2.0, -3.0, 1.0],
        vec![1.0, 2.0, -3.0],
        vec![4.0, -1.0, -2.0],
    ];

    if let Ok((l_matrix, u_matrix)) = lu_decomposition_crout(&matrix) {
        println!("L Matrix:");
        for row in &l_matrix {
            println!("{:?}", row);
        }

        println!("U Matrix:");
        for row in &u_matrix {
            println!("{:?}", row);
        }
    } else {
        println!("LU decomposition failed.");
    }

    let L_matrix = vec![
        vec![2.0, 0.0, 0.0],
        vec![1.0, 3.5, 0.0],
        vec![4.0, 5.0, 1.0],
    ];

    let b_vector = vec![1.0, 4.0, 8.0];

    foward_substitution(&L_matrix, &b_vector);

    let U_matrix = vec![
        vec![1.0, -1.5, 0.5],
        vec![0.0, 1.0, -1.0],
        vec![0.0, 0.0, 1.0],
    ];

    let y_vector = vec![0.5, 1.0, 1.0];
    back_substitution(&U_matrix, &y_vector);

    //let num = numerical_integration::NumericalIntegration::gauss_rule(6, 1);
    let mut num = numerical_integration::GaussRule::new(1, 1);
    num.gauss_rule();
    let weight = num.get_weights();
    println!("Weight in main: {:?}", weight);
    let abscissas = num.get_abscissas();
    println!("Abscissas in main: {:?}", abscissas);
}
