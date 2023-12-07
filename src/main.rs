pub mod matrix;
pub mod vector;
pub mod lin_algebra;
pub mod numerical_integration;
pub mod basis_functions;
pub mod processor;

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
    
    let mut num = numerical_integration::GaussRule::new(1, 1);
    num.gauss_rule();
    let weight = num.weights();
    println!("Weight in main: {:?}", weight);
    let abscissas = num.abscissas();
    println!("Abscissas in main: {:?}", abscissas);


}
