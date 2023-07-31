use crate::matrix::Matrix;
use crate::vector::Vector;

//This file is only for tests and atempts for implement algorithms

pub fn atempt_lu_decomposition(A: &mut Matrix) {
    let n: usize = A.rows;
    let mut L: Matrix = Matrix::identity(n);
    let mut U: Matrix = Matrix::zeros(n, n);

    for k in 0..n {
        let akk = A.get_value(k, k);

        U.set_value(k, k, akk);

        for i in k + 1..n {
            let aik = A.get_value(i, k);
            let lik = aik / akk;
            L.set_value(i, k, lik);
            let uki = A.get_value(k, i);
            U.set_value(k, i, uki);
        }
        for i in k + 1..n {
            for j in k + 1..n {
                let aij = A.get_value(i, j);
                A.set_value(i, j, aij - L.get_value(i, k) * U.get_value(k, j));
            }
        }
    }
    println!("L = {:?}", L);
    println!("U = {:?}", U);
    println!("A = {:?}", A);
}
