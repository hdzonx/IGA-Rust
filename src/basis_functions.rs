use crate::matrix::Matrix;
use crate::vector::Vector;

pub struct BSpline {
    polinomial_order: usize,
    knot_vec: Vector,
}

impl BSpline {
    pub fn new(polinomial_order: usize, knot_vec: Vector) -> BSpline {
        BSpline {
            polinomial_order: polinomial_order,
            knot_vec: knot_vec,
        }
    }

    fn b_spline_matrix(&self, mut displacement: f64) -> Matrix {
        let rows: usize = self.knot_vec.n_rows() - 1;
        let cols: usize = self.polinomial_order + 1;
        let mut bs_matrix = Matrix::zeros(rows, cols);

        let base: f64 = 10.0;
        if displacement == 1.0 {
            displacement = 1.0 - base.powf(10.0);
        }

        for i in 0..rows {
            if displacement >= self.knot_vec.get_value(i)
                && displacement < self.knot_vec.get_value(i + 1)
            {
                bs_matrix.set_value(i, 0, 1.0);
            }
        }

        for k in 1..cols {
            for i in 0..rows - k {
                //let pod = k;
                let mut c1: f64 = self.knot_vec.get_value(i + k) - self.knot_vec.get_value(i);

                if c1 != 0.0 {
                    let old_c1: f64 = c1;
                    c1 = (displacement - self.knot_vec.get_value(i)) / old_c1;
                }

                let mut c2: f64 =
                    self.knot_vec.get_value(i + k + 1) - self.knot_vec.get_value(i + 1);

                if c2 != 0.0 {
                    let old_c2 = c2;
                    c2 = (self.knot_vec.get_value(i + k + 1) - displacement) / old_c2;
                }

                let value =
                    c1 * bs_matrix.get_value(i, k - 1) + c2 * bs_matrix.get_value(i + 1, k - 1);

                bs_matrix.set_value(i, k, value);
            }
        }
        bs_matrix
    }

    fn b_spline_vector(&self, bs_matrix: &Matrix) -> Vector {
        let kno_vec_len = self.knot_vec.n_rows();
        let polin_order = self.polinomial_order;
        let n = kno_vec_len - polin_order - 1;

        let mut bs_vector = Vector::zeros(n);

        for i in 0..n {
            bs_vector.set_value(i, bs_matrix.get_value(i, polin_order));
        }
        bs_vector
    }
}

//////////////////////////
//
//Tests
//
/////////////////////////

#[cfg(test)]

mod tests {
    use crate::basis_functions::BSpline;
    use crate::matrix::Matrix;
    use crate::vector::Vector;
    #[test]
    fn test_0_bs_matrix() {
        let polinomial_order = 2;
        let displacement = 0.03471592209999999;

        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::zeros(6, 3);
        bs_mat_correct.set_value(0, 2, 0.8659570925890132);
        bs_mat_correct.set_value(1, 1, 0.9305681558000001);
        bs_mat_correct.set_value(1, 2, 0.13163251691648037);
        bs_mat_correct.set_value(2, 0, 1.0);
        bs_mat_correct.set_value(2, 1, 0.06943184419999998);
        bs_mat_correct.set_value(2, 2, 0.0024103904945065356);

        let bs_matrix_calc = BSpline::new(polinomial_order, knot_vector);
        let calc = bs_matrix_calc.b_spline_matrix(displacement);

        assert_eq!(bs_mat_correct, calc);
    }
    #[test]
    fn test_1_bs_matrix() {
        let polinomial_order = 2;
        let displacement = 0.16500473910000002;

        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::zeros(6, 3);
        bs_mat_correct.set_value(0, 2, 0.44888729930183624);
        bs_mat_correct.set_value(1, 1, 0.6699905218);
        bs_mat_correct.set_value(1, 2, 0.4966595728472456);
        bs_mat_correct.set_value(2, 0, 1.0);
        bs_mat_correct.set_value(2, 1, 0.33000947820000004);
        bs_mat_correct.set_value(2, 2, 0.05445312785091815);

        let bs_matrix_calc = BSpline::new(polinomial_order, knot_vector);
        let calc = bs_matrix_calc.b_spline_matrix(displacement);
        assert_eq!(bs_mat_correct, calc);
    }

    #[test]
    fn test_2_bs_matrix() {
        let polinomial_order = 2;
        let displacement = 0.8349952609;

        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::zeros(6, 3);
        bs_mat_correct.set_value(1, 2, 0.05445312785091815);
        bs_mat_correct.set_value(2, 1, 0.33000947820000004);
        bs_mat_correct.set_value(2, 2, 0.4966595728472456);
        bs_mat_correct.set_value(3, 0, 1.0);
        bs_mat_correct.set_value(3, 1, 0.6699905218);
        bs_mat_correct.set_value(3, 2, 0.44888729930183624);

        let bs_matrix_calc = BSpline::new(polinomial_order, knot_vector);
        let calc = bs_matrix_calc.b_spline_matrix(displacement);
        assert_eq!(bs_mat_correct, calc);
    }

    #[test]
    fn test_0_bs_vector() {
        let polinomial_order = 2;
        let displacement = 0.03471592209999999;

        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::zeros(6, 3);
        bs_mat_correct.set_value(0, 2, 0.8659570925890132);
        bs_mat_correct.set_value(1, 1, 0.9305681558000001);
        bs_mat_correct.set_value(1, 2, 0.13163251691648037);
        bs_mat_correct.set_value(2, 0, 1.0);
        bs_mat_correct.set_value(2, 1, 0.06943184419999998);
        bs_mat_correct.set_value(2, 2, 0.0024103904945065356);

        let mut bs_vector_correct = Vector::zeros(4);
        bs_vector_correct.set_value(0, 0.8659570925890132);
        bs_vector_correct.set_value(1, 0.13163251691648037);
        bs_vector_correct.set_value(2, 0.0024103904945065356);

        let bs_vec_calc: Vector = BSpline::new(polinomial_order, knot_vector).b_spline_vector(&bs_mat_correct);

        assert_eq!(bs_vector_correct,bs_vec_calc);
    }
}
