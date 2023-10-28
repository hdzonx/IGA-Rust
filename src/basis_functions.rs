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

    fn basis_fn_num(&self) -> usize {
        let kno_vec_len = self.knot_vec.n_rows();
        let polin_order = self.polinomial_order;
        let n = kno_vec_len - polin_order - 1;
        n
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

    fn b_spline_vector(&self, bs_matrix: &Matrix, basis_fn_num: usize) -> Vector {
        let kno_vec_len = self.knot_vec.n_rows();
        let polin_order = self.polinomial_order;

        let mut bs_vector = Vector::zeros(basis_fn_num);

        for i in 0..basis_fn_num {
            bs_vector.set_value(i, bs_matrix.get_value(i, polin_order));
        }
        bs_vector
    }

    fn bspline_frst_deriv(&self, bs_matrix: &Matrix, basis_fn_num: usize) -> Vector {
        let mut vec_frst_deriv = Vector::zeros(basis_fn_num);

        for i in 0..basis_fn_num {
            let mut c1: f64 =
                self.knot_vec.get_value(i + self.polinomial_order) - self.knot_vec.get_value(i);
            if c1 != 0.0 {
                let old_c1 = c1;
                let val_pol = self.polinomial_order as f64;
                c1 = val_pol / old_c1;
            }

            let mut c2 = self.knot_vec.get_value(i + self.polinomial_order + 1)
                - self.knot_vec.get_value(i + 1);
            if c2 != 0.0 {
                let old_c2 = c2;
                let val_pol = self.polinomial_order as f64;
                c2 = val_pol / old_c2;
            }
            let val = c1 * bs_matrix.get_value(i, self.polinomial_order - 1)
                - c2 * bs_matrix.get_value(i + 1, self.polinomial_order - 1);
            vec_frst_deriv.set_value(i, val);
        }
        vec_frst_deriv
    }

    fn bspline_sec_derive(&self, bs_matrix: &Matrix, basis_fn_num: usize) -> Vector {
        let mut sec_deriv_vector = Vector::zeros(basis_fn_num);
        let p = self.polinomial_order;
        let val_pol = self.polinomial_order as f64;
        for i in 0..basis_fn_num {
            let mut c1 = (self.knot_vec.get_value(i + p) - self.knot_vec.get_value(i))
                * (self.knot_vec.get_value(i + p - 1) - self.knot_vec.get_value(i));

            if c1 != 0.0 {
                c1 = val_pol * (val_pol - 1.0) / c1;
            }

            let mut c2a = (self.knot_vec.get_value(i + p) - self.knot_vec.get_value(i + 1))
                * (self.knot_vec.get_value(i + p) - self.knot_vec.get_value(i));
            if c2a != 0.0 {
                c2a = val_pol * (val_pol - 1.0) / c2a;
            }

            let mut c2b = (self.knot_vec.get_value(i + p) - self.knot_vec.get_value(i + 1))
                * (self.knot_vec.get_value(i + p + 1) - self.knot_vec.get_value(i + 1));
            if c2b != 0.0 {
                c2b = val_pol * (val_pol - 1.0) / c2b;
            }

            let c2 = c2a + c2b;

            let mut c3 = (self.knot_vec.get_value(i + p + 1) - self.knot_vec.get_value(i + 1))
                * (self.knot_vec.get_value(i + p + 1) - self.knot_vec.get_value(i + 2));
            if c3 != 0.0 {
                c3 = val_pol * (val_pol - 1.0) / c3;
            }
            let mut result = c1 * bs_matrix.get_value(i, p - 2)
                - c2 * bs_matrix.get_value(i + 1, p - 2)
                + c3 * bs_matrix.get_value(i + 2, p - 2);

            sec_deriv_vector.set_value(i, result);
        }
        sec_deriv_vector
    }

    fn subreg_matrix(&self, sub_region_num: usize) -> Matrix {
        let mut subreg_matrix = Matrix::zeros(sub_region_num, 2);
        let mut subreg_vec = vec![0.0];

        for i in 1..self.knot_vec.n_rows() {
            if self.knot_vec.get_value(i) > self.knot_vec.get_value(i - 1) {
                let val = self.knot_vec.get_value(i);
                subreg_vec.push(val);
            }
        }
        for i in 0..subreg_matrix.n_rows() {
            for j in 1..subreg_matrix.n_cols() {
                subreg_matrix.set_value(i, j - 1, subreg_vec[i]);
                subreg_matrix.set_value(i, j, subreg_vec[i + 1]);
            }
        }
        subreg_matrix
    }
}

/////////////////////////
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

        let bs_vec_calc: Vector =
            BSpline::new(polinomial_order, knot_vector).b_spline_vector(&bs_mat_correct, 4);

        assert_eq!(bs_vector_correct, bs_vec_calc);
    }

    #[test]
    fn subreg_test_0() {
        let polinomial_order = 2;
        let subregion_num = 2;
        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut subreg_matrix_correct = Matrix::zeros(subregion_num, 2);
        subreg_matrix_correct.set_value(0, 0, 0.0);
        subreg_matrix_correct.set_value(0, 1, 0.5);
        subreg_matrix_correct.set_value(1, 0, 0.5);
        subreg_matrix_correct.set_value(1, 1, 1.0);

        let subreg_cal = BSpline::new(polinomial_order, knot_vector).subreg_matrix(subregion_num);
        assert_eq!(subreg_matrix_correct, subreg_cal);
    }

    #[test]
    fn b_spline_frst_derive_test_0() {
        let mut bs_mat = Matrix::zeros(6, 3);
        bs_mat.set_value(0, 2, 0.8659570925890132);
        bs_mat.set_value(1, 1, 0.9305681558000001);
        bs_mat.set_value(1, 2, 0.13163251691648037);
        bs_mat.set_value(2, 0, 1.0);
        bs_mat.set_value(2, 1, 0.06943184419999998);
        bs_mat.set_value(2, 2, 0.0024103904945065356);

        let basis_fn_num: usize = 4;
        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let derive_bs_matrix_calc = BSpline::new(2, knot_vector);
        let calc_val = derive_bs_matrix_calc.bspline_frst_deriv(&bs_mat, basis_fn_num);

        let mut vec_correct_bspline_deriv = Vector::zeros(4);
        vec_correct_bspline_deriv.set_value(0, -3.7222726232000003);
        vec_correct_bspline_deriv.set_value(1, 3.5834089348000004);
        vec_correct_bspline_deriv.set_value(2, 0.13886368839999996);
        vec_correct_bspline_deriv.set_value(3, 0.0);

        assert_eq!(vec_correct_bspline_deriv, calc_val);
    }

    #[test]
    fn b_spline_frst_derive_test_1() {
        let mut bs_mat = Matrix::zeros(6, 3);
        bs_mat.set_value(0, 2, 0.44888729930183624);
        bs_mat.set_value(1, 1, 0.6699905218);
        bs_mat.set_value(1, 2, 0.4966595728472456);
        bs_mat.set_value(2, 0, 1.0);
        bs_mat.set_value(2, 1, 0.33000947820000004);
        bs_mat.set_value(2, 2, 0.05445312785091815);

        let basis_fn_num: usize = 4;
        let mut knot_vector = Vector::zeros(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let derive_bs_matrix_calc = BSpline::new(2, knot_vector);
        let calc_val = derive_bs_matrix_calc.bspline_frst_deriv(&bs_mat, basis_fn_num);

        let mut vec_correct_bspline_deriv = Vector::zeros(4);
        vec_correct_bspline_deriv.set_value(0, -2.6799620872);
        vec_correct_bspline_deriv.set_value(1, 2.0199431307999998);
        vec_correct_bspline_deriv.set_value(2, 0.6600189564000001);
        vec_correct_bspline_deriv.set_value(3, 0.0);

        assert_eq!(vec_correct_bspline_deriv, calc_val);
    }
}
