use std::thread::panicking;

use crate::matrix::Matrix;
use crate::util::{self, Util};
use crate::vector::Vector;
pub struct BasisFunctions {
    polinomial_order: usize,
    knot_vec: Vector,
}

impl BasisFunctions {
    pub fn new(polinomial_order: usize, knot_vec: Vector) -> BasisFunctions {
        BasisFunctions {
            polinomial_order: polinomial_order,
            knot_vec: knot_vec,
        }
    }

    fn basis_func_numbers(&self) -> usize {
        let kno_vec_len = self.knot_vec.n_rows();
        let polin_order = self.polinomial_order;
        let n = kno_vec_len - polin_order - 1;
        n
    }
    //Calculate the number of sub-Regions based on knot vector
    fn calc_sub_region_num(&self) -> u32 {
        let mut count: u32 = 0;
        for i in 1..self.knot_vec.n_rows() {
            if self.knot_vec.get_value(i - 1) != self.knot_vec.get_value(i) {
                count += 1;
            }
        }
        count
    }

    //Basis function of B-Spline in matricial format
    pub fn b_spline_matrix(&self, mut displacement: f64) -> Matrix {
        let rows: usize = self.knot_vec.n_rows() - 1;
        let cols: usize = self.polinomial_order + 1;
        let mut bs_matrix = Matrix::new(rows, cols);

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

    pub fn b_spline_vector(&self, bs_matrix: &Matrix, basis_fn_num: usize) -> Vector {
        let polin_order = self.polinomial_order;

        let mut bs_vector = Vector::new(basis_fn_num);

        for i in 0..basis_fn_num {
            bs_vector.set_value(i, bs_matrix.get_value(i, polin_order));
        }
        bs_vector
    }

    fn bspline_frst_deriv(&self, bs_matrix: &Matrix, basis_fn_num: usize) -> Vector {
        let mut vec_frst_deriv = Vector::new(basis_fn_num);

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

    fn bspline_secnd_deriv(&self, bs_matrix: &Matrix, basis_fn_num: usize) -> Vector {
        let pol_order = self.polinomial_order as f64;
        let mut vec_secnd_deriv = Vector::new(basis_fn_num);
        for i in 0..basis_fn_num {
            println!("iteração {i}");
            let mut c1 = (self.knot_vec.get_value(i + self.polinomial_order)
                - self.knot_vec.get_value(i))
                * (self.knot_vec.get_value(i + self.polinomial_order - 1)
                    - self.knot_vec.get_value(i));
            if c1 != 0.0 {
                c1 = pol_order * (pol_order - 1.0) / c1;
            }
            println!("c1 = {c1}");

            let mut c2a = (self.knot_vec.get_value(i + self.polinomial_order)
                - self.knot_vec.get_value(i + 1))
                * (self.knot_vec.get_value(i + self.polinomial_order) - self.knot_vec.get_value(i));
            if c2a != 0.0 {
                c2a = pol_order * (pol_order - 1.0) / c2a;
            }
            println!("c2a = {c2a}");
            let mut c2b = (self.knot_vec.get_value(i + self.polinomial_order)
                - self.knot_vec.get_value(i + 1))
                * (self.knot_vec.get_value(i + self.polinomial_order + 1)
                    - self.knot_vec.get_value(i + 1));
            if c2b != 0.0 {
                c2b = pol_order * (pol_order - 1.0) / c2b;
            }
            println!("c2b = {c2b}");
            let c2 = c2a + c2b;
            println!("c2 = {c2}");

            let mut c3 = (self.knot_vec.get_value(i + self.polinomial_order + 1)
                - self.knot_vec.get_value(i + 1))
                * (self.knot_vec.get_value(i + self.polinomial_order + 1)
                    - self.knot_vec.get_value(i + 2));
            if c3 != 0.0 {
                c3 = pol_order * (pol_order - 1.0) / c3;
            }
            println!("c3 = {c3}");

            let val = c1 * bs_matrix.get_value(i, self.polinomial_order - 2)
                - c2 * bs_matrix.get_value(i + 1, self.polinomial_order - 2)
                + c3 * bs_matrix.get_value(i + 2, self.polinomial_order - 2);

            vec_secnd_deriv.set_value(i, val);
            println!("val = {val}");
        }
        vec_secnd_deriv
    }

    pub fn subreg_matrix(&self, sub_region_num: usize) -> Matrix {
        let mut subreg_matrix = Matrix::new(sub_region_num, 2);
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

    pub fn nurbs_vector(&self, nurbs_weight: &Vec<f64>, bspline_vector: &Vector) -> Vector {
        let control_points_num = self.basis_func_numbers();
        let mut nurbs_vector = Vector::new(control_points_num);
        for m in 0..control_points_num {
            let mut nurbs_den = 0.0;

            for n in 0..control_points_num {
                if nurbs_weight[n] <= 0.0 {
                    panic!(" NURBS weight values must be greater or equal to 1")
                }
                nurbs_den += bspline_vector.get_value(n) * nurbs_weight[n];
            }
            let nurbs_num = bspline_vector.get_value(m) * nurbs_weight[m];
            nurbs_vector.set_value(m, nurbs_num / nurbs_den);
        }
        nurbs_vector
    }
    fn weight_dot_vector(&self, weight: &Vec<f64>, vector: &Vector) -> f64 {
        let mut val = 0.0;
        for i in 0..weight.len() {
            val += weight[i] * vector.get_value(i);
        }
        val
    }
    //first derivative for NURBS basis functions
    fn nurbs_frst_deriv(
        &self,
        nurbs_weight: &Vec<f64>,
        bspline: &Vector,
        bspline_frst_deriv: &Vector,
    ) -> Vector {
        let mut nurbs_ft_deriv = Vector::new(nurbs_weight.len());

        for i in 0..nurbs_weight.len() {
            let w = self.weight_dot_vector(nurbs_weight, bspline);
            let w_deriv = self.weight_dot_vector(nurbs_weight, bspline_frst_deriv);
            let val = nurbs_weight[i]
                * (w * bspline_frst_deriv.get_value(i) - w_deriv * bspline.get_value(i) / (w * w));
            nurbs_ft_deriv.set_value(i, val);
        }
        nurbs_ft_deriv
    }
    //second derivative for NURBS basis functions
    fn nurbs_secnd_deriv(
        &self,
        weight_vector: Vec<f64>,
        bspline: Vector,
        bspline_frst_deriv: Vector,
        bspline_sec_derive: Vector,
        nurbs_vector: Vector,
        nurbs_frst_deriv: Vector,
    ) -> Vector {
        let order = self.polinomial_order;
        let mut nurbs_sec_deriv_vec = Vector::new(weight_vector.len());
        let w = self.weight_dot_vector(&weight_vector, &bspline);
        //The value 2 bellow is due the second derivative
        for i in 0..2 {
            let a_ik = weight_vector[i] * bspline_sec_derive.get_value(i);

            let mut sum = 0.0;

            for j in 0..order {
                let w_deriv: f64;
                let nurbs_i: f64;
                //The value 2 bellow is due the second derivative
                let binomial = Util::newton_binomial(2, j as u64);
                if j == 0 {
                    w_deriv = self.weight_dot_vector(&weight_vector, &bspline_frst_deriv);
                    nurbs_i = nurbs_frst_deriv.get_value(i);
                } else if j == 1 {
                    w_deriv = self.weight_dot_vector(&weight_vector, &bspline_sec_derive);
                    nurbs_i = nurbs_vector.get_value(i);
                } else {
                    panic!(" Not implemented for derivative greater 2");
                }
                sum += (binomial as f64) * w_deriv * nurbs_i;
            }

            let val = (a_ik - sum) / w;
            nurbs_sec_deriv_vec.set_value(i, val);
        }
        nurbs_sec_deriv_vec
    }
}

/////////////////////////
//
//Tests
//
/////////////////////////

#[cfg(test)]

mod tests {
    use crate::basis_functions::BasisFunctions;
    use crate::matrix::Matrix;
    use crate::vector::Vector;

    #[test]
    fn first_deriv_nurbs_0() {
        let mut bs_mat = Matrix::new(6, 3);
        bs_mat.set_value(0, 2, 0.8659570925890132);
        bs_mat.set_value(1, 1, 0.9305681558000001);
        bs_mat.set_value(1, 2, 0.13163251691648037);
        bs_mat.set_value(2, 0, 1.0);
        bs_mat.set_value(2, 1, 0.06943184419999998);
        bs_mat.set_value(2, 2, 0.0024103904945065356);

        let basis_fn_num: usize = 4;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let basis_function = BasisFunctions::new(2, knot_vector);
        let bspline = &basis_function.b_spline_vector(&bs_mat, basis_fn_num);
        let bspline_derive = basis_function.bspline_frst_deriv(&bs_mat, basis_fn_num);

        let mut vec_correct_bspline_deriv = Vector::new(basis_fn_num);
        vec_correct_bspline_deriv.set_value(0, -3.722272623200001);
        vec_correct_bspline_deriv.set_value(1, 3.5834089348000013);
        vec_correct_bspline_deriv.set_value(2, 0.1388636884);
        vec_correct_bspline_deriv.set_value(3, 0.0);

        //If weight values equals to 1, bsplines derivatives equals nurbs derivatives
        let nurbs_weight: &Vec<f64> = &vec![1.0, 1.0, 1.0, 1.0];
        let nurbs_fisrt_derive =
            basis_function.nurbs_frst_deriv(nurbs_weight, bspline, &bspline_derive);

        assert_eq!(vec_correct_bspline_deriv, nurbs_fisrt_derive);
    }

    #[test]
    fn partition_unit_test_nurbs_0() {
        let nurbs_weight: &Vec<f64> = &vec![1.0, 1.0, 1.0, 1.0];

        let polinomial_order: usize = 2;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs = Vector::new(4);
        bs.set_value(0, 0.0);
        bs.set_value(1, 0.625);
        bs.set_value(2, 0.125);
        bs.set_value(3, 0.25);

        let basis_fn = BasisFunctions::new(polinomial_order, knot_vector);
        let nurbs = basis_fn.nurbs_vector(nurbs_weight, &bs);
        let mut val = 0.0;
        for i in 0..nurbs.n_rows() {
            val += nurbs.get_value(i);
        }
        assert_eq!(val, 1.0);
    }
    #[test]
    fn partition_unit_test_nurbs_1() {
        let nurbs_weight: &Vec<f64> = &vec![1.0, 2.0, 3.0, 4.0];

        let polinomial_order: usize = 2;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs = Vector::new(4);
        bs.set_value(0, 0.);
        bs.set_value(1, 0.625);
        bs.set_value(2, 0.125);
        bs.set_value(3, 0.25);

        let basis_fn = BasisFunctions::new(polinomial_order, knot_vector);
        let nurbs = basis_fn.nurbs_vector(nurbs_weight, &bs);
        let mut val = 0.0;
        for i in 0..nurbs.n_rows() {
            val += nurbs.get_value(i);
        }
        assert_eq!(val, 1.0);
    }

    #[test]
    fn partition_unit_test_nurbs_2() {
        let nurbs_weight: &Vec<f64> = &vec![1.0, 4.0, 4.0, 1.0];

        let polinomial_order: usize = 2;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs = Vector::new(4);
        bs.set_value(0, 0.);
        bs.set_value(1, 0.625);
        bs.set_value(2, 0.125);
        bs.set_value(3, 0.25);

        let basis_fn = BasisFunctions::new(polinomial_order, knot_vector);
        let nurbs = basis_fn.nurbs_vector(nurbs_weight, &bs);
        let mut val = 0.0;
        for i in 0..nurbs.n_rows() {
            val += nurbs.get_value(i);
        }
        assert_eq!(val, 1.0);
    }

    #[test]
    fn partition_unit_test_nurbs_3() {
        let nurbs_weight: &Vec<f64> = &vec![1.0, 1.0, 1.0, 1.0];

        let polinomial_order: usize = 2;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs = Vector::new(4);
        bs.set_value(0, 0.);
        bs.set_value(1, 0.625);
        bs.set_value(2, 0.125);
        bs.set_value(3, 0.25);

        let basis_fn = BasisFunctions::new(polinomial_order, knot_vector);
        let nurbs = basis_fn.nurbs_vector(nurbs_weight, &bs);
        let mut val = 0.0;
        for i in 0..nurbs.n_rows() {
            val += nurbs.get_value(i);
        }
        assert_eq!(val, 1.0);
    }
    #[test]
    fn calc_subregion_test0() {
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let polinomial_order = 2;

        let bs_matrix_calc = BasisFunctions::new(polinomial_order, knot_vector);
        let sub_region_num_calculated = bs_matrix_calc.calc_sub_region_num();
        assert_eq!(2, sub_region_num_calculated);
    }

    #[test]
    fn calc_subregion_test1() {
        let mut knot_vector = Vector::new(11);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.0);
        knot_vector.set_value(4, 0.25);
        knot_vector.set_value(5, 0.5);
        knot_vector.set_value(6, 0.75);
        knot_vector.set_value(7, 1.0);
        knot_vector.set_value(8, 1.0);
        knot_vector.set_value(9, 1.0);
        knot_vector.set_value(10, 1.0);

        let polinomial_order = 3;

        let bs_matrix_calc = BasisFunctions::new(polinomial_order, knot_vector);
        let sub_region_num_calculated = bs_matrix_calc.calc_sub_region_num();
        assert_eq!(4, sub_region_num_calculated);
    }

    #[test]
    fn calc_subregion_test3() {
        let mut knot_vector = Vector::new(15);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.0);
        knot_vector.set_value(4, 0.0);
        knot_vector.set_value(5, 0.125);
        knot_vector.set_value(6, 0.25);
        knot_vector.set_value(7, 0.50);
        knot_vector.set_value(8, 0.75);
        knot_vector.set_value(9, 0.875);
        knot_vector.set_value(10, 1.0);
        knot_vector.set_value(11, 1.0);
        knot_vector.set_value(12, 1.0);
        knot_vector.set_value(13, 1.0);
        knot_vector.set_value(14, 1.0);

        let polinomial_order = 4;

        let bs_matrix_calc = BasisFunctions::new(polinomial_order, knot_vector);
        let sub_region_num_calculated = bs_matrix_calc.calc_sub_region_num();
        assert_eq!(6, sub_region_num_calculated);
    }

    #[test]
    fn test_0_bs_matrix() {
        let polinomial_order = 2;
        let displacement = 0.03471592209999999;

        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::new(6, 3);
        bs_mat_correct.set_value(0, 2, 0.8659570925890132);
        bs_mat_correct.set_value(1, 1, 0.9305681558000001);
        bs_mat_correct.set_value(1, 2, 0.13163251691648037);
        bs_mat_correct.set_value(2, 0, 1.0);
        bs_mat_correct.set_value(2, 1, 0.06943184419999998);
        bs_mat_correct.set_value(2, 2, 0.0024103904945065356);

        let bs_matrix_calc = BasisFunctions::new(polinomial_order, knot_vector);
        let calc = bs_matrix_calc.b_spline_matrix(displacement);

        assert_eq!(bs_mat_correct, calc);
    }
    #[test]
    fn test_1_bs_matrix() {
        let polinomial_order = 2;
        let displacement = 0.16500473910000002;

        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::new(6, 3);
        bs_mat_correct.set_value(0, 2, 0.44888729930183624);
        bs_mat_correct.set_value(1, 1, 0.6699905218);
        bs_mat_correct.set_value(1, 2, 0.4966595728472456);
        bs_mat_correct.set_value(2, 0, 1.0);
        bs_mat_correct.set_value(2, 1, 0.33000947820000004);
        bs_mat_correct.set_value(2, 2, 0.05445312785091815);

        let bs_matrix_calc = BasisFunctions::new(polinomial_order, knot_vector);
        let calc = bs_matrix_calc.b_spline_matrix(displacement);
        assert_eq!(bs_mat_correct, calc);
    }

    #[test]
    fn test_2_bs_matrix() {
        let polinomial_order = 2;
        let displacement = 0.8349952609;

        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::new(6, 3);
        bs_mat_correct.set_value(1, 2, 0.05445312785091815);
        bs_mat_correct.set_value(2, 1, 0.33000947820000004);
        bs_mat_correct.set_value(2, 2, 0.4966595728472456);
        bs_mat_correct.set_value(3, 0, 1.0);
        bs_mat_correct.set_value(3, 1, 0.6699905218);
        bs_mat_correct.set_value(3, 2, 0.44888729930183624);

        let bs_matrix_calc = BasisFunctions::new(polinomial_order, knot_vector);
        let calc = bs_matrix_calc.b_spline_matrix(displacement);
        assert_eq!(bs_mat_correct, calc);
    }

    #[test]
    fn test_0_bs_vector() {
        let polinomial_order = 2;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut bs_mat_correct = Matrix::new(6, 3);
        bs_mat_correct.set_value(0, 2, 0.8659570925890132);
        bs_mat_correct.set_value(1, 1, 0.9305681558000001);
        bs_mat_correct.set_value(1, 2, 0.13163251691648037);
        bs_mat_correct.set_value(2, 0, 1.0);
        bs_mat_correct.set_value(2, 1, 0.06943184419999998);
        bs_mat_correct.set_value(2, 2, 0.0024103904945065356);

        let mut bs_vector_correct = Vector::new(4);
        bs_vector_correct.set_value(0, 0.8659570925890132);
        bs_vector_correct.set_value(1, 0.13163251691648037);
        bs_vector_correct.set_value(2, 0.0024103904945065356);

        let bs_vec_calc: Vector =
            BasisFunctions::new(polinomial_order, knot_vector).b_spline_vector(&bs_mat_correct, 4);

        assert_eq!(bs_vector_correct, bs_vec_calc);
    }

    #[test]
    fn subreg_test_0() {
        let polinomial_order = 2;
        let subregion_num = 2;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let mut subreg_matrix_correct = Matrix::new(subregion_num, 2);
        subreg_matrix_correct.set_value(0, 0, 0.0);
        subreg_matrix_correct.set_value(0, 1, 0.5);
        subreg_matrix_correct.set_value(1, 0, 0.5);
        subreg_matrix_correct.set_value(1, 1, 1.0);

        let subreg_cal =
            BasisFunctions::new(polinomial_order, knot_vector).subreg_matrix(subregion_num);
        assert_eq!(subreg_matrix_correct, subreg_cal);
    }

    #[test]
    fn b_spline_frst_derive_test_0() {
        let mut bs_mat = Matrix::new(6, 3);
        bs_mat.set_value(0, 2, 0.8659570925890132);
        bs_mat.set_value(1, 1, 0.9305681558000001);
        bs_mat.set_value(1, 2, 0.13163251691648037);
        bs_mat.set_value(2, 0, 1.0);
        bs_mat.set_value(2, 1, 0.06943184419999998);
        bs_mat.set_value(2, 2, 0.0024103904945065356);

        let basis_fn_num: usize = 4;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let derive_bs_matrix_calc = BasisFunctions::new(2, knot_vector);
        let calc_val = derive_bs_matrix_calc.bspline_frst_deriv(&bs_mat, basis_fn_num);

        let mut vec_correct_bspline_deriv = Vector::new(4);
        vec_correct_bspline_deriv.set_value(0, -3.7222726232000003);
        vec_correct_bspline_deriv.set_value(1, 3.5834089348000004);
        vec_correct_bspline_deriv.set_value(2, 0.13886368839999996);
        vec_correct_bspline_deriv.set_value(3, 0.0);

        assert_eq!(vec_correct_bspline_deriv, calc_val);
    }

    #[test]
    fn b_spline_frst_derive_test_1() {
        let mut bs_mat = Matrix::new(6, 3);
        bs_mat.set_value(0, 2, 0.44888729930183624);
        bs_mat.set_value(1, 1, 0.6699905218);
        bs_mat.set_value(1, 2, 0.4966595728472456);
        bs_mat.set_value(2, 0, 1.0);
        bs_mat.set_value(2, 1, 0.33000947820000004);
        bs_mat.set_value(2, 2, 0.05445312785091815);

        let basis_fn_num: usize = 4;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let derive_bs_matrix_calc = BasisFunctions::new(2, knot_vector);
        let calc_val = derive_bs_matrix_calc.bspline_frst_deriv(&bs_mat, basis_fn_num);

        let mut vec_correct_bspline_deriv = Vector::new(4);
        vec_correct_bspline_deriv.set_value(0, -2.6799620872);
        vec_correct_bspline_deriv.set_value(1, 2.0199431307999998);
        vec_correct_bspline_deriv.set_value(2, 0.6600189564000001);
        vec_correct_bspline_deriv.set_value(3, 0.0);

        assert_eq!(vec_correct_bspline_deriv, calc_val);
    }

    #[test]
    fn b_spline_secnd_deriv_test_0() {
        let mut bs_mat = Matrix::new(6, 3);
        bs_mat.set_value(0, 2, 0.8659570925890132);
        bs_mat.set_value(1, 1, 0.9305681558000001);
        bs_mat.set_value(1, 2, 0.13163251691648037);
        bs_mat.set_value(2, 0, 1.0);
        bs_mat.set_value(2, 1, 0.06943184419999998);
        bs_mat.set_value(2, 2, 0.0024103904945065356);

        let basis_fn_num: usize = 4;
        let mut knot_vector = Vector::new(7);
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let derive_bs_matrix_calc = BasisFunctions::new(2, knot_vector);
        let calc_val = derive_bs_matrix_calc.bspline_secnd_deriv(&bs_mat, basis_fn_num);

        let mut vec_correct_bspline_deriv = Vector::new(4);
        vec_correct_bspline_deriv.set_value(0, 8.0);
        vec_correct_bspline_deriv.set_value(1, -12.0);
        vec_correct_bspline_deriv.set_value(2, 4.0);
        vec_correct_bspline_deriv.set_value(3, 0.0);

        assert_eq!(vec_correct_bspline_deriv, calc_val);
    }
}
