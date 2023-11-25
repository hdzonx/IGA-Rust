use crate::matrix::Matrix;
use crate::numerical_integration::GaussRule;
use crate::vector::Vector;

struct Processor {}

impl Processor {
    fn nurbs_force_vector(
        subregion_Matrix: &Matrix,
        nurbs_weight: &Vec<f64>,
        gauss: GaussRule,
    ) {
        let subregion_num = subregion_Matrix.n_rows(); //subregion matrix has n rows and 2 columns
        let gauss_abscissas = gauss.abscissas();
        let gauss_weight = gauss.weights();
        let gauss_point_numbers = gauss.abscissas().len();
        let control_points_num = nurbs_weight.len();

        for i in 0..subregion_num {
            let subreg_initial = subregion_Matrix.get_value(i, 0);
            let subreg_final = subregion_Matrix.get_value(i, 1);
            for j in 0..gauss_point_numbers {
                let displacement = (1.0 - gauss_abscissas[j]) * subreg_initial / 2.0
                    + (1.0 + gauss_abscissas[j]) * subreg_final / 2.0;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::numerical_integration;
    use crate::processor;

    #[test]
    fn test_force_vector_0() {
        let mut gauss = numerical_integration::GaussRule::new(1, 1);
        gauss.gauss_rule();
        let nurbs_weight = [1.0, 1.0, 1.0, 1.0];
        
    }
}
