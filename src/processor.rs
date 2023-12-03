use crate::basis_functions::BSpline;
use crate::matrix::Matrix;
use crate::numerical_integration::GaussRule;
use crate::vector::Vector;

struct Processor {}

impl Processor {
    fn nurbs_force_vector(
        subregion_Matrix: &Matrix,// mudar subregion_matrix para subregion number e refatorar 
        nurbs_weight: &Vec<f64>,
        gauss: GaussRule,
        polinomial_order: usize,
        knot_vector: Vector,
    ) {
        let subregion_num = subregion_Matrix.n_rows(); //subregion matrix has n rows and 2 columns
        let gauss_abscissas = gauss.abscissas();
        let gauss_weight = gauss.weights();
        let gauss_point_numbers = gauss.abscissas().len();
        let control_points_num = nurbs_weight.len();

        let bspline_function = BSpline::new(polinomial_order, knot_vector);

        let mut nurbs_vector = Vector::zeros(control_points_num);

        for i in 0..subregion_num {
            let subreg_initial = subregion_Matrix.get_value(i, 0);
            let subreg_final = subregion_Matrix.get_value(i, 1);
            for j in 0..gauss_point_numbers {
                let displacement = (1.0 - gauss_abscissas[j]) * subreg_initial / 2.0
                    + (1.0 + gauss_abscissas[j]) * subreg_final / 2.0;

                let bspline_matrix = bspline_function.b_spline_matrix(displacement);
                let bspline_vector =
                    bspline_function.b_spline_vector(&bspline_matrix, control_points_num);

                for m in 0..control_points_num {
                    let mut nurbs_den = 0.0;

                    for n in 0..control_points_num {
                        nurbs_den = bspline_vector.get_value(n) * nurbs_weight[n];
                    }
                    let nurbs_num = bspline_vector.get_value(m);
                    nurbs_vector.set_value(m, nurbs_num / nurbs_den);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::basis_functions;
    use crate::numerical_integration;
    use crate::processor;

    #[test]
    fn test_force_vector_0() {
        let mut gauss = numerical_integration::GaussRule::new(1, 1);
        gauss.gauss_rule();
        let nurbs_weight = [1.0, 1.0, 1.0, 1.0];

    }
}
