use crate::basis_functions::BasisFunctions;
use crate::matrix::Matrix;
use crate::numerical_integration::GaussRule;
use crate::vector::Vector;

struct Processor {}

impl Processor {
    fn force_vector(
        subregion_num: usize,
        nurbs_weight: &Vec<f64>,
        gauss: GaussRule,
        polinomial_order: usize,
        knot_vector: Vector,
        beam_lenght: f64,
        load: f64,
    ) ->Vector {
        let gauss_abscissas = gauss.abscissas();
        let gauss_weight = gauss.weights();
        let gauss_point_numbers = gauss.abscissas().len();
        let control_points_num = nurbs_weight.len();

        let basis_function = BasisFunctions::new(polinomial_order, knot_vector);

        let subregion_matrix = basis_function.subreg_matrix(subregion_num); //subregion matrix has n rows and 2 columns
        //let mut nurbs_vector = Vector::new(control_points_num);
        let mut force_vector = Vector::new(control_points_num);

        for i in 0..subregion_num {
            let subreg_initial = subregion_matrix.get_value(i, 0);
            let subreg_final = subregion_matrix.get_value(i, 1);

            for j in 0..gauss_point_numbers {
                let displacement = (1.0 - gauss_abscissas[j]) * subreg_initial / 2.0
                    + (1.0 + gauss_abscissas[j]) * subreg_final / 2.0;
                let bspline_matrix = basis_function.b_spline_matrix(displacement);
                let bspline_vector =
                basis_function.b_spline_vector(&bspline_matrix, control_points_num);
                let nurbs = basis_function.nurbs_vector(&nurbs_weight, &bspline_vector);
                println!("b_spline vector = {:?}", bspline_vector);
                println!("nurbs vector = {:?}", nurbs);

                let dxdxsi = (subreg_final - subreg_initial) * beam_lenght / 2.0; //verificar a compatibilidade desta função com nurbs
                let scalar = load * dxdxsi * gauss_weight[j];
                let d_force_vector = nurbs.scalar_by_vector(scalar);
                println!("dforce vector = {:?}", d_force_vector);

                force_vector.add_vector(d_force_vector);
            }
    
        }
        
        println!("force vector = {:?}", force_vector);
        force_vector
    }

}

#[cfg(test)]
mod tests {
    use crate::basis_functions;
    use crate::numerical_integration;
    use crate::processor;
    use crate::vector;

    #[test]
    fn test_force_vector_0() {
        let mut gauss = numerical_integration::GaussRule::new(1, 1);
        gauss.gauss_rule();
        let nurbs_weight: &Vec<f64> = &vec![1.0, 1.0, 1.0, 1.0];
        let subregion_num: usize = 2;
        let polinomial_order: usize = 2;
        let mut knot_vector = vector::Vector::new(7);
        let beam_lenght: f64 = 2.0;
        let load: f64 = 25.0;
        knot_vector.set_value(0, 0.0);
        knot_vector.set_value(1, 0.0);
        knot_vector.set_value(2, 0.0);
        knot_vector.set_value(3, 0.5);
        knot_vector.set_value(4, 1.0);
        knot_vector.set_value(5, 1.0);
        knot_vector.set_value(6, 1.0);

        let force = processor::Processor::force_vector(
            subregion_num,
            nurbs_weight,
            gauss,
            polinomial_order,
            knot_vector,
            beam_lenght,
            load,
        );
    }
}
