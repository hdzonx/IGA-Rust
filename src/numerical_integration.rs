#[derive(PartialEq, Debug)]
pub struct NumericalIntegration {
    gauss_points: u32,
    weigth: Vec<f64>,
    abscissas: Vec<f64>,
}

impl NumericalIntegration {
    pub fn gauss_rule(gauss_point_numbers: u32, dimension_problem: i32) {
        if dimension_problem < 1 || dimension_problem > 3 {
            panic!(" not implemented for this dimesion.");
        }
        if gauss_point_numbers < 1 || gauss_point_numbers > 6 {
            panic!("not implmented for this number of Gauss integration points.")
        }
        let mut total_num_of_integr_points: u32 = 1;

        for i in 0..dimension_problem {
            total_num_of_integr_points *= gauss_point_numbers;
        }

        values_gauss_rule(gauss_point_numbers);
    }
}

fn values_gauss_rule(gauss_point_numbers: u32) {
    //  use crate::numerical_integration::Weights::one_point;
    //gauss_weight(one_point);
    let values = vec![
        NumericalIntegration {
            gauss_points: 1,
            weigth: vec![2.0],
            abscissas: vec![0.0],
        },
        NumericalIntegration {
            gauss_points: 2,
            weigth: vec![1.0, 1.0],
            abscissas: vec![-0.57735027, 0.57735027],
        },
        NumericalIntegration {
            gauss_points: 3,
            weigth: vec![0.55555555555, 0.8888888888888, 0.55555555555],
            abscissas: vec![-0.774596669, 0.0, 0.774596669],
        },
        NumericalIntegration {
            gauss_points: 4,
            weigth: vec![0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451],
            abscissas: vec![-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116],
        },
        NumericalIntegration {
            gauss_points: 5,
            weigth: vec![0.2369268851, 0.4786286705, 0.56888888889, 0.4786286705, 0.2369268851],
            abscissas: vec![-0.9061798459, -0.5384693101, 0.0, 0.5384693101, -0.9061798459],
        },
        NumericalIntegration {
            gauss_points: 6,
            weigth: vec![0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924],
            abscissas: vec![-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142],
        },
    ];

    let search_my_gauss_values = search_gauss_values(values, gauss_point_numbers);

    println!("gauss values = {:?}", search_my_gauss_values);
}

fn search_gauss_values(
    values: Vec<NumericalIntegration>,
    gauss_point_numbers: u32,
) -> Vec<NumericalIntegration> {
    values
        .into_iter()
        .filter(|s| s.gauss_points == gauss_point_numbers)
        .collect()
}
