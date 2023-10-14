pub struct GaussRule {
    gauss_point_numbers: usize,
    dimension_problem: usize,
    weigth: Vec<f64>,
    abscissas: Vec<f64>,
}

impl GaussRule {
    pub fn new(gauss_point_numbers: usize, dimension_problem: usize) -> GaussRule {
        GaussRule {
            gauss_point_numbers,
            dimension_problem,
            weigth: vec![0.0; gauss_point_numbers],
            abscissas: vec![0.0; gauss_point_numbers],
        }
    }
    pub fn gauss_rule(&mut self) {
        let dimension_problem = self.dimension_problem;
        let gauss_point_numbers = self.gauss_point_numbers;

        if dimension_problem < 1 || dimension_problem > 3 {
            panic!(" not implemented for this dimesion.");
        }
        if gauss_point_numbers < 1 || gauss_point_numbers > 6 {
            panic!("not implmented for this number of Gauss integration points.")
        }

        let gauss = NumericalIntegration::values_gauss_rule(gauss_point_numbers);

        let vec_weight: Vec<_> = gauss.iter().map(|p| p.weigth.clone()).collect();
        let weight = vec_weight.get(0).unwrap().to_vec();
        self.weigth = weight;

        let vec_abscissas: Vec<_> = gauss.iter().map(|p| p.abscissas.clone()).collect();
        let abscissas = vec_abscissas.get(0).unwrap().to_vec();
        self.abscissas = abscissas;
    }

    pub fn get_weights(&self) -> Vec<f64> {
        self.weigth.clone()
    }

    pub fn get_abscissas(&self) -> Vec<f64> {
        self.abscissas.clone()
    }
}

#[derive(PartialEq, Debug)]
struct NumericalIntegration {
    gauss_points: usize,
    weigth: Vec<f64>,
    abscissas: Vec<f64>,
}
impl NumericalIntegration {
    fn values_gauss_rule(gauss_point_numbers: usize) -> Vec<NumericalIntegration> {
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
                weigth: vec![
                    0.2369268851,
                    0.4786286705,
                    0.56888888889,
                    0.4786286705,
                    0.2369268851,
                ],
                abscissas: vec![
                    -0.9061798459,
                    -0.5384693101,
                    0.0,
                    0.5384693101,
                    -0.9061798459,
                ],
            },
            NumericalIntegration {
                gauss_points: 6,
                weigth: vec![
                    0.1713244924,
                    0.3607615730,
                    0.4679139346,
                    0.4679139346,
                    0.3607615730,
                    0.1713244924,
                ],
                abscissas: vec![
                    -0.9324695142,
                    -0.6612093865,
                    -0.2386191861,
                    0.2386191861,
                    0.6612093865,
                    0.9324695142,
                ],
            },
        ];

        let my_gauss_values = Self::search_gauss_values(values, gauss_point_numbers);
        println!("gauss values = {:?}", my_gauss_values);

        my_gauss_values
    }

    fn search_gauss_values(
        values: Vec<NumericalIntegration>,
        gauss_point_numbers: usize,
    ) -> Vec<NumericalIntegration> {
        values
            .into_iter()
            .filter(|s| s.gauss_points == gauss_point_numbers)
            .collect()
    }
}

//////////////////////////////////////////
/////////////////////////////////////////
//Tests
///////////////////////////////////////////
///////////////////////////////////////////
#[cfg(test)]
mod tests {
    use crate::numerical_integration;
    #[test]
    fn test_new_value0() {
        let mut num = numerical_integration::GaussRule::new(1, 1);
        num.gauss_rule();
        assert_eq!(vec![2.0], num.get_weights());
        assert_eq!(vec![0.0], num.get_abscissas());
    }
    #[test]
    fn test_new_value1() {
        let mut num = numerical_integration::GaussRule::new(2, 1);
        num.gauss_rule();
        assert_eq!(vec![1.0, 1.0], num.get_weights());
        assert_eq!(vec![-0.57735027, 0.57735027], num.get_abscissas());
    }
    #[test]
    fn test_new_value3() {
        let mut num = numerical_integration::GaussRule::new(3, 1);
        num.gauss_rule();
        assert_eq!(
            vec![0.55555555555, 0.8888888888888, 0.55555555555],
            num.get_weights()
        );
        assert_eq!(vec![-0.774596669, 0.0, 0.774596669], num.get_abscissas());
    }
    #[test]
    fn test_new_value4() {
        let mut num = numerical_integration::GaussRule::new(4, 1);
        num.gauss_rule();
        assert_eq!(
            vec![0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451],
            num.get_weights()
        );
        assert_eq!(
            vec![-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116],
            num.get_abscissas()
        );
    }
    #[test]
    fn test_new_value5() {
        let mut num = numerical_integration::GaussRule::new(5, 1);
        num.gauss_rule();
        assert_eq!(
            vec![
                0.2369268851,
                0.4786286705,
                0.56888888889,
                0.4786286705,
                0.2369268851,
            ],
            num.get_weights()
        );
        assert_eq!(
            vec![
                -0.9061798459,
                -0.5384693101,
                0.0,
                0.5384693101,
                -0.9061798459,
            ],
            num.get_abscissas()
        );
    }
    #[test]
    fn test_new_value6() {
        let mut num = numerical_integration::GaussRule::new(6, 1);
        num.gauss_rule();
        assert_eq!(
            vec![
                0.1713244924,
                0.3607615730,
                0.4679139346,
                0.4679139346,
                0.3607615730,
                0.1713244924,
            ],
            num.get_weights()
        );
        assert_eq!(
            vec![
                -0.9324695142,
                -0.6612093865,
                -0.2386191861,
                0.2386191861,
                0.6612093865,
                0.9324695142,
            ],
            num.get_abscissas()
        );
    }
}
