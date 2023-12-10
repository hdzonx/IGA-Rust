use std::fmt::Debug;
#[derive(Debug, PartialEq)]
pub struct Vector {
    rows: usize,
    value: Vec<f64>,
}

impl Vector {
    pub fn n_rows(&self) -> usize {
        self.rows
    }

    pub fn vector_values(&self) -> Vec<f64> {
        self.value.clone()
    }

    pub fn set_value(&mut self, row: usize, val: f64) {
        self.value[row] = val;
    }

    pub fn get_value(&self, row: usize) -> f64 {
        self.value[row]
    }

    pub fn new(rows: usize) -> Vector {
        Vector {
            rows,
            value: vec![0.0; rows],
        }
    }

    pub fn add_vector(&self, other: &Vector) -> Vector {
        if self.rows != other.rows {
            panic!("add vector of different dimension is not possible")
        }
        let mut new_vector = Vector::new(self.rows);
        for i in 0..self.rows {
            new_vector.value[i] = self.value[i] + other.value[i];
        }
        new_vector
    }

    pub fn subtract_vector(&self, other: &Vector) -> Vector {
        if self.rows != other.rows {
            panic!("subtract vector of different dimension is not possible")
        }
        let mut new_vector = Vector::new(self.rows);
        for i in 0..self.rows {
            new_vector.value[i] = self.value[i] - other.value[i];
        }
        new_vector
    }

    pub fn dot(&self, other: &Vector) -> f64 {
        if self.rows != other.rows {
            panic!("dot vector of different dimension is not possible")
        }
        let mut dot_value: f64 = 0.0;
        for i in 0..self.rows {
            dot_value += self.value[i] + other.value[i];
        }
        dot_value
    }
    pub fn scalar_by_vector(&self, scalar: f64) -> Vector {
        let mut new_vector = Vector::new(self.rows);

        for i in 0..self.rows {
            let val = scalar * self.value[i];
            new_vector.set_value(i, val);
        }
        new_vector
    }
}

#[cfg(test)]
mod tests {
    use crate::vector;
    #[test]
    fn test_new_scalar_by_vector_test_0() {
        let mut old_vector = vector::Vector::new(4);
        old_vector.set_value(0, 1.);
        old_vector.set_value(1, 2.);
        old_vector.set_value(2, 3.);
        old_vector.set_value(3, 4.);

        let scalar = 3.0;

        let mut correct_vector = vector::Vector::new(4);
        correct_vector.set_value(0, 3.);
        correct_vector.set_value(1, 6.);
        correct_vector.set_value(2, 9.);
        correct_vector.set_value(3, 12.);

        let new_vector = old_vector.scalar_by_vector(scalar);
        assert_eq!(correct_vector, new_vector)
    }
    #[test]
    fn test_new_scalar_by_vector_test_1() {
        let mut old_vector = vector::Vector::new(4);
        old_vector.set_value(0, 0.);
        old_vector.set_value(1, 2.);
        old_vector.set_value(2, 3.);
        old_vector.set_value(3, -4.);

        let scalar = 3.0;

        let mut correct_vector = vector::Vector::new(4);
        correct_vector.set_value(0, 0.);
        correct_vector.set_value(1, 6.);
        correct_vector.set_value(2, 9.);
        correct_vector.set_value(3, -12.);

        let new_vector = old_vector.scalar_by_vector(scalar);
        assert_eq!(correct_vector, new_vector)
    }
}
