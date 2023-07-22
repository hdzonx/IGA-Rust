use std::fmt::Debug;
#[derive(Debug)]
pub struct Vector {
    pub rows: usize,
    pub value: Vec<f64>,
}

impl Vector {
    pub fn set_value(&mut self, row: usize, data: f64) {
        self.value[row] = data;
    }

    pub fn get_value(&self, row: usize) -> f64 {
        self.value[row]
    }

    pub fn zeros(rows: usize) -> Vector {
        Vector {
            rows,
            value: vec![0.0; rows],
        }
    }

    pub fn add_vector(&self, other: &Vector) -> Vector {
        if self.rows != other.rows {
            panic!("add vector of different dimension is not possible")
        }
        let mut new_vector = Vector::zeros(self.rows);
        for i in 0..self.rows {
            new_vector.value[i] = self.value[i] + other.value[i];
        }
        new_vector
    }

    pub fn subtract_vector(&self, other: &Vector) -> Vector {
        if self.rows != other.rows {
            panic!("subtract vector of different dimension is not possible")
        }
        let mut new_vector = Vector::zeros(self.rows);
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


}
