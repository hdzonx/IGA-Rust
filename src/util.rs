struct Util {}

impl Util {
    pub fn newton_binomial(num: u128, den: u128) -> u128 {
        if den > num {
            panic!(" num must be greater or equal den")
        }
        let k = Util::factorial(num);
        let j = Util::factorial(den);
        let z = num - den;
        let binomial = k / (j * (Util::factorial(z)));
        binomial
    }

    fn factorial(num: u128) -> u128 {
        let mut val = 1;
        if num != 0 {
            val = num;
            for i in 1..num {
                val *= num - i;
            }
        }

        val
    }
}

#[cfg(test)]
mod tests {
    use crate::util;

    #[test]
    fn factorial_0() {
        let fac = util::Util::factorial(4);
        assert_eq!(24, fac);
    }

    #[test]
    fn factorial_1() {
        let fac = util::Util::factorial(0);
        assert_eq!(1, fac);
    }

    #[test]
    fn factorial_2() {
        let fac = util::Util::factorial(8);
        assert_eq!(40320, fac);
    }

    #[test]
    fn factorial_3() {
        let fac = util::Util::factorial(3);
        assert_eq!(6, fac);
    }

    #[test]
    fn factorial_4() {
        let fac = util::Util::factorial(5);
        assert_eq!(120, fac);
    }
    #[test]
    fn factorial_5() {
        let fac = util::Util::factorial(10);
        assert_eq!(3628800, fac);
    }
    #[test]
    fn binomial_0() {
        let bin = util::Util::newton_binomial(3, 3);
        assert_eq!(1, bin);
    }
    #[test]
    fn binomial_1() {
        let bin = util::Util::newton_binomial(1, 1);
        assert_eq!(1, bin);
    }
    #[test]
    fn binomial_2() {
        let bin = util::Util::newton_binomial(5, 3);
        assert_eq!(10, bin);
    }
    #[test]
    fn binomial_3() {
        let bin = util::Util::newton_binomial(5, 2);
        assert_eq!(10, bin);
    }
    #[test]
    fn binomial_4() {
        let bin = util::Util::newton_binomial(6, 5);
        assert_eq!(6, bin);
    }
    #[test]
    fn binomial_5() {
        let bin = util::Util::newton_binomial(0, 0);
        assert_eq!(1, bin);
    }
}
