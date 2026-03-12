use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

use crate::error::{KnotError, Result};

/// Integer-coefficient polynomial in one variable (t).
/// `coeffs[i]` is the coefficient of t^i.
/// Trailing zeros are trimmed after each operation.
#[derive(Clone, Debug)]
pub struct Polynomial {
    pub coeffs: Vec<i64>,
}

impl Polynomial {
    pub fn zero() -> Self {
        Polynomial { coeffs: vec![] }
    }

    pub fn one() -> Self {
        Polynomial { coeffs: vec![1] }
    }

    pub fn t() -> Self {
        Polynomial { coeffs: vec![0, 1] }
    }

    pub fn one_minus_t() -> Self {
        Polynomial {
            coeffs: vec![1, -1],
        }
    }

    pub fn from_constant(c: i64) -> Self {
        if c == 0 {
            Self::zero()
        } else {
            Polynomial { coeffs: vec![c] }
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Lowest degree with nonzero coefficient.
    pub fn ldegree(&self) -> usize {
        for (i, &c) in self.coeffs.iter().enumerate() {
            if c != 0 {
                return i;
            }
        }
        0
    }

    fn trim(&mut self) {
        while self.coeffs.last() == Some(&0) {
            self.coeffs.pop();
        }
    }

    /// Normalize: divide by t^ldegree (shift coefficients so constant term is nonzero),
    /// then ensure leading coefficient is positive.
    pub fn normalize(&self) -> Polynomial {
        if self.is_zero() {
            return Self::zero();
        }
        let ld = self.ldegree();
        let shifted: Vec<i64> = self.coeffs[ld..].to_vec();
        let mut p = Polynomial { coeffs: shifted };
        p.trim();
        if p.coeffs.last().map_or(false, |&c| c < 0) {
            p = -p;
        }
        p
    }

    /// Negate all coefficients.
    pub fn negate(&self) -> Polynomial {
        Polynomial {
            coeffs: self.coeffs.iter().map(|&c| -c).collect(),
        }
    }

    /// Exact polynomial long division: a / b where b divides a exactly.
    /// Returns quotient. Panics if not exact (remainder != 0).
    pub fn exact_div(a: &Polynomial, b: &Polynomial) -> Polynomial {
        if b.is_zero() {
            panic!("polynomial division by zero");
        }
        if a.is_zero() {
            return Polynomial::zero();
        }

        let mut rem = a.coeffs.clone();
        let b_deg = b.coeffs.len() - 1;
        let b_lead = *b.coeffs.last().unwrap();

        if rem.len() < b.coeffs.len() {
            // a has lower degree than b — should be zero if division is exact
            assert!(rem.iter().all(|&c| c == 0), "non-exact polynomial division");
            return Polynomial::zero();
        }

        let q_len = rem.len() - b_deg;
        let mut q = vec![0i64; q_len];

        for i in (0..q_len).rev() {
            let idx = i + b_deg;
            assert!(
                rem[idx] % b_lead == 0,
                "non-exact polynomial division at degree {idx}: {rem:?} / {b_lead}"
            );
            let coeff = rem[idx] / b_lead;
            q[i] = coeff;
            for j in 0..=b_deg {
                rem[i + j] -= coeff * b.coeffs[j];
            }
        }

        debug_assert!(
            rem.iter().all(|&c| c == 0),
            "non-exact polynomial division: remainder {rem:?}"
        );

        let mut result = Polynomial { coeffs: q };
        result.trim();
        result
    }

    /// Bareiss fraction-free determinant of a matrix of polynomials.
    /// The matrix is n×n, stored as `matrix[row][col]`.
    pub fn determinant(matrix: &[Vec<Polynomial>]) -> Polynomial {
        let n = matrix.len();
        if n == 0 {
            return Polynomial::one();
        }
        if n == 1 {
            return matrix[0][0].clone();
        }

        let mut m: Vec<Vec<Polynomial>> = matrix.to_vec();
        let mut sign = 1i64;
        let mut prev_pivot = Polynomial::one();

        for k in 0..n {
            // Find a nonzero pivot in column k, rows k..n
            let mut pivot_row = None;
            for i in k..n {
                if !m[i][k].is_zero() {
                    pivot_row = Some(i);
                    break;
                }
            }
            let pivot_row = match pivot_row {
                Some(r) => r,
                None => return Polynomial::zero(),
            };

            if pivot_row != k {
                m.swap(k, pivot_row);
                sign = -sign;
            }

            let pivot = m[k][k].clone();

            for i in (k + 1)..n {
                for j in (k + 1)..n {
                    // m[i][j] = (pivot * m[i][j] - m[i][k] * m[k][j]) / prev_pivot
                    let num = &pivot * &m[i][j] - &m[i][k] * &m[k][j];
                    m[i][j] = Polynomial::exact_div(&num, &prev_pivot);
                }
                m[i][k] = Polynomial::zero();
            }

            prev_pivot = pivot;
        }

        let det = m[n - 1][n - 1].clone();
        if sign < 0 {
            -det
        } else {
            det
        }
    }
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        let a = &self.coeffs;
        let b = &other.coeffs;
        let len = a.len().max(b.len());
        for i in 0..len {
            let ca = a.get(i).copied().unwrap_or(0);
            let cb = b.get(i).copied().unwrap_or(0);
            if ca != cb {
                return false;
            }
        }
        true
    }
}

impl Eq for Polynomial {}

impl std::hash::Hash for Polynomial {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // Hash the trimmed representation
        let mut trimmed = self.coeffs.clone();
        while trimmed.last() == Some(&0) {
            trimmed.pop();
        }
        trimmed.hash(state);
    }
}

impl Add for &Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: &Polynomial) -> Polynomial {
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = vec![0i64; len];
        for (i, &c) in self.coeffs.iter().enumerate() {
            coeffs[i] += c;
        }
        for (i, &c) in rhs.coeffs.iter().enumerate() {
            coeffs[i] += c;
        }
        let mut result = Polynomial { coeffs };
        result.trim();
        result
    }
}

impl Sub for &Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: &Polynomial) -> Polynomial {
        let len = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = vec![0i64; len];
        for (i, &c) in self.coeffs.iter().enumerate() {
            coeffs[i] += c;
        }
        for (i, &c) in rhs.coeffs.iter().enumerate() {
            coeffs[i] -= c;
        }
        let mut result = Polynomial { coeffs };
        result.trim();
        result
    }
}

impl Mul for &Polynomial {
    type Output = Polynomial;
    fn mul(self, rhs: &Polynomial) -> Polynomial {
        if self.is_zero() || rhs.is_zero() {
            return Polynomial::zero();
        }
        let len = self.coeffs.len() + rhs.coeffs.len() - 1;
        let mut coeffs = vec![0i64; len];
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in rhs.coeffs.iter().enumerate() {
                coeffs[i + j] += a * b;
            }
        }
        let mut result = Polynomial { coeffs };
        result.trim();
        result
    }
}

impl Neg for Polynomial {
    type Output = Polynomial;
    fn neg(self) -> Polynomial {
        self.negate()
    }
}

// Convenience owned-value impls
impl Add for Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: Polynomial) -> Polynomial {
        &self + &rhs
    }
}

impl Sub for Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: Polynomial) -> Polynomial {
        &self - &rhs
    }
}

impl Mul for Polynomial {
    type Output = Polynomial;
    fn mul(self, rhs: Polynomial) -> Polynomial {
        &self * &rhs
    }
}

impl Add<&Polynomial> for Polynomial {
    type Output = Polynomial;
    fn add(self, rhs: &Polynomial) -> Polynomial {
        &self + rhs
    }
}

impl Sub<&Polynomial> for Polynomial {
    type Output = Polynomial;
    fn sub(self, rhs: &Polynomial) -> Polynomial {
        &self - rhs
    }
}

impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut first = true;
        for (i, &c) in self.coeffs.iter().enumerate() {
            if c == 0 {
                continue;
            }
            if !first && c > 0 {
                write!(f, "+")?;
            }
            match i {
                0 => write!(f, "{c}")?,
                1 => {
                    if c == 1 {
                        write!(f, "t")?;
                    } else if c == -1 {
                        write!(f, "-t")?;
                    } else {
                        write!(f, "{c}*t")?;
                    }
                }
                _ => {
                    if c == 1 {
                        write!(f, "t^{i}")?;
                    } else if c == -1 {
                        write!(f, "-t^{i}")?;
                    } else {
                        write!(f, "{c}*t^{i}")?;
                    }
                }
            }
            first = false;
        }
        Ok(())
    }
}

/// Parse a polynomial string like "-1+3*t-t^2" or "1-t+t^2-t^3+t^4".
/// Handles: implicit coefficients (t = 1*t), optional * between coeff and t,
/// ^ for exponent, constant terms, sparse polynomials.
pub fn parse_polynomial(s: &str) -> Result<Polynomial> {
    let s = s.trim();
    if s.is_empty() {
        return Ok(Polynomial::zero());
    }

    // Tokenize into terms by splitting on + or - while keeping the sign
    let mut terms: Vec<String> = Vec::new();
    let mut current = String::new();

    for (i, ch) in s.chars().enumerate() {
        if (ch == '+' || ch == '-') && i > 0 {
            if !current.is_empty() {
                terms.push(current.clone());
                current.clear();
            }
        }
        current.push(ch);
    }
    if !current.is_empty() {
        terms.push(current);
    }

    let mut coeffs: Vec<i64> = Vec::new();

    for term in &terms {
        let term = term.trim();
        if term.is_empty() {
            continue;
        }

        let (coeff, degree) = parse_term(term)?;

        while coeffs.len() <= degree {
            coeffs.push(0);
        }
        coeffs[degree] += coeff;
    }

    let mut p = Polynomial { coeffs };
    p.trim();
    Ok(p)
}

fn parse_term(term: &str) -> Result<(i64, usize)> {
    let term = term.trim();

    // Check if term contains 't'
    if let Some(t_pos) = term.find('t') {
        // Extract coefficient part (before 't')
        let coeff_part = &term[..t_pos];
        let coeff = match coeff_part.trim_end_matches('*') {
            "" | "+" => 1,
            "-" => -1,
            s => s
                .trim_end_matches('*')
                .parse::<i64>()
                .map_err(|e| KnotError::PolynomialParse(format!("bad coefficient '{s}': {e}")))?,
        };

        // Extract degree part (after 't')
        let degree_part = &term[t_pos + 1..];
        let degree = if degree_part.is_empty() {
            1
        } else if degree_part.starts_with('^') {
            degree_part[1..]
                .parse::<usize>()
                .map_err(|e| KnotError::PolynomialParse(format!("bad exponent: {e}")))?
        } else {
            return Err(KnotError::PolynomialParse(format!(
                "unexpected after t: '{degree_part}'"
            )));
        };

        Ok((coeff, degree))
    } else {
        // Constant term
        let c = term
            .parse::<i64>()
            .map_err(|e| KnotError::PolynomialParse(format!("bad constant '{term}': {e}")))?;
        Ok((c, 0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_arithmetic() {
        let a = Polynomial::t(); // t
        let b = Polynomial::one(); // 1
        let sum = &a + &b; // 1 + t
        assert_eq!(sum.coeffs, vec![1, 1]);

        let diff = &a - &b; // -1 + t
        assert_eq!(diff.coeffs, vec![-1, 1]);

        let prod = &a * &b; // t
        assert_eq!(prod.coeffs, vec![0, 1]);
    }

    #[test]
    fn test_mul() {
        // (1 - t) * (1 + t) = 1 - t^2
        let a = Polynomial {
            coeffs: vec![1, -1],
        };
        let b = Polynomial { coeffs: vec![1, 1] };
        let prod = &a * &b;
        assert_eq!(prod.coeffs, vec![1, 0, -1]);
    }

    #[test]
    fn test_normalize() {
        // t - t^2 → shift by ldegree=1 → [1, -1] → leading coeff is -1 → negate → [-1, 1]
        let p = Polynomial {
            coeffs: vec![0, 1, -1],
        };
        let n = p.normalize();
        assert_eq!(n.coeffs, vec![-1, 1]);
    }

    #[test]
    fn test_normalize_negative_leading() {
        // -1 + t → normalize → 1 - t (negate because leading coeff is negative... wait:
        // leading coeff of -1+t is 1 (positive), so no negation. Let me test: -t + t^2
        let p = Polynomial {
            coeffs: vec![0, -1, 1],
        };
        let n = p.normalize();
        // ldegree=1, shift: [-1, 1], leading is 1 (positive) → [-1, 1]
        assert_eq!(n.coeffs, vec![-1, 1]);
    }

    #[test]
    fn test_parse_trefoil() {
        let p = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(p.coeffs, vec![-1, 1, -1]);
    }

    #[test]
    fn test_parse_figure_eight() {
        let p = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(p.coeffs, vec![-1, 3, -1]);
    }

    #[test]
    fn test_parse_constant() {
        let p = parse_polynomial("1").unwrap();
        assert_eq!(p.coeffs, vec![1]);
    }

    #[test]
    fn test_parse_complex() {
        let p = parse_polynomial("1-t+t^2-t^3+t^4").unwrap();
        assert_eq!(p.coeffs, vec![1, -1, 1, -1, 1]);
    }

    #[test]
    fn test_parse_with_star() {
        let p = parse_polynomial("2-3*t+2*t^2").unwrap();
        assert_eq!(p.coeffs, vec![2, -3, 2]);
    }

    #[test]
    fn test_exact_div() {
        // (1 - t^2) / (1 - t) = (1 + t)
        let a = Polynomial {
            coeffs: vec![1, 0, -1],
        };
        let b = Polynomial {
            coeffs: vec![1, -1],
        };
        let q = Polynomial::exact_div(&a, &b);
        assert_eq!(q.coeffs, vec![1, 1]);
    }

    #[test]
    fn test_determinant_2x2() {
        // det([[1, t], [1-t, 1]]) = 1 - t*(1-t) = 1 - t + t^2
        let one = Polynomial::one();
        let t = Polynomial::t();
        let one_minus_t = Polynomial::one_minus_t();

        let matrix = vec![vec![one.clone(), t.clone()], vec![one_minus_t, one]];

        let det = Polynomial::determinant(&matrix);
        let expected = parse_polynomial("1-t+t^2").unwrap();
        assert_eq!(det, expected);
    }

    #[test]
    fn test_determinant_1x1() {
        let p = Polynomial {
            coeffs: vec![3, -2, 1],
        };
        let matrix = vec![vec![p.clone()]];
        let det = Polynomial::determinant(&matrix);
        assert_eq!(det, p);
    }

    #[test]
    fn test_eq() {
        let a = Polynomial {
            coeffs: vec![1, 2, 3],
        };
        let b = Polynomial {
            coeffs: vec![1, 2, 3, 0, 0],
        };
        assert_eq!(a, b);
    }

    #[test]
    fn test_hash_consistent() {
        use std::collections::HashMap;
        let a = Polynomial {
            coeffs: vec![1, 2, 3],
        };
        let b = Polynomial {
            coeffs: vec![1, 2, 3, 0],
        };
        let mut map = HashMap::new();
        map.insert(a, "test");
        assert_eq!(map.get(&b), Some(&"test"));
    }

    #[test]
    fn test_display() {
        let p = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(format!("{p}"), "-1+3*t-t^2");
    }

    #[test]
    fn test_parse_sparse() {
        // 1-t+t^3-t^5+t^6
        let p = parse_polynomial("1-t+t^3-t^5+t^6").unwrap();
        assert_eq!(p.coeffs, vec![1, -1, 0, 1, 0, -1, 1]);
    }

    #[test]
    fn test_parse_negative_leading_constant() {
        let p = parse_polynomial("1+t-3*t^2+t^3+t^4").unwrap();
        assert_eq!(p.coeffs, vec![1, 1, -3, 1, 1]);
    }
}
