use std::collections::HashMap;
use std::io::BufRead;

use crate::error::{KnotError, Result};
use crate::polynomial::{parse_polynomial, Polynomial};

/// Table mapping Alexander polynomials to knot type names.
/// Stores both `poly -> name` and `-poly -> name` for O(1) lookup,
/// since the Alexander polynomial is defined only up to sign.
pub struct AlexanderTable {
    table: HashMap<Polynomial, String>,
}

impl AlexanderTable {
    /// Parse the table from a reader (file or string).
    /// Format: each line is `knot_name\tpolynomial_string` or space-separated.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        let mut table = HashMap::new();

        for line in reader.lines() {
            let line = line.map_err(KnotError::Io)?;
            let line = line.trim().to_string();
            if line.is_empty() {
                continue;
            }

            // Split by whitespace (tab or spaces)
            let parts: Vec<&str> = line.splitn(2, |c: char| c == '\t' || c == ' ').collect();
            if parts.len() < 2 {
                continue;
            }

            let knot_name = parts[0].trim().to_string();
            let poly_str = parts[1].trim();

            let poly = parse_polynomial(poly_str)?;
            let neg_poly = poly.negate();

            table.insert(poly, knot_name.clone());
            if !table.contains_key(&neg_poly) {
                table.insert(neg_poly, knot_name);
            }
        }

        Ok(AlexanderTable { table })
    }

    /// Load from a file path.
    pub fn from_file(path: &str) -> Result<Self> {
        let file = std::fs::File::open(path).map_err(KnotError::Io)?;
        let reader = std::io::BufReader::new(file);
        Self::from_reader(reader)
    }

    /// Look up a polynomial. Tries both the polynomial and its negation.
    pub fn lookup(&self, poly: &Polynomial) -> Option<&str> {
        if let Some(name) = self.table.get(poly) {
            return Some(name.as_str());
        }
        let neg = poly.negate();
        self.table.get(&neg).map(|s| s.as_str())
    }

    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn is_empty(&self) -> bool {
        self.table.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_table() {
        let data = "1\t1\n3_1\t-1+t-t^2\n4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        assert!(table.len() > 0);

        // Lookup trefoil
        let trefoil = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(table.lookup(&trefoil), Some("3_1"));

        // Lookup negated trefoil (should also work)
        let neg_trefoil = trefoil.negate();
        assert_eq!(table.lookup(&neg_trefoil), Some("3_1"));

        // Lookup unknot
        let unknot = parse_polynomial("1").unwrap();
        assert_eq!(table.lookup(&unknot), Some("1"));
    }

    #[test]
    fn test_figure_eight() {
        let data = "4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(table.lookup(&poly), Some("4_1"));
    }
}
