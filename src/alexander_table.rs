use std::collections::HashMap;
use std::io::BufRead;

use crate::error::{KnotError, Result};
use crate::polynomial::{parse_polynomial, Polynomial};

/// Table mapping Alexander polynomials to knot type names.
///
/// Since the Alexander polynomial is not a complete invariant, one polynomial
/// may correspond to multiple knot types (e.g. `8_3` and `10_1` share
/// `4-9*t+4*t^2`). The table stores all candidates per polynomial.
///
/// - `lookup()` returns the simplest candidate (lowest crossing number) — backward compatible.
/// - `lookup_all()` returns all candidates for advanced disambiguation.
#[derive(Debug)]
pub struct AlexanderTable {
    table: HashMap<Polynomial, Vec<String>>,
}

/// Parse "3_1" -> (3, 1) for sorting by crossing number.
fn knot_sort_key(name: &str) -> (u32, u32) {
    let parts: Vec<&str> = name.split('_').collect();
    if parts.len() == 2 {
        if let (Ok(c), Ok(i)) = (parts[0].parse(), parts[1].parse()) {
            return (c, i);
        }
    }
    // unknot "1" or unparseable — sort first
    (0, 0)
}

impl AlexanderTable {
    /// Parse the table from a reader (file or string).
    /// Format: each line is `knot_name\tpolynomial_string` or space-separated.
    ///
    /// In strict mode (default), malformed lines cause an error with line number.
    /// In non-strict mode, malformed lines are silently skipped.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self> {
        Self::from_reader_impl(reader, true)
    }

    /// Parse with explicit strict flag.
    pub fn from_reader_lenient<R: BufRead>(reader: R) -> Result<Self> {
        Self::from_reader_impl(reader, false)
    }

    fn from_reader_impl<R: BufRead>(reader: R, strict: bool) -> Result<Self> {
        let mut table: HashMap<Polynomial, Vec<String>> = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line.map_err(KnotError::Io)?;
            let line = line.trim().to_string();
            if line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.splitn(2, ['\t', ' ']).collect();
            if parts.len() < 2 {
                if strict {
                    return Err(KnotError::DataParse(format!(
                        "line {}: expected 'knot_name polynomial', got '{line}'",
                        line_num + 1
                    )));
                }
                continue;
            }

            let knot_name = parts[0].trim().to_string();
            let poly_str = parts[1].trim();
            let poly = parse_polynomial(poly_str)?;

            // Insert under both poly and -poly (Alexander polynomial defined up to sign)
            table
                .entry(poly.clone())
                .or_default()
                .push(knot_name.clone());

            let neg_poly = poly.negate();
            if neg_poly != poly {
                table.entry(neg_poly).or_default().push(knot_name);
            }
        }

        // Sort each candidate list by crossing number (simplest first)
        for candidates in table.values_mut() {
            candidates.sort_by_key(|n| knot_sort_key(n));
            candidates.dedup();
        }

        Ok(AlexanderTable { table })
    }

    /// Load from a file path.
    pub fn from_file(path: &str) -> Result<Self> {
        let file = std::fs::File::open(path).map_err(KnotError::Io)?;
        let reader = std::io::BufReader::new(file);
        Self::from_reader(reader)
    }

    /// Built-in table covering all knots up to 9 crossings (86 entries).
    /// No external file needed for common knot identification.
    pub fn builtin() -> Self {
        Self::from_reader(std::io::Cursor::new(BUILTIN_TABLE))
            .expect("built-in table is known-valid")
    }

    /// Merge another table into this one. Entries from `other` are appended
    /// to existing candidate lists (no duplicates).
    pub fn merge(&mut self, other: Self) {
        for (poly, names) in other.table {
            let entry = self.table.entry(poly).or_default();
            for name in names {
                if !entry.contains(&name) {
                    entry.push(name);
                }
            }
            entry.sort_by_key(|n| knot_sort_key(n));
        }
    }

    /// Load built-in table, then merge a user-provided file on top.
    pub fn builtin_with_file(path: &str) -> Result<Self> {
        let mut table = Self::builtin();
        let extra = Self::from_file(path)?;
        table.merge(extra);
        Ok(table)
    }

    /// Look up the simplest (lowest crossing number) knot type for a polynomial.
    /// Tries both the polynomial and its negation.
    pub fn lookup(&self, poly: &Polynomial) -> Option<&str> {
        if let Some(names) = self.table.get(poly) {
            if let Some(first) = names.first() {
                return Some(first.as_str());
            }
        }
        let neg = poly.negate();
        self.table
            .get(&neg)
            .and_then(|names| names.first())
            .map(|s| s.as_str())
    }

    /// Look up all candidate knot types for a polynomial, sorted by crossing number.
    /// Returns an empty slice if not found.
    pub fn lookup_all(&self, poly: &Polynomial) -> &[String] {
        if let Some(names) = self.table.get(poly) {
            if !names.is_empty() {
                return names;
            }
        }
        let neg = poly.negate();
        self.table.get(&neg).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Number of distinct polynomials in the table (including negations).
    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn is_empty(&self) -> bool {
        self.table.is_empty()
    }
}

/// Built-in Alexander polynomial table: all prime knots up to 9 crossings.
/// Format: knot_name<TAB>polynomial, one per line.
const BUILTIN_TABLE: &str = "\
1\t1
3_1\t-1+t-t^2
4_1\t-1+3*t-t^2
5_1\t1-t+t^2-t^3+t^4
5_2\t2-3*t+2*t^2
6_1\t-2+5*t-2*t^2
6_2\t1-3*t+3*t^2-3*t^3+t^4
6_3\t1-3*t+5*t^2-3*t^3+t^4
7_1\t1-t+t^2-t^3+t^4-t^5+t^6
7_2\t3-5*t+3*t^2
7_3\t2-3*t+3*t^2-3*t^3+2*t^4
7_4\t4-7*t+4*t^2
7_5\t2-4*t+5*t^2-4*t^3+2*t^4
7_6\t1-5*t+7*t^2-5*t^3+t^4
7_7\t1-5*t+9*t^2-5*t^3+t^4
8_1\t3-7*t+3*t^2
8_2\t1-3*t+3*t^2-3*t^3+3*t^4-3*t^5+t^6
8_3\t4-9*t+4*t^2
8_4\t2-5*t+5*t^2-5*t^3+2*t^4
8_5\t1-3*t+4*t^2-5*t^3+4*t^4-3*t^5+t^6
8_6\t2-6*t+7*t^2-6*t^3+2*t^4
8_7\t1-3*t+5*t^2-5*t^3+5*t^4-3*t^5+t^6
8_8\t2-6*t+9*t^2-6*t^3+2*t^4
8_9\t1-3*t+5*t^2-7*t^3+5*t^4-3*t^5+t^6
8_10\t1-3*t+6*t^2-7*t^3+6*t^4-3*t^5+t^6
8_11\t2-7*t+9*t^2-7*t^3+2*t^4
8_12\t1-7*t+13*t^2-7*t^3+t^4
8_13\t2-7*t+11*t^2-7*t^3+2*t^4
8_14\t2-8*t+11*t^2-8*t^3+2*t^4
8_15\t3-8*t+11*t^2-8*t^3+3*t^4
8_16\t1-4*t+8*t^2-9*t^3+8*t^4-4*t^5+t^6
8_17\t1-4*t+8*t^2-11*t^3+8*t^4-4*t^5+t^6
8_18\t1-5*t+10*t^2-13*t^3+10*t^4-5*t^5+t^6
8_19\t1-t+t^3-t^5+t^6
8_20\t1-2*t+3*t^2-2*t^3+t^4
8_21\t1-4*t+5*t^2-4*t^3+t^4
9_1\t1-t+t^2-t^3+t^4-t^5+t^6-t^7+t^8
9_2\t4-7*t+4*t^2
9_3\t2-3*t+3*t^2-3*t^3+3*t^4-3*t^5+2*t^6
9_4\t3-5*t+5*t^2-5*t^3+3*t^4
9_5\t6-11*t+6*t^2
9_6\t2-4*t+5*t^2-5*t^3+5*t^4-4*t^5+2*t^6
9_7\t3-7*t+9*t^2-7*t^3+3*t^4
9_8\t2-8*t+11*t^2-8*t^3+2*t^4
9_9\t2-4*t+6*t^2-7*t^3+6*t^4-4*t^5+2*t^6
9_10\t4-8*t+9*t^2-8*t^3+4*t^4
9_11\t1-5*t+7*t^2-7*t^3+7*t^4-5*t^5+t^6
9_12\t2-9*t+13*t^2-9*t^3+2*t^4
9_13\t4-9*t+11*t^2-9*t^3+4*t^4
9_14\t2-9*t+15*t^2-9*t^3+2*t^4
9_15\t2-10*t+15*t^2-10*t^3+2*t^4
9_16\t2-5*t+8*t^2-9*t^3+8*t^4-5*t^5+2*t^6
9_17\t1-5*t+9*t^2-9*t^3+9*t^4-5*t^5+t^6
9_18\t4-10*t+13*t^2-10*t^3+4*t^4
9_19\t2-10*t+17*t^2-10*t^3+2*t^4
9_20\t1-5*t+9*t^2-11*t^3+9*t^4-5*t^5+t^6
9_21\t2-11*t+17*t^2-11*t^3+2*t^4
9_22\t1-5*t+10*t^2-11*t^3+10*t^4-5*t^5+t^6
9_23\t4-11*t+15*t^2-11*t^3+4*t^4
9_24\t1-5*t+10*t^2-13*t^3+10*t^4-5*t^5+t^6
9_25\t3-12*t+17*t^2-12*t^3+3*t^4
9_26\t1-5*t+11*t^2-13*t^3+11*t^4-5*t^5+t^6
9_27\t1-5*t+11*t^2-15*t^3+11*t^4-5*t^5+t^6
9_28\t1-5*t+12*t^2-15*t^3+12*t^4-5*t^5+t^6
9_29\t1-5*t+12*t^2-15*t^3+12*t^4-5*t^5+t^6
9_30\t1-5*t+12*t^2-17*t^3+12*t^4-5*t^5+t^6
9_31\t1-5*t+13*t^2-17*t^3+13*t^4-5*t^5+t^6
9_32\t1-6*t+14*t^2-17*t^3+14*t^4-6*t^5+t^6
9_33\t1-6*t+14*t^2-19*t^3+14*t^4-6*t^5+t^6
9_34\t1-6*t+16*t^2-23*t^3+16*t^4-6*t^5+t^6
9_35\t7-13*t+7*t^2
9_36\t1-5*t+8*t^2-9*t^3+8*t^4-5*t^5+t^6
9_37\t2-11*t+19*t^2-11*t^3+2*t^4
9_38\t5-14*t+19*t^2-14*t^3+5*t^4
9_39\t3-14*t+21*t^2-14*t^3+3*t^4
9_40\t1-7*t+18*t^2-23*t^3+18*t^4-7*t^5+t^6
9_41\t3-12*t+19*t^2-12*t^3+3*t^4
9_42\t1-2*t+t^2-2*t^3+t^4
9_43\t1-3*t+2*t^2-t^3+2*t^4-3*t^5+t^6
9_44\t1-4*t+7*t^2-4*t^3+t^4
9_45\t1-6*t+9*t^2-6*t^3+t^4
9_46\t2-5*t+2*t^2
9_47\t1-4*t+6*t^2-5*t^3+6*t^4-4*t^5+t^6
9_48\t1-7*t+11*t^2-7*t^3+t^4
9_49\t3-6*t+7*t^2-6*t^3+3*t^4
";

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_table() {
        let data = "1\t1\n3_1\t-1+t-t^2\n4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        assert!(!table.is_empty());

        let trefoil = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(table.lookup(&trefoil), Some("3_1"));

        let neg_trefoil = trefoil.negate();
        assert_eq!(table.lookup(&neg_trefoil), Some("3_1"));

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

    #[test]
    fn test_ambiguous_polynomial() {
        // Two knots share the same polynomial
        let data = "8_3\t4-9*t+4*t^2\n10_1\t4-9*t+4*t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("4-9*t+4*t^2").unwrap();

        // lookup returns simplest (lowest crossing number)
        assert_eq!(table.lookup(&poly), Some("8_3"));

        // lookup_all returns both, sorted by crossing number
        let all = table.lookup_all(&poly);
        assert_eq!(all, &["8_3", "10_1"]);
    }

    #[test]
    fn test_lookup_all_unique() {
        let data = "3_1\t-1+t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("-1+t-t^2").unwrap();

        assert_eq!(table.lookup_all(&poly), &["3_1"]);
    }

    #[test]
    fn test_lookup_all_not_found() {
        let data = "3_1\t-1+t-t^2\n";
        let table = AlexanderTable::from_reader(data.as_bytes()).unwrap();
        let poly = parse_polynomial("999+t").unwrap();

        assert!(table.lookup_all(&poly).is_empty());
    }

    #[test]
    fn test_strict_mode_rejects_bad_line() {
        let data = "3_1\t-1+t-t^2\nbadline\n4_1\t-1+3*t-t^2\n";
        let result = AlexanderTable::from_reader(data.as_bytes());
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("line 2"));
    }

    #[test]
    fn test_lenient_mode_skips_bad_line() {
        let data = "3_1\t-1+t-t^2\nbadline\n4_1\t-1+3*t-t^2\n";
        let table = AlexanderTable::from_reader_lenient(data.as_bytes()).unwrap();
        let trefoil = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(table.lookup(&trefoil), Some("3_1"));
        let fig8 = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(table.lookup(&fig8), Some("4_1"));
    }

    #[test]
    fn test_builtin_table() {
        let table = AlexanderTable::builtin();
        // Should contain unknot, trefoil, figure-eight, and 9_49
        let unknot = parse_polynomial("1").unwrap();
        assert_eq!(table.lookup(&unknot), Some("1"));

        let trefoil = parse_polynomial("-1+t-t^2").unwrap();
        assert_eq!(table.lookup(&trefoil), Some("3_1"));

        let fig8 = parse_polynomial("-1+3*t-t^2").unwrap();
        assert_eq!(table.lookup(&fig8), Some("4_1"));

        let k9_49 = parse_polynomial("3-6*t+7*t^2-6*t^3+3*t^4").unwrap();
        assert_eq!(table.lookup(&k9_49), Some("9_49"));
    }

    #[test]
    fn test_merge_tables() {
        let mut base = AlexanderTable::builtin();
        let extra_data = "10_1\t4-9*t+4*t^2\n";
        let extra = AlexanderTable::from_reader(extra_data.as_bytes()).unwrap();
        base.merge(extra);

        // 8_3 and 10_1 share the same polynomial; 8_3 should come first
        let poly = parse_polynomial("4-9*t+4*t^2").unwrap();
        assert_eq!(base.lookup(&poly), Some("8_3"));
        let all = base.lookup_all(&poly);
        assert!(all.contains(&"10_1".to_string()));
    }
}
