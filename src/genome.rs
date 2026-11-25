use clap::ValueEnum;
use std::collections::HashMap;
use std::io::{self, BufRead};

#[derive(Clone, Copy, ValueEnum)]
pub enum Vendor {
    #[value(name = "23andme")]
    TwentyThreeAndMe,
    Wegene,
    Ancestry,
    Ftdna,
    Ftdna2,
    Myheritage,
}

// Index only the markers requested by the models. Zero bytes denote a
// reference marker for which the vendor file has supplied no usable call.
#[derive(Default)]
pub struct Genome {
    rs: HashMap<u32, [u8; 2]>,
    other: HashMap<Box<str>, [u8; 2]>,
}

fn rs_number(id: &str) -> Option<u32> {
    let digits = id.strip_prefix("rs")?;
    // Preserve exact identifier matching for aliases such as rs01 or rs+1.
    if !matches!(digits.as_bytes().first(), Some(b'1'..=b'9')) {
        return None;
    }
    digits.parse().ok()
}

impl Genome {
    pub fn include(&mut self, id: &str) {
        if let Some(number) = rs_number(id) {
            self.rs.entry(number).or_default();
        } else if !self.other.contains_key(id) {
            self.other.insert(id.into(), [0; 2]);
        }
    }

    fn record(&mut self, id: &str, call: [u8; 2]) {
        let slot = match rs_number(id) {
            Some(number) => self.rs.get_mut(&number),
            None => self.other.get_mut(id),
        };
        if let Some(slot) = slot {
            *slot = call;
        }
    }

    pub fn get(&self, id: &str) -> Option<&[u8; 2]> {
        match rs_number(id) {
            Some(number) => self.rs.get(&number),
            None => self.other.get(id),
        }
    }
}

fn base(value: u8) -> bool {
    matches!(value, b'A' | b'C' | b'G' | b'T')
}

// These vendor exports store one SNP per line. Quoted CSV fields can
// contain delimiters (e.g. in an identifier or a vendor annotation).
fn field<'a>(line: &mut &'a str, delimiter: char) -> Option<&'a str> {
    if line.is_empty() {
        return None;
    }
    if let Some(quoted) = line.strip_prefix('"') {
        let mut quotes = quoted.match_indices('"').peekable();
        while let Some((end, _)) = quotes.next() {
            if quotes.next_if(|(next, _)| *next == end + 1).is_some() {
                continue;
            }
            let tail = &quoted[end + 1..];
            *line = if tail.is_empty() {
                ""
            } else {
                tail.strip_prefix(delimiter)?
            };
            return Some(&quoted[..end]);
        }
        return None;
    }
    let (value, tail) = line.split_once(delimiter).unwrap_or((line, ""));
    *line = tail;
    Some(value)
}

fn parse_row(mut row: &str, vendor: Vendor) -> Option<(&str, [u8; 2])> {
    let delimiter = match vendor {
        Vendor::Ftdna | Vendor::Ftdna2 | Vendor::Myheritage => ',',
        _ => '\t',
    };
    if row.ends_with(delimiter) {
        return None;
    }
    let id = field(&mut row, delimiter)?;
    field(&mut row, delimiter)?; // chromosome
    field(&mut row, delimiter)?; // position
    let call = field(&mut row, delimiter)?.as_bytes();
    let alleles = if matches!(vendor, Vendor::Ancestry | Vendor::Ftdna2) {
        let second = field(&mut row, delimiter)?.as_bytes();
        let ([a], [b]) = (call, second) else {
            return None;
        };
        [*a, *b]
    } else {
        match call {
            [a] => [*a, *a], // Upstream counts haploid calls twice.
            [a, b] => [*a, *b],
            _ => return None,
        }
    };
    if id.is_empty() || !row.is_empty() || !base(alleles[1]) {
        return None;
    }
    Some((id, alleles))
}

pub fn read(mut reader: impl BufRead, vendor: Vendor, genome: &mut Genome) -> io::Result<()> {
    let mut valid_input = false;
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let bytes = line.strip_prefix(b"\xef\xbb\xbf").unwrap_or(&line);
        if bytes.starts_with(b"#") {
            continue;
        }
        let row = std::str::from_utf8(bytes).map_err(io::Error::other)?;
        if let Some((id, call)) = parse_row(row.trim_end_matches(['\r', '\n']), vendor) {
            valid_input = true;
            genome.record(id, call);
        }
    }
    if !valid_input {
        return Err(io::Error::other(
            "no usable genotypes; check the input format and --vendor",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vendor_formats_and_missing_calls() {
        for (vendor, row) in [
            (Vendor::Wegene, "s\t1\t2\tAC"),
            (Vendor::TwentyThreeAndMe, "s\t1\t2\tAC"),
            (Vendor::Ancestry, "s\t1\t2\tA\tC"),
            (Vendor::Myheritage, "\"s\",\"1\",\"2\",\"AC\""),
            (Vendor::Ftdna, "s,1,2,AC"),
            (Vendor::Ftdna2, "s,1,2,A,C"),
        ] {
            assert_eq!(parse_row(row, vendor), Some(("s", *b"AC")));
        }
        let rows = b"\xef\xbb\xbf#comment\r\ns\t1\t2\tTT\r\ns\t1\t2\tA\r\ns\t1\t2\t--\r\n";
        let mut genome = Genome::default();
        genome.include("s");
        read(&rows[..], Vendor::Wegene, &mut genome).unwrap();
        assert_eq!(genome.get("s"), Some(b"AA"));
    }

    #[test]
    fn indexes_only_requested_ids_and_keeps_exact_spelling() {
        let mut genome = Genome::default();
        for id in ["rs1", "rs01", "rs4294967296", "probe-1"] {
            genome.include(id);
        }
        let rows = b"rs1\t1\t1\tAG\nrs01\t1\t1\tTT\nrs4294967296\t1\t1\tCG\nprobe-1\t1\t1\tAT\nrs2\t1\t1\tCC\n";
        read(&rows[..], Vendor::Wegene, &mut genome).unwrap();
        assert_eq!(genome.get("rs1"), Some(b"AG"));
        assert_eq!(genome.get("rs01"), Some(b"TT"));
        assert_eq!(genome.get("rs4294967296"), Some(b"CG"));
        assert_eq!(genome.get("probe-1"), Some(b"AT"));
        assert_eq!(genome.get("rs2"), None);
    }
}
