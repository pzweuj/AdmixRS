use crate::genome::Genome;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

pub struct Model {
    pub name: &'static str,
    pub labels: Vec<(&'static str, &'static str)>,
}

pub fn catalog() -> Vec<Model> {
    let mut models = Vec::<Model>::new();
    for line in include_str!("populations.tsv").lines() {
        let mut cols = line.split('\t');
        let name = cols.next().unwrap();
        if models.last().is_none_or(|m| m.name != name) {
            models.push(Model {
                name,
                labels: Vec::new(),
            });
        }
        models
            .last_mut()
            .unwrap()
            .labels
            .push((cols.next().unwrap(), cols.next().unwrap()));
    }
    models
}

pub struct Observations {
    pub frequencies: Vec<f64>,
    // Counts for allele column 3 (F) and column 2 (1 - F), in upstream order.
    pub counts: Vec<[u8; 2]>,
    pub components: usize,
    pub total_markers: usize,
}

impl Model {
    fn open(&self, dir: &Path, suffix: &str) -> io::Result<BufReader<File>> {
        let path = dir.join(format!("{}.{suffix}", self.name));
        File::open(&path)
            .map(BufReader::new)
            .map_err(|e| io::Error::other(format!("{}: {e}", path.display())))
    }

    pub fn index_markers(&self, dir: &Path, genome: &mut Genome) -> io::Result<()> {
        let mut alleles = self.open(dir, "alleles")?;
        let mut line = String::new();
        while alleles.read_line(&mut line)? != 0 {
            if let Some(id) = line.split_ascii_whitespace().next() {
                genome.include(id);
            }
            line.clear();
        }
        Ok(())
    }

    pub fn load(&self, dir: &Path, genome: &Genome) -> io::Result<Observations> {
        let k = self.labels.len();
        let alleles = self.open(dir, "alleles")?;
        let frequencies = self.open(dir, &format!("{k}.F"))?;
        load_rows(alleles, frequencies, k, genome)
    }
}

fn load_rows(
    mut alleles: impl BufRead,
    mut frequencies: impl BufRead,
    k: usize,
    genome: &Genome,
) -> io::Result<Observations> {
    let mut result = Observations {
        frequencies: Vec::new(),
        counts: Vec::new(),
        components: k,
        total_markers: 0,
    };
    let (mut marker, mut row) = (String::new(), String::new());
    loop {
        marker.clear();
        row.clear();
        match (
            alleles.read_line(&mut marker)?,
            frequencies.read_line(&mut row)?,
        ) {
            (0, 0) => break,
            (0, _) | (_, 0) => {
                return Err(io::Error::other(
                    "allele and frequency files have different row counts",
                ));
            }
            _ => {}
        }
        result.total_markers += 1;
        let invalid = || {
            io::Error::other(format!(
                "invalid reference data at row {}",
                result.total_markers
            ))
        };
        let mut fields = marker.split_ascii_whitespace();
        let id = fields.next().ok_or_else(invalid)?;
        let allele = |s: Option<&str>| match s?.as_bytes() {
            [b @ (b'A' | b'C' | b'G' | b'T')] => Some(*b),
            _ => None,
        };
        let minor = allele(fields.next()).ok_or_else(invalid)?;
        let major = allele(fields.next()).ok_or_else(invalid)?;
        if minor == major || fields.next().is_some() {
            return Err(invalid());
        }
        let counts = genome
            .get(id)
            .map(|call| {
                [major, minor]
                    .map(|allele| u8::from(call[0] == allele) + u8::from(call[1] == allele))
            })
            .unwrap_or_default();
        let matched = counts != [0, 0];
        let mut cols = row.split_ascii_whitespace();
        for _ in 0..k {
            let f = cols
                .next()
                .and_then(|v| v.parse::<f64>().ok())
                .filter(|f| (0.0..=1.0).contains(f))
                .ok_or_else(invalid)?;
            if matched {
                result.frequencies.push(f.clamp(1e-12, 1.0 - 1e-12));
            }
        }
        if cols.next().is_some() {
            return Err(invalid());
        }
        if matched {
            result.counts.push(counts);
        }
    }
    if result.counts.is_empty() {
        return Err(io::Error::other("no matching SNP alleles for this model"));
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn orientation_and_row_alignment() {
        let mut genome = Genome::default();
        for id in ["a", "b", "c"] {
            genome.include(id);
        }
        crate::genome::read(
            &b"b\t1\t1\tGG\nc\t1\t2\tAC\n"[..],
            crate::genome::Vendor::Wegene,
            &mut genome,
        )
        .unwrap();
        let markers = &b"a T A\nb A G\nc A G\n"[..];
        let freq = &b".2 .4\n.9 .1\n.3 .7\n"[..];
        let obs = load_rows(markers, freq, 2, &genome).unwrap();
        assert_eq!(obs.counts, [[2, 0], [0, 1]]);
        assert_eq!(obs.frequencies, [0.9, 0.1, 0.3, 0.7]);
        assert_eq!(obs.total_markers, 3);
        assert!(load_rows(markers, &freq[6..], 2, &genome).is_err());
        assert!(load_rows(markers, freq, 3, &genome).is_err());
        assert!(load_rows(markers, freq, 2, &Genome::default()).is_err());
    }
}
