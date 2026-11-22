use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Denominator {
    All,
    Unambiguous,
    Assigned,
}

impl Denominator {
    pub fn parse(value: &str) -> Option<Self> {
        match value {
            "all" => Some(Self::All),
            "unambiguous" => Some(Self::Unambiguous),
            "assigned" => Some(Self::Assigned),
            _ => None,
        }
    }
}

pub type Taxonomy = HashMap<i64, i64>;

fn io_error(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

pub fn load_taxonomy<R: BufRead>(mut reader: R) -> io::Result<Taxonomy> {
    let mut taxonomy = HashMap::new();
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            continue;
        }
        let taxid: i64 = fields[0]
            .parse()
            .map_err(|_| io_error(format!("Invalid taxID in taxonomy line: {}", line.trim())))?;
        let parent: i64 = fields[1].parse().map_err(|_| {
            io_error(format!(
                "Invalid parent taxID in taxonomy line: {}",
                line.trim()
            ))
        })?;
        if parent != taxid {
            taxonomy.insert(taxid, parent);
        }
    }
    Ok(taxonomy)
}

pub fn read_taxonomy(taxdb: &Path) -> io::Result<Taxonomy> {
    let path = taxdb;
    let file = std::fs::File::open(&path)
        .map_err(|e| io::Error::new(e.kind(), format!("Error opening {}: {e}", path.display())))?;
    load_taxonomy(io::BufReader::new(file))
}

pub fn find_taxdb(supplied: Option<&str>) -> Result<PathBuf, String> {
    let supplied = match supplied {
        Some(value) => value.to_string(),
        None => std::env::var("KRAKEN_DEFAULT_DB")
            .map_err(|_| "Must specify DB with either --db or $KRAKEN_DEFAULT_DB".to_string())?,
    };

    let taxdb = if supplied.contains('/') {
        PathBuf::from(&supplied)
    } else if Path::new(&supplied).is_file() {
        PathBuf::from(&supplied)
    } else {
        let db_path = std::env::var("KRAKEN_DB_PATH").unwrap_or_else(|_| ".".to_string());
        let mut normalized = db_path;
        if normalized.starts_with(':') {
            normalized.insert(1, '.');
        }
        if normalized.ends_with(':') {
            normalized.push('.');
        }
        normalized = normalized.replace("::", ":.:");

        let mut found = None;
        for dir in normalized.split(':') {
            let candidate = Path::new(dir).join(&supplied);
            if candidate.is_file() {
                found = Some(candidate);
                break;
            }
        }
        found.ok_or_else(|| {
            let printed = std::env::var("KRAKEN_DB_PATH")
                .map(|v| format!("\"{v}\""))
                .unwrap_or_else(|_| "undefined".to_string());
            format!("unable to find {supplied} in $KRAKEN_DB_PATH ({printed})")
        })?
    };

    if !taxdb.is_file() {
        return Err(format!("taxDB file \"{}\" does not exist", taxdb.display()));
    }
    Ok(taxdb)
}

fn add_hit_sums(hit_sums: &mut HashMap<i64, u64>, taxonomy: &Taxonomy, taxid: i64, count: u64) {
    let mut current = taxid;
    loop {
        *hit_sums.entry(current).or_default() += count;
        match taxonomy.get(&current) {
            Some(&parent) => current = parent,
            None => break,
        }
    }
}

pub fn filter_line(
    line: &str,
    taxonomy: &Taxonomy,
    threshold: f64,
    denominator: Denominator,
) -> io::Result<String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Ok(line.to_string());
    }

    let seqid = fields[1];
    let called_taxon: i64 = fields[2].parse().map_err(|_| {
        io_error(format!(
            "Invalid called taxon {} in Kraken output line",
            fields[2]
        ))
    })?;
    let len = fields[3];
    let hit_list = fields[4];

    let mut hit_counts: HashMap<&str, u64> = HashMap::new();
    for hit in hit_list.split_whitespace() {
        let (taxid, count) = hit.split_once(':').unwrap_or((hit, "1"));
        if let Ok(count) = count.parse::<u64>() {
            *hit_counts.entry(taxid).or_default() += count;
        }
    }

    let mut hit_sums: HashMap<i64, u64> = HashMap::new();
    let mut total_kmers = 0u64;
    let mut total_unambig = 0u64;
    let mut total_hits = 0u64;

    for (&taxid, &count) in &hit_counts {
        total_kmers += count;
        if taxid != "A" {
            total_unambig += count;
            if let Ok(taxid) = taxid.parse::<i64>() {
                if taxid > 0 {
                    total_hits += count;
                    add_hit_sums(&mut hit_sums, taxonomy, taxid, count);
                }
            }
        }
    }

    let denominator_value = match denominator {
        Denominator::All => total_kmers,
        Denominator::Unambiguous => total_unambig,
        Denominator::Assigned => total_hits,
    };

    let mut pct = 0.0;
    let mut new_taxon = called_taxon;
    while taxonomy.contains_key(&new_taxon) {
        pct = if denominator_value == 0 {
            0.0
        } else {
            *hit_sums.get(&new_taxon).unwrap_or(&0) as f64 / denominator_value as f64
        };
        if pct >= threshold - 1e-5 {
            break;
        }
        new_taxon = taxonomy[&new_taxon];
    }

    let code = if new_taxon > 0 { "C" } else { "U" };
    Ok(format!(
        "{code}\t{seqid}\t{new_taxon}\t{len}\tP={pct:.3}\t{hit_list}"
    ))
}

pub fn filter<R: BufRead, W: Write>(
    mut reader: R,
    mut writer: W,
    taxonomy: &Taxonomy,
    threshold: f64,
    denominator: Denominator,
) -> io::Result<()> {
    let mut line = Vec::new();
    while reader.read_until(b'\n', &mut line)? != 0 {
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        let text = std::str::from_utf8(&line)
            .map_err(|_| io_error("Kraken output contains invalid UTF-8"))?;
        writeln!(
            writer,
            "{}",
            filter_line(text, taxonomy, threshold, denominator)?
        )?;
        line.clear();
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn taxonomy() -> Taxonomy {
        HashMap::from([(2, 1), (3, 1), (4, 2)])
    }

    #[test]
    fn loads_taxonomy_and_climbs_to_threshold() {
        let taxonomy = load_taxonomy("1\t1\n2\t1\n3\t1\n4\t2\n".as_bytes()).unwrap();
        let line = "C\tread1\t4\t100\t2:5 1:5";
        let out = filter_line(line, &taxonomy, 0.5, Denominator::Unambiguous).unwrap();
        assert_eq!(out, "C\tread1\t2\t100\tP=0.500\t2:5 1:5");
    }

    #[test]
    fn assigned_denominator_excludes_unassigned() {
        let taxonomy = taxonomy();
        let line = "C\tread1\t4\t100\t4:5 A:5";
        let out = filter_line(line, &taxonomy, 0.5, Denominator::Assigned).unwrap();
        assert_eq!(out, "C\tread1\t4\t100\tP=1.000\t4:5 A:5");
    }

    #[test]
    fn climbs_called_taxon_until_proportion_is_met() {
        let taxonomy = taxonomy();
        let line = "C\tread1\t4\t100\t2:3 3:7";
        let out = filter_line(line, &taxonomy, 0.8, Denominator::All).unwrap();
        assert_eq!(out, "C\tread1\t1\t100\tP=0.300\t2:3 3:7");
    }

    #[test]
    fn malformed_line_is_passed_through() {
        let out = filter_line("short\tline", &taxonomy(), 0.5, Denominator::All).unwrap();
        assert_eq!(out, "short\tline");
    }
}
