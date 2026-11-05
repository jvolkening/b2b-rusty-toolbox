//! Tabular output (`-outfmt 6` / `-outfmt 7`).

use crate::blastdb::{Database, Subject};
use crate::fasta::{complement, decode_base, Query};
use crate::search::{QHsp, SubjectHit};
use crate::stats::{bit_score, fmt_bits, fmt_evalue};
use std::io::{self, Write};

pub const VERSION_BANNER: &str = "BLASTN 2.17.0+";

const FIELDS: &[(&str, &str)] = &[
    ("qseqid", "query id"),
    ("qgi", "query gi"),
    ("qacc", "query acc."),
    ("qaccver", "query acc.ver"),
    ("qlen", "query length"),
    ("sseqid", "subject id"),
    ("sallseqid", "subject ids"),
    ("sgi", "subject gi"),
    ("sallgi", "subject gis"),
    ("sacc", "subject acc."),
    ("saccver", "subject acc.ver"),
    ("sallacc", "subject accs."),
    ("slen", "subject length"),
    ("qstart", "q. start"),
    ("qend", "q. end"),
    ("sstart", "s. start"),
    ("send", "s. end"),
    ("qseq", "query seq"),
    ("sseq", "subject seq"),
    ("evalue", "evalue"),
    ("bitscore", "bit score"),
    ("score", "score"),
    ("length", "alignment length"),
    ("pident", "% identity"),
    ("nident", "identical"),
    ("mismatch", "mismatches"),
    ("positive", "positives"),
    ("gapopen", "gap opens"),
    ("gaps", "gaps"),
    ("ppos", "% positives"),
    ("frames", "query/sbjct frames"),
    ("qframe", "query frame"),
    ("sframe", "sbjct frame"),
    ("btop", "BTOP"),
    ("staxid", "subject tax id"),
    ("staxids", "subject tax ids"),
    ("stitle", "subject title"),
    ("salltitles", "subject titles"),
    ("sstrand", "subject strand"),
    ("qcovs", "% query coverage per subject"),
    ("qcovhsp", "% query coverage per hsp"),
];

const STD: &str = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore";

/// Parse the argument of `-outfmt` into (format number, field names).
pub fn parse_outfmt(spec: &str) -> Result<(u32, Vec<&'static str>), String> {
    let mut parts = spec.split_whitespace();
    let n: u32 = parts
        .next()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| format!("invalid -outfmt '{spec}'"))?;
    if n != 6 && n != 7 {
        return Err(format!("only -outfmt 6 and 7 are supported (got {n})"));
    }
    let mut names = Vec::new();
    let rest: Vec<&str> = parts.collect();
    let rest = if rest.is_empty() { vec!["std"] } else { rest };
    for r in rest {
        let list: Vec<&str> = if r == "std" { STD.split(' ').collect() } else { vec![r] };
        for w in list {
            match FIELDS.iter().find(|(k, _)| *k == w) {
                Some((k, _)) => {
                    if !names.contains(k) {
                        names.push(*k)
                    }
                }
                None => return Err(format!("unsupported -outfmt field '{w}'")),
            }
        }
    }
    Ok((n, names))
}

fn header_of(name: &str) -> &'static str {
    FIELDS.iter().find(|(k, _)| *k == name).map(|x| x.1).unwrap()
}

struct SubjInfo {
    sseqid: String,
    sacc: String,
    saccver: String,
    sallseqid: String,
    sallacc: String,
    stitle: String,
    salltitles: String,
    taxids: String,
}

fn first_token(s: &str) -> &str {
    s.split_whitespace().next().unwrap_or("")
}

fn join_ids(ids: &[crate::asn1::SeqId]) -> String {
    let mut out = String::new();
    for (i, id) in ids.iter().enumerate() {
        let l = id.label();
        if i + 1 < ids.len() {
            out.push_str(l.trim_end_matches('|'));
            out.push('|');
        } else {
            out.push_str(l);
        }
    }
    out
}

fn subject_info(db: &Database, gid: u64) -> SubjInfo {
    let (v, oid) = db.locate(gid);
    let dls = db.volumes[v].deflines(oid);
    let parsed = dls.first().map(|d| d.ids.first().map(|i| !i.is_blast_ordinal()).unwrap_or(false));
    let na = || "N/A".to_string();
    let taxids = {
        let t: Vec<String> = dls.iter().filter_map(|d| d.taxid).filter(|&t| t > 0).map(|t| t.to_string()).collect();
        if t.is_empty() { "0".to_string() } else { t.join(";") }
    };
    if dls.is_empty() {
        let id = format!("gnl|BL_ORD_ID|{gid}");
        return SubjInfo {
            sseqid: id.clone(),
            sacc: id.clone(),
            saccver: id.clone(),
            sallseqid: id.clone(),
            sallacc: id,
            stitle: na(),
            salltitles: na(),
            taxids,
        };
    }
    if parsed == Some(true) {
        let d = &dls[0];
        let id0 = d.ids.iter().find(|i| !i.is_gi()).unwrap_or(&d.ids[0]);
        let title_of = |d: &crate::asn1::DefLine| if d.title.is_empty() { na() } else { d.title.clone() };
        SubjInfo {
            sseqid: join_ids(&d.ids),
            sacc: id0.accession().to_string(),
            saccver: id0.accession_version().to_string(),
            sallseqid: dls.iter().map(|d| join_ids(&d.ids)).collect::<Vec<_>>().join(";"),
            sallacc: dls
                .iter()
                .filter_map(|d| d.ids.iter().find(|i| !i.is_gi()).or(d.ids.first()))
                .map(|i| i.accession().to_string())
                .collect::<Vec<_>>()
                .join(";"),
            stitle: title_of(d),
            salltitles: dls.iter().map(title_of).collect::<Vec<_>>().join("<>"),
            taxids,
        }
    } else {
        let tok = |d: &crate::asn1::DefLine| first_token(&d.title).to_string();
        let t0 = tok(&dls[0]);
        let t0 = if t0.is_empty() { format!("gnl|BL_ORD_ID|{gid}") } else { t0 };
        SubjInfo {
            sseqid: t0.clone(),
            sacc: t0.clone(),
            saccver: t0.clone(),
            sallseqid: dls.iter().map(tok).collect::<Vec<_>>().join(";"),
            sallacc: dls.iter().map(tok).collect::<Vec<_>>().join(";"),
            stitle: if dls[0].title.is_empty() { na() } else { dls[0].title.clone() },
            salltitles: dls.iter().map(|d| d.title.clone()).collect::<Vec<_>>().join("<>"),
            taxids,
        }
    }
}

/// Aligned strings (query, subject) in plus-query orientation.
fn alignment_strings(q: &Query, h: &QHsp, sub: &Subject) -> (Vec<u8>, Vec<u8>) {
    let qs: &[u8] = if h.minus { &q.rc } else { &q.fwd };
    let (mut qi, mut si) = (h.h.qs, h.h.ss);
    let mut a = Vec::new();
    let mut b = Vec::new();
    for &(op, n) in &h.h.ops {
        for _ in 0..n {
            match op {
                0 => {
                    a.push(decode_base(qs[qi]));
                    b.push(decode_base(sub.get(si)));
                    qi += 1;
                    si += 1;
                }
                1 => {
                    a.push(b'-');
                    b.push(decode_base(sub.get(si)));
                    si += 1;
                }
                _ => {
                    a.push(decode_base(qs[qi]));
                    b.push(b'-');
                    qi += 1;
                }
            }
        }
    }
    if h.minus {
        a.reverse();
        b.reverse();
        for c in a.iter_mut().chain(b.iter_mut()) {
            if *c != b'-' {
                let code = crate::fasta::encode_base(*c).unwrap_or(19);
                *c = decode_base(complement(code));
            }
        }
    }
    (a, b)
}

fn btop(a: &[u8], b: &[u8]) -> String {
    let mut out = String::new();
    let mut run = 0usize;
    for (&x, &y) in a.iter().zip(b) {
        if x == y && x != b'-' {
            run += 1;
        } else {
            if run > 0 {
                out.push_str(&run.to_string());
                run = 0;
            }
            out.push(x as char);
            out.push(y as char);
        }
    }
    if run > 0 {
        out.push_str(&run.to_string());
    }
    out
}

fn coverage(mut ranges: Vec<(usize, usize)>) -> usize {
    ranges.sort();
    let (mut total, mut end) = (0usize, 0usize);
    for (s, e) in ranges {
        let s = s.max(end);
        if e > s {
            total += e - s;
            end = e;
        }
    }
    total
}

fn pct(n: usize, d: usize) -> u32 {
    if d == 0 {
        return 0;
    }
    // BLAST rounds to nearest below 99% and truncates above it, so that only
    // a full-length alignment reports 100
    let p = 100.0 * n as f64 / d as f64;
    (if p < 99.0 { p + 0.5 } else { p }) as u32
}

pub struct Printer<'a> {
    pub fields: Vec<&'static str>,
    pub fmt7: bool,
    pub db_arg: &'a str,
    pub db: &'a Database,
}

impl Printer<'_> {
    pub fn print_query<W: Write>(&self, out: &mut W, q: &Query, hits: &[SubjectHit]) -> io::Result<()> {
        if self.fmt7 {
            writeln!(out, "# {VERSION_BANNER}")?;
            writeln!(out, "# Query: {}", q.title)?;
            writeln!(out, "# Database: {}", self.db_arg)?;
            if !hits.is_empty() {
                let h: Vec<&str> = self.fields.iter().map(|f| header_of(f)).collect();
                writeln!(out, "# Fields: {}", h.join(", "))?;
            }
            let n: usize = hits.iter().map(|h| h.hsps.len()).sum();
            writeln!(out, "# {n} hits found")?;
        }
        let qlen = q.len();
        let need_seq = self.fields.iter().any(|f| matches!(*f, "qseq" | "sseq" | "btop"));
        for hit in hits {
            let info = subject_info(self.db, hit.gid);
            let (v, oid) = self.db.locate(hit.gid);
            let sub = if need_seq { Some(self.db.volumes[v].subject(oid)) } else { None };
            let cov = coverage(
                hit.hsps
                    .iter()
                    .map(|h| if h.minus { (qlen - h.h.qe, qlen - h.h.qs) } else { (h.h.qs, h.h.qe) })
                    .collect(),
            );
            for h in &hit.hsps {
                let mut cols: Vec<String> = Vec::with_capacity(self.fields.len());
                let (qstart, qend, sstart, send) = if h.minus {
                    (qlen - h.h.qe + 1, qlen - h.h.qs, h.h.se, h.h.ss + 1)
                } else {
                    (h.h.qs + 1, h.h.qe, h.h.ss + 1, h.h.se)
                };
                let strings = if need_seq { Some(alignment_strings(q, h, sub.as_ref().unwrap())) } else { None };
                for f in &self.fields {
                    let s = match *f {
                        "qseqid" | "qacc" | "qaccver" => q.id.clone(),
                        "qgi" | "sgi" | "sallgi" => "N/A".to_string(),
                        "qlen" => qlen.to_string(),
                        "sseqid" => info.sseqid.clone(),
                        "sallseqid" => info.sallseqid.clone(),
                        "sacc" => info.sacc.clone(),
                        "saccver" => info.saccver.clone(),
                        "sallacc" => info.sallacc.clone(),
                        "slen" => hit.subject_len.to_string(),
                        "qstart" => qstart.to_string(),
                        "qend" => qend.to_string(),
                        "sstart" => sstart.to_string(),
                        "send" => send.to_string(),
                        "qseq" => String::from_utf8_lossy(&strings.as_ref().unwrap().0).into_owned(),
                        "sseq" => String::from_utf8_lossy(&strings.as_ref().unwrap().1).into_owned(),
                        "evalue" => fmt_evalue(h.evalue),
                        "bitscore" => fmt_bits(bit_score(h.h.score)),
                        "score" => h.h.score.to_string(),
                        "length" => h.h.length.to_string(),
                        "pident" => {
                            if h.h.nident == h.h.length {
                                "100.000".to_string()
                            } else {
                                format!("{:.3}", 100.0 * h.h.nident as f64 / h.h.length as f64)
                            }
                        }
                        "nident" | "positive" => h.h.nident.to_string(),
                        "mismatch" => h.h.mismatches().to_string(),
                        "gapopen" => h.h.gap_opens.to_string(),
                        "gaps" => h.h.gaps.to_string(),
                        "ppos" => format!("{:.2}", 100.0 * h.h.nident as f64 / h.h.length as f64),
                        "frames" => format!("1/{}", if h.minus { -1 } else { 1 }),
                        "qframe" => "1".to_string(),
                        "sframe" => (if h.minus { -1 } else { 1 }).to_string(),
                        "btop" => {
                            let s = strings.as_ref().unwrap();
                            btop(&s.0, &s.1)
                        }
                        "staxid" => info.taxids.split(';').next().unwrap_or("N/A").to_string(),
                        "staxids" => info.taxids.clone(),
                        "stitle" => info.stitle.clone(),
                        "salltitles" => info.salltitles.clone(),
                        "sstrand" => (if h.minus { "minus" } else { "plus" }).to_string(),
                        "qcovs" => pct(cov, qlen).to_string(),
                        "qcovhsp" => pct(h.h.qe - h.h.qs, qlen).to_string(),
                        _ => unreachable!(),
                    };
                    cols.push(s);
                }
                writeln!(out, "{}", cols.join("\t"))?;
            }
        }
        Ok(())
    }
}
