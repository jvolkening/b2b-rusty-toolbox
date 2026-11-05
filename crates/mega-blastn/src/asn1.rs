//! Minimal BER reader for the `Blast-def-line-set` records stored in `.nhr`
//! files. Only the pieces needed to print BLAST-style subject identifiers are
//! decoded (title, Seq-ids, taxid).

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdKind {
    Local,
    Gi,
    General,
    Other,
}

#[derive(Debug, Clone)]
pub struct SeqId {
    kind: IdKind,
    /// Complete FASTA-style label, e.g. `gb|AB123456.1|` or `loc1` (local ids
    /// carry no `lcl|` prefix, matching BLAST's tabular output).
    label: String,
    /// Accession without version (`sacc`).
    acc: String,
    /// Accession with version (`saccver`).
    acc_ver: String,
}

impl SeqId {
    pub fn label(&self) -> &str {
        &self.label
    }
    pub fn accession(&self) -> &str {
        &self.acc
    }
    pub fn accession_version(&self) -> &str {
        &self.acc_ver
    }
    pub fn is_gi(&self) -> bool {
        self.kind == IdKind::Gi
    }
    pub fn is_blast_ordinal(&self) -> bool {
        self.kind == IdKind::General && self.label.starts_with("gnl|BL_ORD_ID|")
    }
}

#[derive(Debug, Clone, Default)]
pub struct DefLine {
    pub title: String,
    pub ids: Vec<SeqId>,
    pub taxid: Option<i64>,
}

#[derive(Debug)]
struct Tlv<'a> {
    class: u8,
    constructed: bool,
    tag: u32,
    content: &'a [u8],
    consumed: usize,
}

fn parse_tlv(buf: &[u8]) -> Option<Tlv<'_>> {
    let b0 = *buf.first()?;
    let class = b0 >> 6;
    let constructed = b0 & 0x20 != 0;
    let mut tag = (b0 & 0x1f) as u32;
    let mut pos = 1;
    if tag == 0x1f {
        tag = 0;
        loop {
            let b = *buf.get(pos)?;
            pos += 1;
            tag = (tag << 7) | (b & 0x7f) as u32;
            if b & 0x80 == 0 {
                break;
            }
        }
    }
    let l0 = *buf.get(pos)?;
    pos += 1;
    if l0 == 0x80 {
        if !constructed {
            return None;
        }
        // indefinite length: children until end-of-contents marker
        let start = pos;
        loop {
            if buf.get(pos..pos + 2)? == [0, 0] {
                return Some(Tlv {
                    class,
                    constructed,
                    tag,
                    content: &buf[start..pos],
                    consumed: pos + 2,
                });
            }
            let child = parse_tlv(&buf[pos..])?;
            pos += child.consumed;
        }
    }
    let len = if l0 & 0x80 == 0 {
        l0 as usize
    } else {
        let n = (l0 & 0x7f) as usize;
        let mut v = 0usize;
        for _ in 0..n {
            v = (v << 8) | *buf.get(pos)? as usize;
            pos += 1;
        }
        v
    };
    let content = buf.get(pos..pos + len)?;
    Some(Tlv { class, constructed, tag, content, consumed: pos + len })
}

fn children<'a>(t: &Tlv<'a>) -> Vec<Tlv<'a>> {
    let mut out = Vec::new();
    let mut pos = 0;
    while pos < t.content.len() {
        match parse_tlv(&t.content[pos..]) {
            Some(c) => {
                pos += c.consumed;
                out.push(c);
            }
            None => break,
        }
    }
    out
}

#[derive(Debug, Clone)]
enum Leaf {
    Int(i64),
    Str(String),
}

impl Leaf {
    fn text(&self) -> String {
        match self {
            Leaf::Int(i) => i.to_string(),
            Leaf::Str(s) => s.clone(),
        }
    }
}

fn leaf(t: &Tlv<'_>) -> Option<Leaf> {
    if t.constructed {
        children(t).iter().find_map(leaf)
    } else if t.class == 0 && t.tag == 2 {
        let mut v: i64 = if t.content.first().is_some_and(|b| b & 0x80 != 0) { -1 } else { 0 };
        for b in t.content {
            v = (v << 8) | *b as i64;
        }
        Some(Leaf::Int(v))
    } else {
        Some(Leaf::Str(String::from_utf8_lossy(t.content).into_owned()))
    }
}

fn leaves(t: &Tlv<'_>, out: &mut Vec<Leaf>) {
    if t.constructed {
        for c in children(t) {
            leaves(&c, out);
        }
    } else if let Some(l) = leaf(t) {
        out.push(l);
    }
}

fn by_ctx_tag<'a>(kids: &'a [Tlv<'a>], tag: u32) -> Option<&'a Tlv<'a>> {
    kids.iter().find(|k| k.class == 2 && k.tag == tag)
}

fn decode_seqid(t: &Tlv<'_>) -> Option<SeqId> {
    if t.class != 2 {
        return None;
    }
    let inner = if t.constructed { children(t).into_iter().next()? } else { return None };
    let prefix = match t.tag {
        0 => "",
        1..=3 => "",
        4 => "gb",
        5 => "emb",
        6 => "pir",
        7 => "sp",
        8 => "pat",
        9 => "ref",
        10 => "gnl",
        11 => "gi",
        12 => "dbj",
        13 => "prf",
        14 => "pdb",
        15 => "tpg",
        16 => "tpe",
        17 => "tpd",
        18 => "gpp",
        19 => "nat",
        _ => return None,
    };
    match t.tag {
        0 => {
            let s = leaf(&inner)?.text();
            Some(SeqId { kind: IdKind::Local, label: s.clone(), acc: s.clone(), acc_ver: s })
        }
        11 => {
            let s = leaf(&inner)?.text();
            Some(SeqId {
                kind: IdKind::Gi,
                label: format!("gi|{s}"),
                acc: s.clone(),
                acc_ver: s,
            })
        }
        10 => {
            let mut v = Vec::new();
            leaves(&inner, &mut v);
            let db = v.first()?.text();
            let tag = v.get(1).map(|l| l.text()).unwrap_or_default();
            Some(SeqId {
                kind: IdKind::General,
                label: format!("gnl|{db}|{tag}"),
                acc: format!("{db}:{tag}"),
                acc_ver: format!("{db}:{tag}"),
            })
        }
        14 => {
            let kids = children(&inner);
            let mol = by_ctx_tag(&kids, 0).and_then(leaf).map(|l| l.text()).unwrap_or_default();
            let chain_id = by_ctx_tag(&kids, 3).and_then(leaf).map(|l| l.text());
            let chain = chain_id.or_else(|| {
                by_ctx_tag(&kids, 1).and_then(leaf).and_then(|l| match l {
                    Leaf::Int(i) => char::from_u32(i as u32).map(|c| c.to_string()),
                    Leaf::Str(s) => Some(s),
                })
            });
            let chain = chain.unwrap_or_default();
            Some(SeqId {
                kind: IdKind::Other,
                label: format!("pdb|{mol}|{chain}"),
                acc: format!("{mol}_{chain}"),
                acc_ver: format!("{mol}_{chain}"),
            })
        }
        1 | 2 | 3 | 8 | 19 => None,
        _ => {
            let kids = children(&inner);
            let text = |tag: u32| by_ctx_tag(&kids, tag).and_then(leaf);
            let name = text(0).map(|l| l.text()).unwrap_or_default();
            let acc = text(1).map(|l| l.text()).unwrap_or_default();
            let version = text(3).and_then(|l| match l {
                Leaf::Int(i) if i > 0 => Some(i),
                _ => None,
            });
            let acc_ver = match version {
                Some(v) if !acc.is_empty() => format!("{acc}.{v}"),
                _ => acc.clone(),
            };
            Some(SeqId {
                kind: IdKind::Other,
                label: format!("{prefix}|{acc_ver}|{name}"),
                acc: if acc.is_empty() { name.clone() } else { acc },
                acc_ver: if acc_ver.is_empty() { name } else { acc_ver },
            })
        }
    }
}

/// Parse a `Blast-def-line-set` blob.
pub fn parse_defline_set(buf: &[u8]) -> Vec<DefLine> {
    let mut out = Vec::new();
    let Some(top) = parse_tlv(buf) else { return out };
    for dl in children(&top) {
        let kids = children(&dl);
        let mut d = DefLine::default();
        if let Some(t) = by_ctx_tag(&kids, 0).and_then(leaf) {
            d.title = t.text();
        }
        if let Some(ids) = by_ctx_tag(&kids, 1) {
            for seq_of in children(ids) {
                for id in children(&seq_of) {
                    if let Some(sid) = decode_seqid(&id) {
                        d.ids.push(sid);
                    }
                }
            }
        }
        if let Some(Leaf::Int(i)) = by_ctx_tag(&kids, 2).and_then(leaf) {
            d.taxid = Some(i);
        }
        out.push(d);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_id_and_title() {
        // 30 80 | 30 80 | a0 80 1a 0e title | a1 80 30 80 a0 80 1a 04 seq0 ...
        let mut b: Vec<u8> = vec![0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x05];
        b.extend(b"hello");
        b.extend([0, 0, 0xa1, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, 0x04]);
        b.extend(b"seq0");
        b.extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        let d = parse_defline_set(&b);
        assert_eq!(d.len(), 1);
        assert_eq!(d[0].title, "hello");
        assert_eq!(d[0].ids[0].label(), "seq0");
    }
}
