//! Reader for NCBI BLAST nucleotide databases (format versions 4 and 5).
//!
//! All files are memory-mapped, so resident memory is bounded by whatever the
//! OS keeps in its page cache rather than by the database size.

use crate::asn1::{self, DefLine};
use memmap2::Mmap;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

fn err<T>(msg: impl Into<String>) -> io::Result<T> {
    Err(io::Error::new(io::ErrorKind::InvalidData, msg.into()))
}

#[inline]
fn be32(b: &[u8], off: usize) -> u32 {
    u32::from_be_bytes([b[off], b[off + 1], b[off + 2], b[off + 3]])
}

/// A run of ambiguous residues inside a sequence.
#[derive(Debug, Clone, Copy)]
struct AmbRun {
    start: u32,
    end: u32,
    code: u8,
}

/// Map an NCBI4na residue mask to the internal base code (0-3 = ACGT,
/// otherwise 4 + mask so that distinct IUPAC codes remain distinct).
#[inline]
pub fn code_from_na4(mask: u8) -> u8 {
    match mask {
        1 => 0,
        2 => 1,
        4 => 2,
        8 => 3,
        m => 4 + m,
    }
}

/// A view of one database sequence.
pub struct Subject<'a> {
    packed: &'a [u8],
    len: usize,
    amb: Vec<AmbRun>,
}

impl<'a> Subject<'a> {
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn packed(&self) -> &'a [u8] {
        self.packed
    }

    pub fn has_ambiguity(&self) -> bool {
        !self.amb.is_empty()
    }

    /// Base code at `pos`: 0-3 for ACGT, >= 4 for ambiguity codes.
    #[inline(always)]
    pub fn get(&self, pos: usize) -> u8 {
        if !self.amb.is_empty() {
            let i = self.amb.partition_point(|r| r.end as usize <= pos);
            if let Some(r) = self.amb.get(i)
                && r.start as usize <= pos {
                    return r.code;
                }
        }
        (self.packed[pos >> 2] >> (6 - 2 * (pos & 3))) & 3
    }
}

pub struct Volume {
    nin: Mmap,
    nsq: Mmap,
    nhr: Option<Mmap>,
    num_oids: u32,
    pub title: String,
    pub length: u64,
    hdr_off: usize,
    seq_off: usize,
    amb_off: usize,
}

fn map(path: &Path) -> io::Result<Mmap> {
    let f = File::open(path)
        .map_err(|e| io::Error::new(e.kind(), format!("{}: {e}", path.display())))?;
    // SAFETY: the databases are treated as read-only inputs.
    unsafe { Mmap::map(&f) }
}

impl Volume {
    fn open(base: &Path) -> io::Result<Volume> {
        let with_ext = |e: &str| {
            let mut s = base.as_os_str().to_owned();
            s.push(".");
            s.push(e);
            PathBuf::from(s)
        };
        let nin = map(&with_ext("nin"))?;
        let nsq = map(&with_ext("nsq"))?;
        let nhr = map(&with_ext("nhr")).ok();
        #[cfg(unix)]
        {
            let _ = nsq.advise(memmap2::Advice::Sequential);
        }
        let b = &nin[..];
        if b.len() < 16 {
            return err(format!("{}: truncated index", base.display()));
        }
        let version = be32(b, 0);
        if version != 4 && version != 5 {
            return err(format!("{}: unsupported BLAST DB version {version}", base.display()));
        }
        if be32(b, 4) != 0 {
            return err(format!("{}: not a nucleotide database", base.display()));
        }
        let mut off = 8;
        if version == 5 {
            off += 4; // volume number
        }
        let read_str = |off: &mut usize| -> io::Result<String> {
            if *off + 4 > b.len() {
                return err("truncated index");
            }
            let n = be32(b, *off) as usize;
            *off += 4;
            if *off + n > b.len() {
                return err("truncated index");
            }
            let s = String::from_utf8_lossy(&b[*off..*off + n]).into_owned();
            *off += n;
            Ok(s)
        };
        let title = read_str(&mut off)?;
        if version == 5 {
            read_str(&mut off)?; // LMDB file name
        }
        read_str(&mut off)?; // date
        if off + 16 > b.len() {
            return err("truncated index");
        }
        let num_oids = be32(b, off);
        let length = u64::from_le_bytes(b[off + 4..off + 12].try_into().unwrap());
        off += 16; // + max sequence length
        let region = 4 * (num_oids as usize + 1);
        if off + 3 * region > b.len() {
            return err("truncated index offsets");
        }
        Ok(Volume {
            nin,
            nsq,
            nhr,
            num_oids,
            title,
            length,
            hdr_off: off,
            seq_off: off + region,
            amb_off: off + 2 * region,
        })
    }

    pub fn num_oids(&self) -> u32 {
        self.num_oids
    }

    pub fn subject(&self, oid: u32) -> Subject<'_> {
        let i = oid as usize;
        let n = &self.nin[..];
        let start = be32(n, self.seq_off + 4 * i) as usize;
        let amb_start = be32(n, self.amb_off + 4 * i) as usize;
        let end = be32(n, self.seq_off + 4 * (i + 1)) as usize;
        let packed = &self.nsq[start..amb_start];
        let len = if packed.is_empty() {
            0
        } else {
            (packed.len() - 1) * 4 + (packed[packed.len() - 1] & 3) as usize
        };
        let mut amb = Vec::new();
        if end > amb_start + 4 {
            let a = &self.nsq[amb_start..end];
            let hdr = be32(a, 0);
            let new_fmt = hdr & 0x8000_0000 != 0;
            let count = (hdr & 0x7fff_ffff) as usize;
            let mut k = 1;
            while k <= count && 4 * k + 4 <= a.len() {
                let w = be32(a, 4 * k);
                let code = code_from_na4(((w >> 28) & 0xf) as u8);
                let (rlen, pos) = if new_fmt {
                    k += 1;
                    if 4 * k + 4 > a.len() {
                        break;
                    }
                    (((w >> 16) & 0xfff) + 1, be32(a, 4 * k))
                } else {
                    (((w >> 24) & 0xf) + 1, w & 0xff_ffff)
                };
                k += 1;
                amb.push(AmbRun { start: pos, end: pos + rlen, code });
            }
            amb.sort_by_key(|r| r.start);
        }
        Subject { packed, len, amb }
    }

    pub fn seq_len(&self, oid: u32) -> usize {
        self.subject(oid).len()
    }

    pub fn deflines(&self, oid: u32) -> Vec<DefLine> {
        let Some(h) = &self.nhr else { return Vec::new() };
        let i = oid as usize;
        let s = be32(&self.nin, self.hdr_off + 4 * i) as usize;
        let e = be32(&self.nin, self.hdr_off + 4 * (i + 1)) as usize;
        if e > h.len() || s > e {
            return Vec::new();
        }
        asn1::parse_defline_set(&h[s..e])
    }
}

pub struct Database {
    pub volumes: Vec<Volume>,
    /// First global OID of each volume.
    starts: Vec<u64>,
    pub total_seqs: u64,
    pub total_len: u64,
}

fn blastdb_dirs() -> Vec<PathBuf> {
    let mut v = vec![PathBuf::from(".")];
    if let Some(p) = std::env::var_os("BLASTDB") {
        v.extend(std::env::split_paths(&p));
    }
    v
}

fn resolve(name: &str) -> io::Result<PathBuf> {
    let p = Path::new(name);
    let has = |base: &Path| {
        let mut s = base.as_os_str().to_owned();
        s.push(".nin");
        let mut a = base.as_os_str().to_owned();
        a.push(".nal");
        Path::new(&s).exists() || Path::new(&a).exists()
    };
    if p.is_absolute() || name.contains('/') {
        if has(p) {
            return Ok(p.to_path_buf());
        }
    } else {
        for d in blastdb_dirs() {
            let c = d.join(name);
            if has(&c) {
                return Ok(c);
            }
        }
    }
    err(format!("BLAST database '{name}' not found (searched cwd and $BLASTDB)"))
}

struct Alias {
    dblist: Vec<String>,
    nseq: Option<u64>,
    length: Option<u64>,
}

fn parse_alias(path: &Path) -> io::Result<Alias> {
    let text = std::fs::read_to_string(path)?;
    let mut a = Alias { dblist: Vec::new(), nseq: None, length: None };
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (key, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        match key {
            "DBLIST" => {
                let mut in_q = false;
                let mut cur = String::new();
                for ch in rest.chars() {
                    match ch {
                        '"' => {
                            if in_q {
                                a.dblist.push(std::mem::take(&mut cur));
                            }
                            in_q = !in_q;
                        }
                        c if in_q => cur.push(c),
                        c if c.is_whitespace() => {
                            if !cur.is_empty() {
                                a.dblist.push(std::mem::take(&mut cur));
                            }
                        }
                        c => cur.push(c),
                    }
                }
                if !cur.is_empty() {
                    a.dblist.push(cur);
                }
            }
            "NSEQ" => a.nseq = rest.trim().parse().ok(),
            "LENGTH" => a.length = rest.trim().parse().ok(),
            "TITLE" | "DBTYPE" | "VERSION" | "CREATE-DATE" | "UPDATE-DATE" | "MEMB_BIT" => {}
            "OIDLIST" | "GILIST" | "SEQIDLIST" | "TAXIDLIST" | "NEGATIVE_SEQIDLIST"
            | "NEGATIVE_TAXIDLIST" | "MASKLIST" | "FIRST_OID" | "LAST_OID" => {
                return err(format!(
                    "{}: alias restriction '{key}' is not supported",
                    path.display()
                ));
            }
            _ => {}
        }
    }
    Ok(a)
}

fn collect(base: &Path, out: &mut Vec<PathBuf>, totals: &mut (Option<u64>, Option<u64>), top: bool) -> io::Result<()> {
    let mut nal = base.as_os_str().to_owned();
    nal.push(".nal");
    let nal = PathBuf::from(nal);
    let mut nin = base.as_os_str().to_owned();
    nin.push(".nin");
    // a volume with an index file wins over an alias of the same name
    if nal.exists() && !Path::new(&nin).exists() {
        let alias = parse_alias(&nal)?;
        if top {
            *totals = (alias.nseq, alias.length);
        }
        let dir = nal.parent().map(Path::to_path_buf).unwrap_or_default();
        for entry in alias.dblist {
            let p = Path::new(&entry);
            let child = if p.is_absolute() { p.to_path_buf() } else { dir.join(p) };
            collect(&child, out, totals, false)?;
        }
        Ok(())
    } else {
        out.push(base.to_path_buf());
        Ok(())
    }
}

impl Database {
    pub fn open(name: &str) -> io::Result<Database> {
        let base = resolve(name)?;
        let mut paths = Vec::new();
        let mut totals = (None, None);
        collect(&base, &mut paths, &mut totals, true)?;
        let mut volumes = Vec::with_capacity(paths.len());
        let mut starts = Vec::with_capacity(paths.len());
        let mut n = 0u64;
        let mut len = 0u64;
        for p in &paths {
            let v = Volume::open(p)?;
            starts.push(n);
            n += v.num_oids as u64;
            len += v.length;
            volumes.push(v);
        }
        Ok(Database {
            volumes,
            starts,
            total_seqs: totals.0.unwrap_or(n),
            total_len: totals.1.unwrap_or(len),
        })
    }

    pub fn num_oids(&self) -> u64 {
        self.starts.last().map(|s| s + self.volumes.last().unwrap().num_oids as u64).unwrap_or(0)
    }

    pub fn volume_start(&self, v: usize) -> u64 {
        self.starts[v]
    }

    /// Resolve a global OID to (volume index, local OID).
    pub fn locate(&self, gid: u64) -> (usize, u32) {
        let v = self.starts.partition_point(|s| *s <= gid) - 1;
        (v, (gid - self.starts[v]) as u32)
    }
}
