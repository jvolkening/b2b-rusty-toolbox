use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::Path;

const DEFAULTS: &str = include_str!("../default_adapters.json");
const CONTEXT_NAMES: [&str; 7] = [
    "native_start_prefix",
    "native_start_suffix",
    "native_end_prefix",
    "native_end_suffix",
    "rapid_old_prefix",
    "rapid_new_prefix",
    "rapid_suffix",
];

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sequence {
    pub name: String,
    pub sequence: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Adapter {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub start: Option<Sequence>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub end: Option<Sequence>,
    #[serde(skip)]
    pub best_start: f64,
    #[serde(skip)]
    pub best_end: f64,
}

impl Adapter {
    pub fn barcode(&self) -> bool {
        self.name.starts_with("Barcode ")
    }

    pub fn barcode_direction(&self) -> &'static str {
        if self.start.as_ref().is_some_and(|s| s.name.contains("_rev")) {
            "reverse"
        } else {
            "forward"
        }
    }

    pub fn barcode_name(&self) -> String {
        let mut names = vec![self.name.as_str()];
        if let Some(seq) = &self.start {
            names.push(&seq.name);
        }
        if let Some(seq) = &self.end {
            names.push(&seq.name);
        }
        names.sort_by_key(|name| name.len());
        names[0].replace(' ', "_")
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    version: u32,
    pub adapters: Vec<Adapter>,
    pub barcode_contexts: BTreeMap<String, String>,
}

fn valid_name(name: &str) -> bool {
    !name.trim().is_empty() && !name.contains(['/', '\\', '\n', '\r', '\0'])
}

fn validate_sequence(seq: &mut Sequence) -> Result<()> {
    if !valid_name(&seq.name) || seq.sequence.is_empty() {
        bail!("adapter sequence requires a valid name and nonempty sequence");
    }
    seq.sequence.make_ascii_uppercase();
    if !seq.sequence.bytes().all(|c| b"ACGTN".contains(&c)) {
        bail!("adapter sequence must contain only A, C, G, T or N");
    }
    Ok(())
}

pub fn load(path: Option<&Path>) -> Result<Config> {
    let text = if let Some(path) = path {
        fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?
    } else {
        DEFAULTS.to_owned()
    };
    let mut config: Config = serde_json::from_str(&text).context("parsing adapter JSON")?;
    if config.version != 1 || config.adapters.is_empty() {
        bail!("adapter configuration must have version 1 and at least one adapter");
    }
    let mut names = HashSet::new();
    for adapter in &mut config.adapters {
        if !valid_name(&adapter.name) || !names.insert(adapter.name.clone()) {
            bail!("invalid or duplicate adapter name: {}", adapter.name);
        }
        if adapter.start.is_none() && adapter.end.is_none() {
            bail!("adapter {} has no sequences", adapter.name);
        }
        if let Some(seq) = &mut adapter.start {
            validate_sequence(seq)?;
        }
        if let Some(seq) = &mut adapter.end {
            validate_sequence(seq)?;
        }
        if adapter.barcode() {
            let (Some(start), Some(_)) = (&adapter.start, &adapter.end) else {
                bail!("barcode {} requires start and end sequences", adapter.name);
            };
            let forward = adapter.name.to_lowercase().contains("(forward)");
            let reverse = adapter.name.to_lowercase().contains("(reverse)");
            if forward == reverse || start.name.contains("_rev") != reverse {
                bail!("barcode {} has an inconsistent orientation", adapter.name);
            }
        }
    }
    if config.barcode_contexts.len() != CONTEXT_NAMES.len()
        || CONTEXT_NAMES
            .iter()
            .any(|key| !config.barcode_contexts.contains_key(*key))
    {
        bail!("barcode_contexts must contain all seven context fields");
    }
    for seq in config.barcode_contexts.values_mut() {
        seq.make_ascii_uppercase();
        if !seq.bytes().all(|c| b"ACGTN".contains(&c)) {
            bail!("barcode contexts must contain only A, C, G, T or N");
        }
    }
    Ok(config)
}

pub fn write_defaults(path: &Path) -> Result<()> {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .with_context(|| format!("creating {}", path.display()))?;
    file.write_all(DEFAULTS.as_bytes())?;
    Ok(())
}

pub fn expand(adapters: &mut Vec<Adapter>, contexts: &BTreeMap<String, String>) {
    let names: HashSet<String> = adapters.iter().map(|a| a.name.clone()).collect();
    let mut extra = Vec::new();
    for n in 1..=96 {
        let reverse = format!("Barcode {n} (reverse)");
        if names.contains("SQK-NSK007")
            && let Some(barcode) = adapters.iter().find(|a| a.name == reverse)
        {
            let start = &barcode.start.as_ref().unwrap().sequence;
            let end = &barcode.end.as_ref().unwrap().sequence;
            extra.push(Adapter {
                name: format!("Native barcoding {n} (full sequence)"),
                start: Some(Sequence {
                    name: format!("NB{n:02}_start"),
                    sequence: format!(
                        "{}{start}{}",
                        contexts["native_start_prefix"], contexts["native_start_suffix"]
                    ),
                }),
                end: Some(Sequence {
                    name: format!("NB{n:02}_end"),
                    sequence: format!(
                        "{}{end}{}",
                        contexts["native_end_prefix"], contexts["native_end_suffix"]
                    ),
                }),
                best_start: 0.0,
                best_end: 0.0,
            });
        }
        let forward = format!("Barcode {n} (forward)");
        if names.contains("Rapid")
            && let Some(barcode) = adapters.iter().find(|a| a.name == forward)
        {
            let variant = if names.contains("RBK004_upstream") {
                Some("new")
            } else if names.contains("SQK-NSK007") {
                Some("old")
            } else {
                None
            };
            if let Some(variant) = variant {
                let start = &barcode.start.as_ref().unwrap().sequence;
                extra.push(Adapter {
                    name: format!("Rapid barcoding {n} (full sequence, {variant})"),
                    start: Some(Sequence {
                        name: format!("RB{n:02}_full"),
                        sequence: format!(
                            "{}{start}{}",
                            contexts[&format!("rapid_{variant}_prefix")],
                            contexts["rapid_suffix"]
                        ),
                    }),
                    end: None,
                    best_start: 0.0,
                    best_end: 0.0,
                });
            }
        }
    }
    adapters.extend(extra);
}
