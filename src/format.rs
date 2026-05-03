use anyhow::{anyhow, Context, Result};
use std::path::Path;
use std::str::FromStr;

/// Supported input/output formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Eigenstrat,
    PackedAncestrymap,
    Ancestrymap,
    Ped,
    PackedPed,
    Tgeno,
}

impl Format {
    pub fn default_output_extensions(self) -> (&'static str, &'static str, &'static str) {
        match self {
            Format::PackedPed => ("bed", "bim", "fam"),
            _ => ("geno", "snp", "ind"),
        }
    }
}

impl FromStr for Format {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        match s.trim().to_ascii_uppercase().as_str() {
            "EIGENSTRAT" => Ok(Format::Eigenstrat),
            "PACKEDANCESTRYMAP" => Ok(Format::PackedAncestrymap),
            "ANCESTRYMAP" => Ok(Format::Ancestrymap),
            "PED" => Ok(Format::Ped),
            "PACKEDPED" => Ok(Format::PackedPed),
            "TGENO" | "TRANSPOSE_GENO" => Ok(Format::Tgeno),
            other => Err(anyhow!("unknown outputformat: {other}")),
        }
    }
}

pub fn infer_input_format(geno: &Path) -> Result<Format> {
    let ext = geno
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    match ext.as_str() {
        "bed" => Ok(Format::PackedPed),
        "ped" => Ok(Format::Ped),
        "tgeno" => Ok(Format::Tgeno),
        "ancestrymapgeno" => Ok(Format::Ancestrymap),
        "packedancestrymapgeno" => Ok(Format::PackedAncestrymap),
        "geno" => sniff_geno_header(geno),
        other => Err(anyhow!(
            "cannot infer input format from extension {other:?} on {}",
            geno.display()
        )),
    }
}

/// Distinguish EIGENSTRAT text from PACKEDANCESTRYMAP binary by first bytes.
/// PAM starts with ASCII "GENO ", EIGENSTRAT starts with '0'/'1'/'2'/'9'.
fn sniff_geno_header(path: &Path) -> Result<Format> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)
        .with_context(|| format!("open {} for format sniff", path.display()))?;
    let mut buf = [0u8; 6];
    let n = f.read(&mut buf)?;
    if n >= 6 && &buf[..6] == b"TGENO " {
        Ok(Format::Tgeno)
    } else if n >= 5 && &buf[..5] == b"GENO " {
        Ok(Format::PackedAncestrymap)
    } else if n >= 1 && matches!(buf[0], b'0' | b'1' | b'2' | b'9') {
        Ok(Format::Eigenstrat)
    } else {
        Err(anyhow!(
            "cannot identify .geno format of {}: first bytes = {:?}",
            path.display(),
            &buf[..n]
        ))
    }
}
