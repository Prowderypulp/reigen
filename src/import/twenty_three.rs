use super::vendor::{is_header_or_comment_line, DtcCall, VendorParser};
use anyhow::{anyhow, Result};

pub struct TwentyThreeAndMeParser;

impl VendorParser for TwentyThreeAndMeParser {
    fn parse_line(&mut self, line: &str) -> Result<Option<DtcCall>> {
        let line = line.trim();
        if is_header_or_comment_line(line) {
            return Ok(None);
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 {
            return Ok(None);
        }

        let rsid = parts[0].to_string();
        let chrom = parts[1].to_string();
        let pos = parts[2]
            .parse::<u64>()
            .map_err(|_| anyhow!("invalid pos: {}", parts[2]))?;
        let geno = parts[3];

        if geno == "--" || geno == "???" || geno.is_empty() {
            return Ok(None);
        }

        let a1 = geno.chars().nth(0).unwrap_or(' ');
        let a2 = geno.chars().nth(1).unwrap_or(a1); // Handle haploid as homozygous

        Ok(Some(DtcCall {
            rsid,
            chrom,
            pos,
            a1,
            a2,
        }))
    }
}
