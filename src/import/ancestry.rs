use super::vendor::{is_header_or_comment_line, DtcCall, VendorParser};
use anyhow::{anyhow, Result};

pub struct AncestryParser;

impl VendorParser for AncestryParser {
    fn parse_line(&mut self, line: &str) -> Result<Option<DtcCall>> {
        let line = line.trim();
        if is_header_or_comment_line(line) {
            return Ok(None);
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            return Ok(None);
        }

        let rsid = parts[0].to_string();
        let chrom = parts[1].to_string();
        let pos = parts[2]
            .parse::<u64>()
            .map_err(|_| anyhow!("invalid pos: {}", parts[2]))?;
        let a1_str = parts[3];
        let a2_str = parts[4];

        if a1_str == "0" || a2_str == "0" {
            return Ok(None);
        }

        let a1 = a1_str.chars().next().unwrap_or(' ');
        let a2 = a2_str.chars().next().unwrap_or(' ');

        Ok(Some(DtcCall {
            rsid,
            chrom,
            pos,
            a1,
            a2,
        }))
    }
}
