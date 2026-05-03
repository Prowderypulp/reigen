use super::vendor::{is_header_or_comment_line, parse_csv_fields, DtcCall, VendorParser};
use anyhow::{anyhow, Result};

pub struct FtdnaParser;

impl VendorParser for FtdnaParser {
    fn parse_line(&mut self, line: &str) -> Result<Option<DtcCall>> {
        let line = line.trim();
        if is_header_or_comment_line(line) {
            return Ok(None);
        }

        let parts = parse_csv_fields(line)?;

        if parts.len() < 4 {
            return Ok(None);
        }

        let rsid = parts[0].clone();
        let chrom = parts[1].clone();
        let pos_str = &parts[2];
        if pos_str.is_empty() {
            return Ok(None);
        }
        let pos = pos_str
            .parse::<u64>()
            .map_err(|_| anyhow!("invalid pos: {}", pos_str))?;
        let geno = &parts[3];

        if geno == "--" || geno.is_empty() {
            return Ok(None);
        }

        let a1 = geno.chars().nth(0).unwrap_or(' ');
        let a2 = geno.chars().nth(1).unwrap_or(a1);

        Ok(Some(DtcCall {
            rsid,
            chrom,
            pos,
            a1,
            a2,
        }))
    }
}
