use anyhow::Result;

#[derive(Debug, Clone)]
pub struct DtcCall {
    /// Vendor rsID. Kept only for diagnostics / future matching fallback;
    /// reconciliation keys on (chrom, pos) instead.
    #[allow(dead_code)]
    pub rsid: String,
    pub chrom: String,
    pub pos: u64,
    pub a1: char,
    pub a2: char,
}

pub trait VendorParser {
    /// Parse a single line from the kit file.
    /// Returns Some(call) if the line is a valid genotype call,
    /// or None if it's a comment/header/filtered-out record.
    fn parse_line(&mut self, line: &str) -> Result<Option<DtcCall>>;
}

pub fn is_header_or_comment_line(line: &str) -> bool {
    let s = line.trim();
    if s.is_empty() || s.starts_with('#') {
        return true;
    }

    let lower = s.to_ascii_lowercase();
    if lower.starts_with("rsid")
        || lower.starts_with("\"rsid\"")
        || lower.starts_with("snpid")
        || lower.starts_with("\"snpid\"")
    {
        return true;
    }

    // Fallback heuristic for vendor header rows that are not '#' prefixed.
    lower.contains("chromosome")
        && lower.contains("position")
        && (lower.contains("genotype") || lower.contains("allele1") || lower.contains("result"))
}

pub fn parse_csv_fields(line: &str) -> Result<Vec<String>> {
    let mut out = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut in_quotes = false;

    while let Some(ch) = chars.next() {
        if in_quotes {
            if ch == '"' {
                if chars.peek() == Some(&'"') {
                    field.push('"');
                    chars.next();
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(ch);
            }
        } else {
            match ch {
                '"' => in_quotes = true,
                ',' => {
                    out.push(field.trim().to_string());
                    field.clear();
                }
                _ => field.push(ch),
            }
        }
    }

    if in_quotes {
        anyhow::bail!("unterminated quoted CSV field");
    }

    out.push(field.trim().to_string());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::{is_header_or_comment_line, parse_csv_fields};

    #[test]
    fn csv_parses_quoted_fields_with_commas() {
        let fields = parse_csv_fields("\"rs1\",\"1\",\"123\",\"A,G\"").unwrap();
        assert_eq!(fields, vec!["rs1", "1", "123", "A,G"]);
    }

    #[test]
    fn csv_parses_escaped_quote() {
        let fields = parse_csv_fields("\"a\"\"b\",2").unwrap();
        assert_eq!(fields, vec!["a\"b", "2"]);
    }

    #[test]
    fn detects_unprefixed_header_lines() {
        assert!(is_header_or_comment_line(
            "rsid\tchromosome\tposition\tgenotype"
        ));
        assert!(is_header_or_comment_line("RSID,CHROMOSOME,POSITION,RESULT"));
        assert!(!is_header_or_comment_line("rs123\t1\t12345\tAG"));
    }
}
