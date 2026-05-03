use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use reigen::format::Format;
use reigen::merge::{run_merge, MergeArgs};
use std::fs::File;
use std::io::{BufWriter, Write};
use tempfile::tempdir;

fn pack_record(gs: &[u8]) -> Vec<u8> {
    let n = gs.len();
    let mut out = vec![0u8; (n * 2 + 7) / 8];
    for (i, &g) in gs.iter().enumerate() {
        let two = match g {
            0 => 0b00,
            1 => 0b01,
            2 => 0b10,
            _ => 0b11,
        };
        let byte = i / 4;
        let shift = 6 - 2 * (i % 4);
        out[byte] |= two << shift;
    }
    out
}

fn create_dataset(dir: &std::path::Path, prefix: &str, nind: usize, nsnp: usize, snp_start: usize) {
    let mut ind = BufWriter::new(File::create(dir.join(format!("{}.ind", prefix))).unwrap());
    for i in 0..nind {
        writeln!(ind, "{}_S{:03} U Pop", prefix, i).unwrap();
    }
    ind.flush().unwrap();

    let mut snp = BufWriter::new(File::create(dir.join(format!("{}.snp", prefix))).unwrap());
    for i in 0..nsnp {
        let idx = snp_start + i;
        writeln!(snp, "rs{:06} 1 0.0 {} A C", idx, idx * 100 + 1).unwrap();
    }
    snp.flush().unwrap();

    let mut geno = BufWriter::new(File::create(dir.join(format!("{}.geno", prefix))).unwrap());
    let rlen = std::cmp::max(48, (nind * 2 + 7) / 8);
    let header = format!("GENO {} {} 0 0", nind, nsnp);
    let mut hbuf = vec![0u8; rlen];
    hbuf[..header.len()].copy_from_slice(header.as_bytes());
    geno.write_all(&hbuf).unwrap();

    let gs: Vec<u8> = (0..nind).map(|_| 0).collect(); // all homozygous ref
    let rec = pack_record(&gs);
    let mut rbuf = vec![0u8; rlen];
    rbuf[..rec.len()].copy_from_slice(&rec);

    for _ in 0..nsnp {
        geno.write_all(&rbuf).unwrap();
    }
    geno.flush().unwrap();
}

fn bench_merge(c: &mut Criterion) {
    let dir = tempdir().unwrap();

    let sizes = vec![1000, 10000];

    for nsnp in sizes {
        let ds1 = dir.path().join(format!("ds1_{}", nsnp));
        let ds2 = dir.path().join(format!("ds2_{}", nsnp));

        create_dataset(dir.path(), &format!("ds1_{}", nsnp), 100, nsnp, 0);
        create_dataset(dir.path(), &format!("ds2_{}", nsnp), 100, nsnp, nsnp / 2);

        c.bench_with_input(
            BenchmarkId::new("merge_intersection", nsnp),
            &nsnp,
            |b, &_n| {
                b.iter(|| {
                    let args = MergeArgs {
                        inputs: vec![
                            ds1.to_str().unwrap().to_string(),
                            ds2.to_str().unwrap().to_string(),
                        ],
                        out_format: Format::PackedAncestrymap,
                        out_prefix: None,
                        out_geno: Some(dir.path().join("out.geno")),
                        out_snp: Some(dir.path().join("out.snp")),
                        out_ind: Some(dir.path().join("out.ind")),
                        allow_ambiguous: false,
                        intersection: true,
                        flip_strand: true,
                        strict_ids: false,
                        numchrom: 22,
                        no_familynames: false,
                    };
                    run_merge(args).unwrap();
                })
            },
        );
    }
}

criterion_group!(benches, bench_merge);
criterion_main!(benches);
