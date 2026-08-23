use anyhow::Result;
use clap::Parser;
use rayon::prelude::*;
use rust_annotale::{find_first_ltp, is_dna, open_sequence_reader, translate};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "TALE Repeat Differences", long_about = None)]
struct Args {
    /// TALE sequences, complete DNA or AS sequences
    #[arg(short = 't', long = "tale-sequences", required = true)]
    tale_sequences: PathBuf,

    /// The output directory
    #[arg(long = "outdir", default_value = ".")]
    outdir: PathBuf,
}

fn extract_repeats_aa(sequence: &[u8]) -> Vec<Vec<u8>> {
    // Same DNA-vs-protein contract as build_families and frameshifts:
    // anything outside the DNA alphabet is treated as an amino-acid sequence.
    let aa_seq = if is_dna(sequence) { translate(sequence) } else { sequence.to_vec() };

    let mut repeats = Vec::new();
    if let Some(start) = find_first_ltp(&aa_seq) {
        let mut curr = start;
        while curr + 34 <= aa_seq.len() {
            repeats.push(aa_seq[curr..curr + 34].to_vec());
            curr += 34;
        }
    }

    repeats
}

/// Summed Levenshtein distance over aligned repeat pairs up to the shorter TALE.
fn pairwise_repeat_distance(repeats_a: &[Vec<u8>], repeats_b: &[Vec<u8>]) -> usize {
    let min_len = repeats_a.len().min(repeats_b.len());
    (0..min_len)
        .map(|i| {
            strsim::levenshtein(
                std::str::from_utf8(&repeats_a[i]).unwrap_or(""),
                std::str::from_utf8(&repeats_b[i]).unwrap_or(""),
            )
        })
        .sum()
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.outdir.exists() {
        std::fs::create_dir_all(&args.outdir)?;
    }

    let mut reader = open_sequence_reader(&args.tale_sequences)?;
    let mut tale_data = Vec::new();

    while let Some(rec) = reader.next_record()? {
        let repeats = extract_repeats_aa(&rec.seq);
        tale_data.push((rec.id.clone(), repeats));
    }

    // Compute only the upper triangle (i < j) in parallel, then mirror it —
    // the matrix is symmetric, so this halves the work of the naive version.
    let n = tale_data.len();
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            pairs.push((i, j));
        }
    }
    let distances: Vec<usize> = pairs
        .par_iter()
        .map(|&(i, j)| pairwise_repeat_distance(&tale_data[i].1, &tale_data[j].1))
        .collect();

    let mut matrix = vec![vec![0usize; n]; n];
    for (idx, &(i, j)) in pairs.iter().enumerate() {
        let dist = distances[idx];
        matrix[i][j] = dist;
        matrix[j][i] = dist;
    }

    let out_path = args.outdir.join("repdiff_matrix.tsv");
    let mut writer = BufWriter::new(File::create(&out_path)?);

    write!(writer, "TALE")?;
    for (id_b, _) in &tale_data {
        write!(writer, "\t{}", id_b)?;
    }
    writeln!(writer)?;

    for (row_idx, id_a) in tale_data.iter().map(|(id, _)| id).enumerate() {
        write!(writer, "{}", id_a)?;
        for dist in &matrix[row_idx] {
            write!(writer, "\t{}", dist)?;
        }
        writeln!(writer)?;
    }

    println!("Repeat differences calculated. Output written to {:?}", out_path);
    Ok(())
}
