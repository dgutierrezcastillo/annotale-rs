use anyhow::{Context, Result};
use clap::Parser;
use rust_annotale::{dna_to_rvds, is_dna, open_sequence_reader, parse_rvd_sequence};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "TALE Frameshift and Truncation Scanner", long_about = None)]
struct Args {
    /// Input FASTA file (can be DNA sequences or RVD sequences)
    #[arg(short = 'i', long = "input", required = true)]
    input: PathBuf,

    /// Output file path (optional, writes to stdout if not provided)
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut reader = open_sequence_reader(&args.input)
        .with_context(|| format!("Failed to read input FASTA file {:?}", args.input))?;

    let mut tales: Vec<(String, Vec<String>)> = Vec::new();

    while let Some(rec) = reader.next_record()? {
        let id = rec.id.clone();
        let seq = &rec.seq;
        if is_dna(seq) {
            let rvds = dna_to_rvds(seq);
            tales.push((id, rvds));
        } else {
            let seq_str = String::from_utf8_lossy(seq);
            let rvds = parse_rvd_sequence(&seq_str);
            tales.push((id, rvds));
        }
    }

    // Set up output target
    let mut writer: Box<dyn Write> = match args.output {
        Some(ref path) => {
            let file = File::create(path)
                .with_context(|| format!("Failed to create output file {:?}", path))?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(std::io::stdout()),
    };

    // First scan: Internal Size/Frameshifts
    // Compares all pairs (i, j) where i < j
    for i in 0..tales.len() {
        for j in (i + 1)..tales.len() {
            let (id1, rvds1) = &tales[i];
            let (id2, rvds2) = &tales[j];

            if rvds1.len() != rvds2.len() && rvds1.len() >= 4 && rvds2.len() >= 4 {
                let prefix_match = rvds1[0..4] == rvds2[0..4];
                let suffix_match = rvds1[rvds1.len() - 4..] == rvds2[rvds2.len() - 4..];

                if prefix_match && suffix_match {
                    writeln!(writer, "{}", id1)?;
                    writeln!(writer, "{}", rvds1.join("-"))?;
                    writeln!(writer, "{}", id2)?;
                    writeln!(writer, "{}", rvds2.join("-"))?;
                    writeln!(writer, "+++++++++++++++++++++++++++++++++++++++++++++++++")?;
                }
            }
        }
    }

    writeln!(writer, "#######################################################")?;

    // Second scan: Truncations
    // Compares all pairs (i, j) where i < j
    for i in 0..tales.len() {
        for j in (i + 1)..tales.len() {
            let (id1, rvds1) = &tales[i];
            let (id2, rvds2) = &tales[j];

            let (long_id, long_rvds, short_id, short_rvds) = if rvds1.len() > rvds2.len() {
                (id1, rvds1, id2, rvds2)
            } else if rvds1.len() < rvds2.len() {
                (id2, rvds2, id1, rvds1)
            } else {
                continue;
            };

            if short_rvds.is_empty() {
                continue;
            }

            let prefix_match = long_rvds[0..short_rvds.len()] == *short_rvds;
            let suffix_match = long_rvds[long_rvds.len() - short_rvds.len()..] == *short_rvds;

            if prefix_match || suffix_match {
                writeln!(writer, "{}", long_id)?;
                writeln!(writer, "{}", long_rvds.join("-"))?;
                writeln!(writer, "{}", short_id)?;
                writeln!(writer, "{}", short_rvds.join("-"))?;
                writeln!(writer, "+++++++++++++++++++++++++++++++++++++++++++++++++")?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}
