use anyhow::Result;
use clap::Parser;
use rust_annotale::{dna_to_rvds, find_first_ltp, open_sequence_reader, translate};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about = "TALE Analysis", long_about = None)]
struct Args {
    /// The DNA sequences of the TALEs
    #[arg(short = 't', long = "tale-sequences", required = true)]
    tale_sequences: PathBuf,

    /// The output directory
    #[arg(long = "outdir", default_value = ".")]
    outdir: PathBuf,
}

fn extract_repeats_and_rvds(sequence: &[u8]) -> (Vec<Vec<u8>>, String) {
    // Repeat parts are 102-bp chunks from the first LTP motif; the RVD
    // string comes from the shared walker so fixes land everywhere at once.
    let mut repeats = Vec::new();
    if let Some(start) = find_first_ltp(&translate(sequence)) {
        let mut curr = start * 3;
        while curr + 102 <= sequence.len() {
            repeats.push(sequence[curr..curr + 102].to_vec());
            curr += 102;
        }
    }
    let rvd_str = dna_to_rvds(sequence).join("-");
    (repeats, rvd_str)
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.outdir.exists() {
        std::fs::create_dir_all(&args.outdir)?;
    }

    let mut reader = open_sequence_reader(&args.tale_sequences)?;

    let dna_path = args.outdir.join("TALE_DNA_parts.fasta");
    let mut dna_writer = BufWriter::new(File::create(dna_path)?);

    let rvd_path = args.outdir.join("TALE_RVDs.fasta");
    let mut rvd_writer = BufWriter::new(File::create(rvd_path)?);

    while let Some(rec) = reader.next_record()? {
        let (repeats, rvd_str) = extract_repeats_and_rvds(&rec.seq);

        writeln!(rvd_writer, ">{}\n{}", rec.id, rvd_str)?;

        for (i, repeat) in repeats.iter().enumerate() {
            writeln!(dna_writer, ">{}_repeat_{}\n{}", rec.id, i + 1, String::from_utf8_lossy(repeat))?;
        }
    }

    println!("Analysis complete. Outputs written to {:?}", args.outdir);
    Ok(())
}
