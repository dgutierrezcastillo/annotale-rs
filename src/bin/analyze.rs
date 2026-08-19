use anyhow::{Context, Result};
use bio::io::fasta;
use clap::Parser;
use rust_annotale::translate;
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
    // A highly simplified regex-like sliding window for TALE repeats.
    // TALE repeats are typically 34 amino acids long (102 bp).
    // The core conserved motif is LTP[D/E]QVVAIAS which is near the start.
    // We will just break the sequence into 102 bp chunks starting from the first "LTP" or similar motif.
    let mut repeats = Vec::new();
    let mut rvd_str = String::new();

    let aa_seq = translate(sequence);
    
    // Find first 'LTP'
    let mut start_idx = None;
    for i in 0..aa_seq.len().saturating_sub(3) {
        if aa_seq[i] == b'L' && aa_seq[i+1] == b'T' && aa_seq[i+2] == b'P' {
            start_idx = Some(i);
            break;
        }
    }

    if let Some(start) = start_idx {
        let mut curr = start * 3;
        while curr + 102 <= sequence.len() {
            let repeat_dna = &sequence[curr..curr+102];
            repeats.push(repeat_dna.to_vec());
            
            let repeat_aa = translate(repeat_dna);
            if repeat_aa.len() >= 13 {
                if !rvd_str.is_empty() { rvd_str.push('-'); }
                rvd_str.push(repeat_aa[12] as char);
                rvd_str.push(repeat_aa[13] as char);
            }
            curr += 102;
        }
    }

    (repeats, rvd_str)
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.outdir.exists() {
        std::fs::create_dir_all(&args.outdir)?;
    }

    let reader = fasta::Reader::from_file(&args.tale_sequences)?;
    
    let dna_path = args.outdir.join("TALE_DNA_parts.fasta");
    let mut dna_writer = BufWriter::new(File::create(dna_path)?);
    
    let rvd_path = args.outdir.join("TALE_RVDs.fasta");
    let mut rvd_writer = BufWriter::new(File::create(rvd_path)?);

    for record in reader.records() {
        let rec = record?;
        let (repeats, rvd_str) = extract_repeats_and_rvds(rec.seq());
        
        writeln!(rvd_writer, ">{}\n{}", rec.id(), rvd_str)?;
        
        for (i, repeat) in repeats.iter().enumerate() {
            writeln!(dna_writer, ">{}_repeat_{}\n{}", rec.id(), i + 1, String::from_utf8_lossy(repeat))?;
        }
    }

    println!("Analysis complete. Outputs written to {:?}", args.outdir);
    Ok(())
}
