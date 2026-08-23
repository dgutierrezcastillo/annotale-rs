use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use flate2::write::GzEncoder;
use flate2::Compression;
use hmmer_pure_rs::alphabet::{Alphabet, AlphabetType};
use hmmer_pure_rs::bg::Bg;
use hmmer_pure_rs::hmmfile;
use hmmer_pure_rs::profile::{profile_config, P7_LOCAL};
use hmmer_pure_rs::sequence::Sequence as HmmSequence;
use hmmer_pure_rs::{Hmm, OProfile, Pipeline, Profile, TopHits};
use rayon::prelude::*;
use rust_annotale::{extract_consensus, open_sequence_reader, revcomp, SeqRecord};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum FilterMode {
    Heuristic,
    Hmm,
    Auto,
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum Preset {
    PacbioHifi,
    PacbioClr,
    Ont,
    None,
}

#[derive(Parser, Debug)]
#[command(
    name = "filter-reads",
    author,
    version,
    about = "Filter long-read sequencing datasets (PacBio/ONT) for TALE repeats"
)]
struct Args {
    #[arg(short, long, help = "Input FASTQ/FASTA file (.gz supported)")]
    input: String,

    #[arg(long, help = "Path to repeats.hmm file")]
    hmm: String,

    #[arg(long, help = "Output path for TALE-containing reads (.gz supported)")]
    output_repeats: String,

    #[arg(long, help = "Output path for non-TALE reads (.gz supported)")]
    output_norepeats: String,

    #[arg(
        short,
        long,
        value_enum,
        default_value_t = FilterMode::Auto,
        help = "Filtering mode"
    )]
    mode: FilterMode,

    #[arg(
        short,
        long,
        value_enum,
        default_value_t = Preset::None,
        help = "Technology preset"
    )]
    preset: Preset,

    #[arg(
        long,
        help = "Minimum number of matching 10-mer fragments (default: HMM length / 30)"
    )]
    min_parts: Option<usize>,

    #[arg(
        short,
        long,
        help = "Score threshold in bits (default: ~26.8 bits for repeats.hmm)"
    )]
    threshold: Option<f32>,

    #[arg(
        long,
        default_value_t = 10000,
        help = "Batch size for streaming reads"
    )]
    batch_size: usize,
}

struct SeqWriter<W: Write> {
    writer: W,
}

impl<W: Write> SeqWriter<W> {
    fn new(writer: W) -> Self {
        Self { writer }
    }

    fn write_record(&mut self, record: &SeqRecord) -> io::Result<()> {
        if let Some(qual) = &record.qual {
            writeln!(self.writer, "@{}", record.id)?;
            writeln!(self.writer, "{}", String::from_utf8_lossy(&record.seq))?;
            writeln!(self.writer, "+")?;
            writeln!(self.writer, "{}", String::from_utf8_lossy(qual))?;
        } else {
            writeln!(self.writer, ">{}", record.id)?;
            writeln!(self.writer, "{}", String::from_utf8_lossy(&record.seq))?;
        }
        Ok(())
    }
}

struct HmmScorer {
    repeats_hmm: Hmm,
    abc: Alphabet,
    bg: Bg,
}

impl HmmScorer {
    fn new(hmm_path: &Path) -> Result<Self> {
        let hmms = hmmfile::read_hmm_file(hmm_path)
            .with_context(|| format!("Failed to read HMM file {}", hmm_path.display()))?;

        let abc = Alphabet::new(AlphabetType::Dna);
        let bg = Bg::new(&abc);

        Ok(Self {
            repeats_hmm: hmms[0].clone(),
            abc,
            bg,
        })
    }

    fn score_window(&self, seq: &[u8]) -> f32 {
        let l_val = seq.len();
        let mut gm = Profile::new(self.repeats_hmm.m, &self.abc);
        profile_config(&self.repeats_hmm, &self.bg, &mut gm, l_val as i32, P7_LOCAL);
        let mut om = OProfile::convert(&gm);

        let mut pli = Pipeline::new();
        pli.new_model(&gm);
        let mut th = TopHits::new();

        let dsq = self.abc.digitize(seq);
        let sq = HmmSequence {
            name: "window".to_string(),
            acc: String::new(),
            desc: String::new(),
            taxid: -1,
            dsq,
            n: l_val,
            l: l_val,
        };

        pli.run(&mut gm, &mut om, &self.bg, &self.repeats_hmm, &sq, &mut th);

        let mut max_score = f32::NEG_INFINITY;
        for hit in &th.hits {
            for domain in &hit.dcl {
                if domain.bitscore > max_score {
                    max_score = domain.bitscore;
                }
            }
        }
        max_score
    }

    fn score_full_read(&self, seq: &[u8], threshold: f32) -> bool {
        let l_val = seq.len();
        if l_val < 100 {
            return false;
        }
        let mut gm = Profile::new(self.repeats_hmm.m, &self.abc);
        profile_config(&self.repeats_hmm, &self.bg, &mut gm, l_val as i32, P7_LOCAL);
        let mut om = OProfile::convert(&gm);

        let mut pli = Pipeline::new();
        pli.new_model(&gm);
        let mut th = TopHits::new();

        let dsq = self.abc.digitize(seq);
        let sq = HmmSequence {
            name: "read".to_string(),
            acc: String::new(),
            desc: String::new(),
            taxid: -1,
            dsq,
            n: l_val,
            l: l_val,
        };

        pli.run(&mut gm, &mut om, &self.bg, &self.repeats_hmm, &sq, &mut th);

        for hit in &th.hits {
            for domain in &hit.dcl {
                if domain.bitscore > threshold {
                    return true;
                }
            }
        }
        false
    }
}


fn open_writer(path: &Path) -> Result<BufWriter<Box<dyn Write>>> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let writer: Box<dyn Write> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };
    Ok(BufWriter::new(writer))
}

/// True when `haystack` contains `needle` as a contiguous substring (byte-exact).
fn contains_fragment(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

fn is_known_fragment(parts: &[Vec<u8>], fragment: &[u8]) -> bool {
    parts.iter().any(|p| p.as_slice() == fragment)
}

fn find_repeats_heuristic(
    seq: &[u8],
    scorer: &HmmScorer,
    parts: &[Vec<u8>],
    frag: usize,
    consensus_len: usize,
    min_parts: usize,
    threshold_bits: f32,
) -> bool {
    // Normalize once per call instead of allocating uppercased Strings on
    // every window slide.
    let seq = seq.to_ascii_uppercase();
    let num_lay = (consensus_len as f64 * 1.1).round() as usize;
    let w = num_lay;

    if seq.len() < w {
        return false;
    }

    let mut max_val = f32::NEG_INFINITY;
    let mut num = -1;

    for j in 0..=(seq.len() - w) {
        let sub = &seq[j..j + w];

        if num == -1 {
            num = 0;
            for part in parts {
                if contains_fragment(sub, part) {
                    num += 1;
                }
            }
        }

        if num > min_parts as i32 {
            let score = scorer.score_window(sub);
            if score > max_val {
                max_val = score;
                if max_val > threshold_bits {
                    return true;
                }
            }
        }

        if j < seq.len() - w {
            let leaving = &seq[j..j + frag];
            let entering = &seq[j + w - frag + 1..j + w + 1];

            if is_known_fragment(parts, leaving) {
                num -= 1;
            }
            if is_known_fragment(parts, entering) {
                num += 1;
            }
        }
    }

    max_val > threshold_bits
}

fn main() -> Result<()> {
    let args = Args::parse();
    let start_time = Instant::now();

    println!("Initializing HMM scorer from {}...", args.hmm);
    let scorer = HmmScorer::new(Path::new(&args.hmm))?;
    let consensus = extract_consensus(Path::new(&args.hmm))?;
    let consensus_len = consensus.len();

    // Setup 10-mer fragments (uppercased, deduplicated)
    let frag = 10;
    let mut parts: Vec<Vec<u8>> = Vec::new();
    for i in 0..(consensus_len / frag) {
        let fragment = consensus[i * frag..(i + 1) * frag]
            .to_ascii_uppercase()
            .into_bytes();
        if !parts.contains(&fragment) {
            parts.push(fragment);
        }
    }

    // Determine thresholds and modes
    let mut filter_mode = args.mode;
    let mut min_parts = args.min_parts.unwrap_or(parts.len() / 3);
    let mut threshold_bits = args.threshold.unwrap_or(26.8); // Equivalent of standard Java threshold

    match args.preset {
        Preset::PacbioHifi => {
            println!("Applying preset: PacBio HiFi (high speed, accurate reads)");
            filter_mode = FilterMode::Heuristic;
            min_parts = 3;
            threshold_bits = 26.8;
        }
        Preset::PacbioClr => {
            println!("Applying preset: PacBio CLR (robust full HMM scanning)");
            filter_mode = FilterMode::Hmm;
            threshold_bits = 12.0; // Lower bitscore threshold for higher error reads
        }
        Preset::Ont => {
            println!("Applying preset: ONT (robust full HMM scanning)");
            filter_mode = FilterMode::Hmm;
            threshold_bits = 10.0; // Lower threshold to accommodate high indel rates
        }
        Preset::None => {
            if filter_mode == FilterMode::Auto {
                filter_mode = FilterMode::Heuristic;
            }
        }
    }

    println!("Configuration:");
    println!("  Mode: {:?}", filter_mode);
    println!("  HMM length: {} bp", consensus_len);
    println!("  Threshold: {:.1} bits", threshold_bits);
    if filter_mode == FilterMode::Heuristic {
        println!("  Heuristic fragments: {} (min matching: {})", parts.len(), min_parts);
    }

    println!("Opening input file: {}...", args.input);
    let mut reader = open_sequence_reader(&args.input)?;

    println!("Opening output files:\n  Repeats: {}\n  No-repeats: {}", args.output_repeats, args.output_norepeats);
    let mut writer_repeats = SeqWriter::new(open_writer(Path::new(&args.output_repeats))?);
    let mut writer_norepeats = SeqWriter::new(open_writer(Path::new(&args.output_norepeats))?);

    let mut total_reads = 0;
    let mut repeat_reads = 0;
    let mut norepeat_reads = 0;

    let mut batch = Vec::with_capacity(args.batch_size);

    println!("Filtering reads in batches of {}...", args.batch_size);

    loop {
        batch.clear();
        for _ in 0..args.batch_size {
            if let Some(record) = reader.next_record()? {
                batch.push(record);
            } else {
                break;
            }
        }

        if batch.is_empty() {
            break;
        }

        total_reads += batch.len();

        // Process batch in parallel using Rayon
        let results: Vec<(SeqRecord, bool)> = batch
            .par_iter()
            .map(|record| {
                let seq = &record.seq;
                let rev = revcomp(seq);
                
                let matches = match filter_mode {
                    FilterMode::Heuristic => {
                        find_repeats_heuristic(seq, &scorer, &parts, frag, consensus_len, min_parts, threshold_bits)
                            || find_repeats_heuristic(&rev, &scorer, &parts, frag, consensus_len, min_parts, threshold_bits)
                    }
                    FilterMode::Hmm | FilterMode::Auto => {
                        scorer.score_full_read(seq, threshold_bits)
                            || scorer.score_full_read(&rev, threshold_bits)
                    }
                };

                (record.clone(), matches)
            })
            .collect();

        // Write results to appropriate files on main thread
        for (record, has_repeats) in results {
            if has_repeats {
                writer_repeats.write_record(&record)?;
                repeat_reads += 1;
            } else {
                writer_norepeats.write_record(&record)?;
                norepeat_reads += 1;
            }
        }

        print!("\rProcessed {} reads... ({} repeats, {} non-repeats)", total_reads, repeat_reads, norepeat_reads);
        io::stdout().flush()?;
    }

    println!("\n\nFiltering complete!");
    println!("Results:");
    println!("  Total reads processed: {}", total_reads);
    println!("  TALE-containing reads: {} ({:.2}%)", repeat_reads, (repeat_reads as f64 / total_reads as f64) * 100.0);
    println!("  Non-TALE reads:        {} ({:.2}%)", norepeat_reads, (norepeat_reads as f64 / total_reads as f64) * 100.0);
    println!("  Total time elapsed:    {:.2?}", start_time.elapsed());

    Ok(())
}
