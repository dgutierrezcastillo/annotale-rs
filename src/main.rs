use anyhow::{Context, Result};
use bio::alphabets::dna::revcomp;
use bio::seq_analysis::orf::{Finder, Orf};
use clap::Parser;
use hmmer_pure_rs::alphabet::{Alphabet, AlphabetType};
use hmmer_pure_rs::bg::Bg;
use hmmer_pure_rs::hmmfile;
use hmmer_pure_rs::profile::{profile_config, P7_LOCAL};
use hmmer_pure_rs::sequence::Sequence as HmmSequence;
use hmmer_pure_rs::{Hmm, OProfile, Pipeline, Profile, TopHits};
use rayon::prelude::*;
use rust_annotale::{open_sequence_reader, SeqRecord, TALERegion, Translator};
use std::path::Path;

const CLUSTER_MAX_GAP_BP: usize = 500;
const BOUNDARY_BUFFER_BP: usize = 1200;
const MIN_ORF_LEN_BP: usize = 300;
const PSEUDO_ORF_THRESHOLD_BP: usize = 900;
const RVD_AA_OFFSET_BP: usize = 36;

#[derive(Parser, Debug)]
#[command(author, version, about = "A Rust implementation of AnnoTALE")]
struct Args {
    /// Input sequence file: FASTA or FASTQ (.gz supported)
    #[arg(short, long, alias = "fasta")]
    input: String,

    #[arg(long)]
    hmm_dir: String,

    #[arg(short, long, default_value_t = 10.0)]
    threshold: f32,

    /// Metagenomic mode: streams large multi-contig/read files in batches, suppresses empty logs, and skips very short fragments
    #[arg(short = 'm', long = "metagenome")]
    metagenome: bool,

    /// Minimum contig/read length in bp to scan (default: 200 bp)
    #[arg(long, default_value_t = 200)]
    min_length: usize,

    /// Batch size for parallel processing in streaming/metagenomic mode
    #[arg(long, default_value_t = 1000)]
    batch_size: usize,
}

struct TALEFinder {
    repeats_hmm: Hmm,
    abc: Alphabet,
    bg: Bg,
    translator: Translator,
}

impl TALEFinder {
    fn new(hmm_dir: &str) -> Result<Self> {
        let repeats_path = Path::new(hmm_dir).join("repeats.hmm");
        let repeats_hmms = hmmfile::read_hmm_file(&repeats_path)
            .with_context(|| format!("Failed to read {}", repeats_path.display()))?;

        let abc = Alphabet::new(AlphabetType::Dna);
        let bg = Bg::new(&abc);

        Ok(Self {
            repeats_hmm: repeats_hmms[0].clone(),
            abc,
            bg,
            translator: Translator::new(),
        })
    }

    fn scan_sequence(&self, record_id: &str, sequence: &[u8]) -> Vec<TALERegion> {
        let mut results = Vec::new();
        results.extend(self.process_strand(record_id, sequence, '+'));

        let rev_seq = revcomp(sequence);
        results.extend(self.process_strand(record_id, &rev_seq, '-'));

        results
    }

    fn process_strand(&self, id: &str, sequence: &[u8], strand: char) -> Vec<TALERegion> {
        let mut raw_matches = self.run_hmmer_raw(id, sequence);
        if raw_matches.is_empty() {
            return Vec::new();
        }

        raw_matches.sort_by_key(|m| m.0);
        let clusters = group_matches_into_clusters(&raw_matches);

        let mut final_tales = Vec::new();
        for (domains, c_score) in clusters {
            let c_start = domains[0].0;
            let c_end = domains.last().unwrap().1;

            let search_start = c_start.saturating_sub(BOUNDARY_BUFFER_BP);
            let search_end = std::cmp::min(sequence.len(), c_end + BOUNDARY_BUFFER_BP);

            let (cds_rel_start, cds_rel_end, is_pseudo) =
                self.refine_cds(sequence, search_start, search_end);
            let final_cds_start = search_start + cds_rel_start;
            let final_cds_end = search_start + cds_rel_end;

            let rvds = if !is_pseudo && (final_cds_end - final_cds_start) > 100 {
                self.extract_rvds(sequence, final_cds_start, &domains)
            } else {
                "N/A".to_string()
            };

            let (actual_start, actual_end) = if strand == '+' {
                (final_cds_start, final_cds_end)
            } else {
                let l = sequence.len();
                (l - final_cds_end, l - final_cds_start)
            };

            final_tales.push(TALERegion {
                strand,
                start: actual_start,
                end: actual_end,
                score: c_score,
                is_pseudo,
                rvds,
            });
        }

        final_tales
    }

    fn run_hmmer_raw(&self, id: &str, sequence: &[u8]) -> Vec<(usize, usize, f32)> {
        let mut matches = Vec::new();
        let l_val = sequence.len();
        let mut gm = Profile::new(self.repeats_hmm.m, &self.abc);
        profile_config(&self.repeats_hmm, &self.bg, &mut gm, l_val as i32, P7_LOCAL);
        let mut om = OProfile::convert(&gm);

        let mut pli = Pipeline::new();
        pli.new_model(&gm);
        let mut th = TopHits::new();

        let dsq = self.abc.digitize(sequence);
        let sq = HmmSequence {
            name: id.to_string(),
            acc: String::new(),
            desc: String::new(),
            dsq,
            n: l_val,
            l: l_val,
        };

        pli.run(&mut gm, &mut om, &self.bg, &self.repeats_hmm, &sq, &mut th);

        for hit in &th.hits {
            for domain in &hit.dcl {
                if domain.bitscore > 10.0 {
                    matches.push((domain.iali as usize, domain.jali as usize, domain.bitscore));
                }
            }
        }
        matches
    }

    fn refine_cds(&self, sequence: &[u8], start: usize, end: usize) -> (usize, usize, bool) {
        let target_seq = &sequence[start..end];
        let start_codons = vec![b"ATG"];
        let stop_codons = vec![b"TGA", b"TAG", b"TAA"];
        let finder = Finder::new(start_codons, stop_codons, MIN_ORF_LEN_BP);

        let mut max_len = 0;
        let mut best_orf: Option<Orf> = None;

        for orf in finder.find_all(target_seq) {
            let len = orf.end - orf.start;
            if len > max_len {
                max_len = len;
                best_orf = Some(orf);
            }
        }

        if let Some(orf) = best_orf {
            let is_pseudo = (orf.end - orf.start) < PSEUDO_ORF_THRESHOLD_BP;
            (orf.start, orf.end, is_pseudo)
        } else {
            (0, 0, true)
        }
    }

    fn extract_rvds(
        &self,
        sequence: &[u8],
        cds_start: usize,
        domains: &[(usize, usize, f32)],
    ) -> String {
        let mut rvd_str = String::new();

        for m in domains {
            let domain_start = m.0;
            if domain_start < cds_start {
                continue;
            }

            let offset = domain_start - cds_start;
            let frame_shift = offset % 3;
            let aligned_start = domain_start - frame_shift;

            let rvd_dna_start = aligned_start + RVD_AA_OFFSET_BP;
            let rvd_dna_end = rvd_dna_start + 6;

            if rvd_dna_end <= sequence.len() {
                let rvd_dna = &sequence[rvd_dna_start..rvd_dna_end];
                let rvd_aa = self.translator.translate(rvd_dna);
                if rvd_aa.len() >= 2 {
                    if !rvd_str.is_empty() {
                        rvd_str.push('-');
                    }
                    rvd_str.push(rvd_aa[0] as char);
                    rvd_str.push(rvd_aa[1] as char);
                }
            }
        }

        if rvd_str.is_empty() {
            "N/A".to_string()
        } else {
            rvd_str
        }
    }
}

fn group_matches_into_clusters(
    raw_matches: &[(usize, usize, f32)],
) -> Vec<(Vec<(usize, usize, f32)>, f32)> {
    let mut clusters = Vec::new();
    if raw_matches.is_empty() {
        return clusters;
    }

    let mut current_cluster = vec![raw_matches[0]];
    let mut c_end = raw_matches[0].1;
    let mut c_score = raw_matches[0].2;

    for &match_item in &raw_matches[1..] {
        if match_item.0 < c_end + CLUSTER_MAX_GAP_BP {
            c_end = std::cmp::max(c_end, match_item.1);
            c_score += match_item.2;
            current_cluster.push(match_item);
        } else {
            clusters.push((current_cluster.clone(), c_score));
            current_cluster = vec![match_item];
            c_end = match_item.1;
            c_score = match_item.2;
        }
    }
    clusters.push((current_cluster, c_score));
    clusters
}

fn main() -> Result<()> {
    let args = Args::parse();

    println!("Initializing TALEFinder with HMMs from {}...", args.hmm_dir);
    let finder = TALEFinder::new(&args.hmm_dir)?;

    println!(
        "Scanning {} for TALE effectors{}...",
        args.input,
        if args.metagenome { " (metagenomic mode)" } else { "" }
    );

    let mut reader = open_sequence_reader(&args.input)?;
    let mut total_scanned = 0;
    let mut total_tales = 0;
    let mut batch = Vec::with_capacity(args.batch_size);

    while let Some(record) = reader.next_record()? {
        if record.seq.len() < args.min_length {
            continue;
        }
        batch.push(record);

        if batch.len() >= args.batch_size {
            total_scanned += batch.len();
            total_tales += process_batch(&batch, &finder, args.metagenome);
            batch.clear();
        }
    }

    if !batch.is_empty() {
        total_scanned += batch.len();
        total_tales += process_batch(&batch, &finder, args.metagenome);
    }

    println!(
        "\nScan complete: {} contigs/sequences scanned, {} TALE effectors found.",
        total_scanned, total_tales
    );

    Ok(())
}

fn process_batch(
    batch: &[SeqRecord],
    finder: &TALEFinder,
    metagenome: bool,
) -> usize {
    let results: Vec<(String, Vec<TALERegion>)> = batch
        .par_iter()
        .map(|record| {
            let id = &record.id;
            let seq = &record.seq;
            let mut matches = finder.scan_sequence(id, seq);
            matches.sort_by_key(|m| m.start);
            (id.clone(), matches)
        })
        .collect();

    let mut tales_found = 0;
    for (id, matches) in results {
        if !matches.is_empty() {
            tales_found += matches.len();
            println!("\nFound {} potential TALE effectors in {}", matches.len(), id);
            print_table(&matches);
        } else if !metagenome {
            println!("No TALE effectors found in {}", id);
        }
    }
    tales_found
}

fn print_table(matches: &[TALERegion]) {
    println!(
        "{:<5} | {:<8} | {:<2} | {:<10} | {:<10} | {:<8} | {:<10}",
        "No.", "Type", "St", "Start", "End", "Score", "RVDs"
    );
    println!(
        "{:-<5}-|-{:-<8}-|-{:-<2}-|-{:-<10}-|-{:-<10}-|-{:-<8}-|-{:-<10}",
        "", "", "", "", "", "", ""
    );

    for (i, region) in matches.iter().enumerate() {
        let status = if region.is_pseudo { "PSEUDO" } else { "CDS" };
        println!(
            "{:5} | {:<8} | {:<2} | {:10} | {:10} | {:8.1} | {}",
            i + 1,
            status,
            region.strand,
            region.start,
            region.end,
            region.score,
            region.rvds
        );
    }
}
