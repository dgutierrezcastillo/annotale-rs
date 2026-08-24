use anyhow::{Context, Result};

use clap::Parser;
use hmmer_pure_rs::alphabet::{Alphabet, AlphabetType};
use hmmer_pure_rs::bg::Bg;
use hmmer_pure_rs::hmmfile;
use hmmer_pure_rs::profile::{profile_config, P7_LOCAL};
use hmmer_pure_rs::sequence::Sequence as HmmSequence;
use hmmer_pure_rs::{Hmm, OProfile, Pipeline, Profile, TopHits};
use rayon::prelude::*;
use rust_annotale::{open_sequence_reader, revcomp, SeqRecord, TALERegion};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

const CLUSTER_MAX_GAP_BP: usize = 500;
const PREFILTER_FRAG_BP: usize = 10;
const TERM_SCAN_STEP: i64 = 5;
/// Minimum bitscore for a repeats.hmm domain to count as one repeat when
/// extracting RVDs. Override with env ANNOTALE_EXTRACT_BITS for tuning.
fn extract_min_domain_bits() -> f32 {
    static VAL: std::sync::OnceLock<f32> = std::sync::OnceLock::new();
    *VAL.get_or_init(|| {
        std::env::var("ANNOTALE_EXTRACT_BITS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(EXTRACT_MIN_DOMAIN_BITS)
    })
}

const EXTRACT_MIN_DOMAIN_BITS: f32 = 10.0;

#[derive(Parser, Debug)]
#[command(author, version, about = "A Rust implementation of AnnoTALE")]
struct Args {
    /// Input sequence file: FASTA or FASTQ (.gz supported)
    #[arg(short, long, alias = "fasta")]
    input: String,

    #[arg(long)]
    hmm_dir: String,

    /// Peak bitscore cutoff override; 0 uses the Java default
    /// (consensus_len * ln(1.3), converted to bits)
    #[arg(short, long, default_value_t = 0.0)]
    threshold: f32,

    /// Write results as a TSV table (contig, strand, start, end, score,
    /// pseudo, RVDs) to this file instead of formatted tables on stdout
    #[arg(short = 'O', long)]
    output: Option<String>,

    /// Metagenomic mode: streams large multi-contig/read files in batches, suppresses empty logs, and skips very short fragments
    #[arg(short = 'm', long = "metagenome")]
    metagenome: bool,

    /// Minimum contig/read length in bp to scan (default: 200 bp)
    #[arg(long, default_value_t = 200)]
    min_length: usize,

    /// Batch size for parallel processing in streaming/metagenomic mode
    #[arg(long, default_value_t = 1000)]
    batch_size: usize,

    /// Disable fast k-mer heuristic pre-filtering
    #[arg(long)]
    no_kmer_filter: bool,

    /// Minimum matching k-mers required to trigger full HMM profile scan
    #[arg(long, default_value_t = 2)]
    min_kmers: usize,
}

/// Scores fixed-length windows against one profile HMM.
///
/// Profiles are built once and reconfigured per window inside
/// `Pipeline::run`, so thousands of windows share one allocation.
struct WindowScorer {
    hmm: Hmm,
    abc: Alphabet,
    bg: Bg,
    gm: Profile,
    om: OProfile,
    pli: Pipeline,
}

impl WindowScorer {
    fn new(hmm: Hmm, window_len: usize) -> Self {
        let abc = Alphabet::new(AlphabetType::Dna);
        let bg = Bg::new(&abc);
        let mut gm = Profile::new(hmm.m, &abc);
        profile_config(&hmm, &bg, &mut gm, window_len as i32, P7_LOCAL);
        let om = OProfile::convert(&gm);
        let mut pli = Pipeline::new();
        pli.new_model(&gm);
        Self {
            hmm,
            abc,
            bg,
            gm,
            om,
            pli,
        }
    }

    /// Best domain alignment of `seq` against the model.
    /// Returns (bitscore, ali_start_0based, ali_end_exclusive); NEG_INFINITY when nothing aligns.
    fn score_window(&mut self, name: &str, seq: &[u8]) -> (f32, usize, usize) {
        let l_val = seq.len();
        if l_val == 0 {
            return (f32::NEG_INFINITY, 0, 0);
        }
        let dsq = self.abc.digitize(seq);
        let sq = HmmSequence {
            name: name.to_string(),
            acc: String::new(),
            desc: String::new(),
            taxid: -1,
            dsq,
            n: l_val,
            l: l_val,
        };
        let mut th = TopHits::new();
        self.pli
            .run(&mut self.gm, &mut self.om, &self.bg, &self.hmm, &sq, &mut th);

        let mut best = (f32::NEG_INFINITY, 0usize, 0usize);
        for hit in &th.hits {
            for domain in &hit.dcl {
                if domain.bitscore > best.0 {
                    best = (
                        domain.bitscore,
                        domain.iali.saturating_sub(1) as usize,
                        domain.jali as usize,
                    );
                }
            }
        }
        best
    }
}

struct TALEFinder {
    consensus_len: usize,
    window_len: usize,
    /// Java peak threshold: consensus_len * ln(1.3) nats, stored in bits.
    peak_threshold_bits: f32,
    /// Consensus 10-bp fragments used by the rolling window prefilter and k-mer gate.
    prefilter_parts: Vec<Vec<u8>>,
    min_kmers: usize,
    use_kmer_filter: bool,
    repeats: Mutex<WindowScorer>,
    starts: Option<Mutex<WindowScorer>>,
    ends: Option<Mutex<WindowScorer>>,
}

impl TALEFinder {
    fn new(
        hmm_dir: &str,
        use_kmer_filter: bool,
        min_kmers: usize,
        threshold_override: f32,
    ) -> Result<Self> {
        let load = |name: &str| -> Result<Hmm> {
            let path = Path::new(hmm_dir).join(name);
            let hmms = hmmfile::read_hmm_file(&path)
                .with_context(|| format!("Failed to read {}", path.display()))?;
            Ok(hmms[0].clone())
        };

        let repeats_hmm = load("repeats.hmm")?;
        let starts_path = Path::new(hmm_dir).join("starts.hmm");
        let ends_path = Path::new(hmm_dir).join("ends.hmm");
        let starts_hmm = if starts_path.exists() {
            Some(load("starts.hmm")?)
        } else {
            eprintln!(
                "warning: {} missing; terminus refinement disabled",
                starts_path.display()
            );
            None
        };
        let ends_hmm = if ends_path.exists() {
            Some(load("ends.hmm")?)
        } else {
            eprintln!(
                "warning: {} missing; terminus refinement disabled",
                ends_path.display()
            );
            None
        };

        let repeats_path = Path::new(hmm_dir).join("repeats.hmm");
        let consensus = rust_annotale::extract_consensus(&repeats_path)?.to_ascii_uppercase();
        let consensus_len = consensus.len();
        if consensus_len == 0 {
            anyhow::bail!("repeats.hmm has no consensus");
        }
        let window_len = (consensus_len as f64 * 1.1).round() as usize;

        let peak_threshold_bits = if threshold_override > 0.0 {
            threshold_override
        } else {
            consensus_len as f32 * 1.3f32.ln() / 2f32.ln()
        };

        Ok(Self {
            consensus_len,
            window_len,
            peak_threshold_bits,
            prefilter_parts: rust_annotale::extract_kmers(&consensus, PREFILTER_FRAG_BP),
            min_kmers,
            use_kmer_filter,
            repeats: Mutex::new(WindowScorer::new(repeats_hmm, window_len)),
            starts: starts_hmm.map(|h| Mutex::new(WindowScorer::new(h, window_len))),
            ends: ends_hmm.map(|h| Mutex::new(WindowScorer::new(h, window_len))),
        })
    }

    #[inline]
    fn passes_kmer_filter(&self, sequence: &[u8]) -> bool {
        if !self.use_kmer_filter || self.prefilter_parts.is_empty() {
            return true;
        }
        let mut matches = 0;
        for kmer in &self.prefilter_parts {
            if sequence.windows(kmer.len()).any(|w| w.eq_ignore_ascii_case(kmer)) {
                matches += 1;
                if matches >= self.min_kmers {
                    return true;
                }
            }
        }
        false
    }

    fn scan_sequence(&self, record_id: &str, sequence: &[u8]) -> Vec<TALERegion> {
        let mut results = Vec::new();

        if self.passes_kmer_filter(sequence) {
            results.extend(self.process_strand(record_id, sequence, '+'));
        }

        let rev_seq = revcomp(sequence);
        if self.passes_kmer_filter(&rev_seq) {
            results.extend(self.process_strand(record_id, &rev_seq, '-'));
        }

        dedup_nested(results)
    }

    /// Java NHMMer.findRepeats port: rolling consensus-fragment prefilter,
    /// HMM-scored windows, greedy peak picking with +-window/2 suppression,
    /// position-sorted clustering with a 500 bp gap.
    fn find_repeat_regions(&self, id: &str, sequence: &[u8]) -> Vec<(usize, usize, f32)> {
        let w = self.window_len;
        let up = sequence.to_ascii_uppercase();
        if up.len() < w {
            return Vec::new();
        }
        let n_windows = up.len() - w + 1;
        let parts = &self.prefilter_parts;
        let n_parts = parts.len();
        if n_parts == 0 {
            return Vec::new();
        }

        let mut vals = vec![0f32; n_windows];
        {
            let mut scorer = self.repeats.lock().unwrap_or_else(|e| e.into_inner());
            let mut num: i32 = -1;
            for (j, v) in vals.iter_mut().enumerate() {
                let sub = &up[j..j + w];
                if num == -1 {
                    num = 0;
                    for part in parts {
                        if contains_fragment(sub, part) {
                            num += 1;
                        }
                    }
                }
                if num as usize > n_parts / 2 {
                    *v = scorer.score_window(id, sub).0;
                }
                if j < up.len() - w {
                    let leaving = &up[j..j + PREFILTER_FRAG_BP];
                    let entering = &up[j + w - PREFILTER_FRAG_BP + 1..j + w + 1];
                    if parts.iter().any(|p| p.as_slice() == leaving) {
                        num -= 1;
                    }
                    if parts.iter().any(|p| p.as_slice() == entering) {
                        num += 1;
                    }
                }
            }
        }

        let peaks = pick_peaks(&mut vals, w / 2, self.peak_threshold_bits);
        if std::env::var("ANNOTALE_DEBUG").is_ok() {
            let scored = vals.iter().filter(|&&v| v != 0.0).count();
            let max_v = vals.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
            eprintln!(
                "[debug] {id} len={} windows={n_windows} scored={scored} max_score={max_v:.1} t={:.1} peaks={}",
                sequence.len(),
                self.peak_threshold_bits,
                peaks.len(),
                id = id,
                n_windows = n_windows,
            );
        }
        cluster_peaks(&peaks, CLUSTER_MAX_GAP_BP, self.consensus_len)
    }

    /// Java NHMMer.getBestTerminus port: walk outward from the repeat cluster
    /// in steps of 5 bp, scoring windows against a flanking profile HMM;
    /// stop early after >10 consecutively decreasing scores; refine the
    /// cluster boundary to the best window's alignment extent.
    ///
    /// Deviation from Java: the refined coordinate comes from the HMMER
    /// domain alignment bounds instead of a Viterbi match-state run.
    fn get_best_terminus(
        &self,
        which: &Mutex<WindowScorer>,
        id: &str,
        sequence: &[u8],
        start: usize,
        end: usize,
        is_start: bool,
    ) -> Option<(usize, usize)> {
        let range = terminus_scan_range(
            is_start,
            start as i64,
            end as i64,
            sequence.len(),
            self.window_len,
            self.consensus_len,
        )?;
        let (from, to) = range;
        let mut scorer = which.lock().unwrap_or_else(|e| e.into_inner());

        let mut best: Option<(f32, i64, usize, usize)> = None;
        let mut prev_added = f32::INFINITY;
        let mut have_prev = false;
        let mut num_lower = 0;

        let mut i = from;
        while (is_start && i >= to) || (!is_start && i < to) {
            let idx = i as usize;
            let (score, ali_s, ali_e) = scorer.score_window(id, &sequence[idx..idx + self.window_len]);
            if have_prev && score < prev_added {
                num_lower += 1;
            } else {
                num_lower = 0;
            }
            if num_lower > 10 {
                break;
            }
            prev_added = score;
            have_prev = true;
            if score > best.map_or(f32::NEG_INFINITY, |b| b.0) {
                best = Some((score, i, ali_s, ali_e));
            }
            if is_start {
                i -= TERM_SCAN_STEP;
            } else {
                i += TERM_SCAN_STEP;
            }
        }

        let (_, pos, ali_s, ali_e) = best?;
        if is_start {
            Some((pos as usize + ali_s, 0))
        } else {
            Some((0, pos as usize + ali_e))
        }
    }

    fn process_strand(&self, id: &str, sequence: &[u8], strand: char) -> Vec<TALERegion> {
        let regions = self.find_repeat_regions(id, sequence);
        if regions.is_empty() {
            return Vec::new();
        }

        let mut final_tales = Vec::new();
        for (c_start, c_end, _peak_score) in regions {
            let mut term_start = c_start;
            let mut term_end = c_end;
            if let Some(scorer) = self.starts.as_ref() {
                if let Some((s, _)) = self.get_best_terminus(scorer, id, sequence, c_start, c_end, true) {
                    term_start = s;
                }
            }
            if let Some(scorer) = self.ends.as_ref() {
                if let Some((_, e)) = self.get_best_terminus(scorer, id, sequence, c_start, c_end, false) {
                    term_end = e;
                }
            }

            // Java refine(): longest ORF between stops within the refined span.
            let (cds_start, cds_end) = self.refine_cds(sequence, term_start, term_end.saturating_sub(1));
            if std::env::var("ANNOTALE_DEBUG").is_ok() {
                eprintln!(
                    "[debug] {id} {strand} cluster=[{c_start},{c_end}) term=[{term_start},{term_end}) cds=[{cds_start},{cds_end})"
                );
            }

            // Degenerate / no-ORF guard (P1): Java's refine always yields a
            // positive-length span; when ours cannot, drop the cluster.
            if cds_end <= cds_start {
                continue;
            }

            // Java pseudo-gene flag: the terminus-refined mRNA span exceeds
            // the CDS by more than one terminal domain.
            let mrna_len = term_end.max(cds_end) - term_start.min(cds_start);
            let cons_aa = self.consensus_len / 3;
            let is_pseudo =
                mrna_len.saturating_sub(cons_aa) > (cds_end - cds_start);

            let rvds = if cds_end - cds_start > 100 {
                self.extract_rvds_hmm(id, sequence, cds_start, cds_end)
            } else {
                "N/A".to_string()
            };

            let (actual_start, actual_end) = if strand == '+' {
                (cds_start, cds_end)
            } else {
                let l = sequence.len();
                (l - cds_end, l - cds_start)
            };

            final_tales.push(TALERegion {
                strand,
                start: actual_start,
                end: actual_end,
                score: (cds_end - cds_start) as f32,
                is_pseudo,
                rvds,
            });
        }

        final_tales
    }

    /// Per-repeat HMM-aligned RVD extraction.
    ///
    /// Scans the CDS with `repeats.hmm` (a 102-nt nucleotide profile HMM
    /// modelling one TALE repeat), takes each domain alignment block, and
    /// reads the RVD diresidue (codons 12/13) from that block alone. Unlike
    /// the fixed-stride walker both tools shipped with, this cannot drift
    /// after indels between repeats: every repeat is anchored by its own
    /// HMM alignment.
    fn extract_rvds_hmm(&self, id: &str, sequence: &[u8], cds_start: usize, cds_end: usize) -> String {
        let seq = &sequence[cds_start..cds_end];
        let mut scorer = self.repeats.lock().unwrap_or_else(|e| e.into_inner());

        let l_val = seq.len();
        let WindowScorer {
            hmm,
            abc,
            bg,
            gm,
            om,
            pli,
        } = &mut *scorer;
        let dsq = abc.digitize(seq);
        let sq = HmmSequence {
            name: id.to_string(),
            acc: String::new(),
            desc: String::new(),
            taxid: -1,
            dsq,
            n: l_val,
            l: l_val,
        };
        let mut th = TopHits::new();
        pli.run(gm, om, bg, hmm, &sq, &mut th);

        let mut blocks: Vec<(usize, usize, f32)> = Vec::new();
        for hit in &th.hits {
            for domain in &hit.dcl {
                if domain.bitscore > extract_min_domain_bits() {
                    blocks.push((
                        domain.iali.saturating_sub(1) as usize,
                        domain.jali as usize,
                        domain.bitscore,
                    ));
                }
            }
        }
        if blocks.is_empty() {
            return "N/A".to_string();
        }
        if std::env::var("ANNOTALE_DEBUG").is_ok() {
            eprintln!(
                "[debug] rvds {id}: {} blocks {:?} cds_len={l}",
                blocks.len(),
                &blocks[..blocks.len().min(40)],
                l = seq.len()
            );
        }
        blocks.sort_by_key(|b| b.0);

        let mut rvds: Vec<String> = Vec::new();
        let mut prev_end = 0usize;

        // Lock the CDS reading frame once by majority vote: the correct
        // phase is the one under which the most blocks open with a
        // translatable LTP anchor. HMMER bounds are arbitrary bp and do not
        // respect codons, but the CDS itself is codon-aligned, so one
        // global phase applies to every block.
        let mut votes = [0usize; 3];
        for &(ali_s, ali_e, _) in &blocks {
            for (phase, vote) in votes.iter_mut().enumerate() {
                let s = ali_s + ((phase + 3 - ali_s % 3) % 3);
                if ali_e > s && !rvds_from_block(&seq[s..ali_e]).is_empty() {
                    *vote += 1;
                }
            }
        }
        let global_phase = votes
            .iter()
            .enumerate()
            .max_by_key(|(_, v)| **v)
            .map(|(i, _)| i)
            .unwrap_or(0);

        for (ali_s, ali_e, _bits) in blocks.clone() {
            let mut s = ali_s.max(prev_end);
            let ext_floor = prev_end;
            if ali_e <= s || ali_e > seq.len() {
                continue;
            }
            prev_end = ali_e;
            s += (global_phase + 3 - s % 3) % 3;
            if ali_e <= s {
                continue;
            }
            let mut block_rvds = rvds_from_block(&seq[s..ali_e]);
            if block_rvds.is_empty() && s > ext_floor {
                // HMMER sometimes trims the repeat's leading anchor residue
                // off a domain boundary; retry with the start pulled back
                // toward the previous block (same reading frame).
                let mut s2 = ext_floor.max(s.saturating_sub(12));
                s2 += (global_phase + 3 - s2 % 3) % 3;
                if s2 < s {
                    block_rvds = rvds_from_block(&seq[s2..ali_e]);
                }
            }
            if std::env::var("ANNOTALE_DEBUG").is_ok() {
                eprintln!(
                    "[debug] rvdblock {id} [{s},{ali_e}) len={} -> {:?}",
                    ali_e - s,
                    block_rvds
                );
            }
            rvds.extend(block_rvds);
        }

        if rvds.is_empty() {
            "N/A".to_string()
        } else {
            rvds.join("-")
        }
    }

    /// Faithful port of Java NHMMer.refine(): longest stop-to-stop ORF over
    /// three frames inside the inclusive span (no ATG requirement), then
    /// snapped to an ATG — searched inside the ORF when it follows another
    /// stop, otherwise by walking codon-wise upstream out of the span — and
    /// extended downstream codon-wise until a stop codon / contig end.
    /// Returns absolute half-open CDS bounds.
    fn refine_cds(&self, sequence: &[u8], start: usize, end_incl: usize) -> (usize, usize) {
        const STOP: u8 = b'*';
        const ATG: &[u8] = b"ATG";
        if end_incl < start || end_incl >= sequence.len() || start >= sequence.len() {
            return (0, 0);
        }
        let span = &sequence[start..=end_incl];

        // Longest inter-stop segment across all three frames.
        let mut max_len_aa = 0usize;
        let mut max_frame = 0usize;
        let mut max_part = 0usize; // index of the winning segment within its frame

        for frame in 0..3usize {
            let n_aa = if frame + 3 <= span.len() { (span.len() - frame) / 3 } else { 0 };
            let mut seg_start = 0usize;
            let mut seg_idx = 0usize;
            for k in 0..=n_aa {
                let is_stop = k < n_aa
                    && rust_annotale::translate_codon(&span[frame + k * 3..frame + k * 3 + 3]) == STOP;
                if k == n_aa || is_stop {
                    let seg_len = k - seg_start;
                    if seg_len > max_len_aa {
                        max_len_aa = seg_len;
                        max_frame = frame;
                        max_part = seg_idx;
                    }
                    seg_start = k + 1;
                    seg_idx += 1;
                }
            }
        }

        // Byte offset of the winning segment relative to `start`, including
        // the stop codon in the length (Java: len = (parts[j].length()+1)*3).
        let mut off = max_frame;
        {
            let n_aa = if max_frame + 3 <= span.len() { (span.len() - max_frame) / 3 } else { 0 };
            let mut seg_start = 0usize;
            let mut seg_idx = 0usize;
            for k in 0..=n_aa {
                let is_stop = k < n_aa
                    && rust_annotale::translate_codon(&span[max_frame + k * 3..max_frame + k * 3 + 3])
                        == STOP;
                if k == n_aa || is_stop {
                    if seg_idx == max_part {
                        break;
                    }
                    off += (k - seg_start + 1) * 3; // segment plus its stop codon
                    seg_start = k + 1;
                    seg_idx += 1;
                }
            }
        }
        let mut len = (max_len_aa + 1) * 3;

        if max_part > 0 {
            // Snap to the first ATG inside the winning ORF.
            let n_codons = len / 3;
            for c in 0..n_codons {
                if &sequence[start + off + c * 3..start + off + c * 3 + 3] == ATG {
                    off += c * 3;
                    len -= c * 3;
                    break;
                }
            }
        } else {
            // ORF touches the span edge: walk upstream codon-wise to an ATG.
            while start + off >= 3
                && &sequence[start + off - 3..start + off] != ATG
            {
                off -= 3;
                len += 3;
            }
        }

        // Extend downstream until the last codon is a stop.
        while start + off + len + 3 <= sequence.len()
            && rust_annotale::translate_codon(&sequence[start + off + len..start + off + len + 3])
                != STOP
        {
            len += 3;
        }
        while start + off + len + 3 > sequence.len() {
            len -= 3;
        }

        (start + off, start + off + len)
    }
}

/// True when `haystack` contains `needle` as a contiguous substring (byte-exact).
fn contains_fragment(haystack: &[u8], needle: &[u8]) -> bool {
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// All RVD diresidues (codons 12/13) from one HMM-aligned block.
///
/// A real TALE repeat opens with the anchor motif `L?{P,Q}` (LTP, LPP,
/// LTQ variants are all observed); the RVD is residues 12/13 (1-based)
/// after that anchor. Block boundaries from the HMM alignment wobble by a
/// few bp around the repeat start, so the anchor is searched within the
/// first few residues. Requiring it near the block start rejects spurious
/// domains over non-canonical terminal repeats, which fixed-stride
/// walkers happily transcribe into garbage diresidues.
fn rvds_from_block(block: &[u8]) -> Vec<String> {
    const ANCHOR_WINDOW: usize = 10;
    const RVD_FIRST_RESIDUE: usize = 11; // 0-based offset of residue 12

    let aa = rust_annotale::translate(block);
    for i in 0..aa.len().min(ANCHOR_WINDOW) {
        if aa[i] == b'L'
            && i + 2 < aa.len()
            && (aa[i + 2] == b'P' || aa[i + 2] == b'Q')
            && aa.len() >= i + RVD_FIRST_RESIDUE + 2
        {
            return vec![format!(
                "{}{}",
                aa[i + RVD_FIRST_RESIDUE] as char,
                aa[i + RVD_FIRST_RESIDUE + 1] as char
            )];
        }
    }
    Vec::new()
}

/// Greedy peak picking: repeatedly take the global maximum, zero out its
/// +-`suppress` neighborhood, stop below `threshold`. Returns (position, score).
fn pick_peaks(vals: &mut [f32], suppress: usize, threshold: f32) -> Vec<(usize, f32)> {
    let mut peaks = Vec::new();
    loop {
        let mut max_idx = usize::MAX;
        let mut max_val = threshold;
        for (i, &v) in vals.iter().enumerate() {
            if v > max_val {
                max_val = v;
                max_idx = i;
            }
        }
        if max_idx == usize::MAX {
            break;
        }
        peaks.push((max_idx, max_val));
        let lo = max_idx.saturating_sub(suppress);
        let hi = std::cmp::min(vals.len(), max_idx + suppress);
        for v in &mut vals[lo..hi] {
            *v = 0.0;
        }
    }
    peaks
}

/// Cluster position-sorted peaks: a gap > `gap_bp` between peak starts starts
/// a new region spanning [first_peak, last_peak + consensus_len].
fn cluster_peaks(
    peaks: &[(usize, f32)],
    gap_bp: usize,
    consensus_len: usize,
) -> Vec<(usize, usize, f32)> {
    let mut sorted: Vec<&(usize, f32)> = peaks.iter().collect();
    sorted.sort_by_key(|p| p.0);

    let mut regions = Vec::new();
    let mut first: Option<usize> = None;
    let mut last_pos = 0usize;
    let mut best_score = f32::NEG_INFINITY;

    for &(pos, score) in sorted {
        match first {
            None => {
                first = Some(pos);
                last_pos = pos;
                best_score = score;
            }
            Some(f) => {
                if pos - last_pos > gap_bp {
                    regions.push((f, last_pos + consensus_len, best_score));
                    first = Some(pos);
                    best_score = score;
                } else {
                    best_score = best_score.max(score);
                }
                last_pos = pos;
            }
        }
    }
    if let Some(f) = first {
        regions.push((f, last_pos + consensus_len, best_score));
    }
    regions
}

/// Remove regions fully contained in another region of the same strand
/// (Java NHMMer post-filter). Larger containing region is kept.
fn dedup_nested(regions: Vec<TALERegion>) -> Vec<TALERegion> {
    let mut keep = vec![true; regions.len()];
    for (i, a) in regions.iter().enumerate() {
        for (j, b) in regions.iter().enumerate() {
            if i != j
                && a.strand == b.strand
                && b.start <= a.start
                && a.end <= b.end
                && !(a.start == b.start && a.end == b.end && i < j)
            {
                keep[i] = false;
            }
        }
    }
    regions
        .into_iter()
        .zip(keep)
        .filter_map(|(r, k)| k.then_some(r))
        .collect()
}

/// Scan range for terminus refinement, ported from Java getBestTerminus.
/// Returns (first_index, limit_index); the caller steps by 5 toward `limit`.
fn terminus_scan_range(
    is_start: bool,
    start: i64,
    end: i64,
    seq_len: usize,
    window_len: usize,
    consensus_len: usize,
) -> Option<(i64, i64)> {
    let w = window_len as i64;
    let tenth = (0.1 * consensus_len as f64).round() as i64;
    if is_start {
        let begin = start - w + tenth;
        let floor = std::cmp::max(0, start - w - 200);
        if begin < 0 || begin < floor {
            None
        } else {
            Some((begin, floor))
        }
    } else {
        let begin = end - tenth;
        let limit = std::cmp::min(end + 200, seq_len as i64 - w + 1);
        if begin >= limit {
            None
        } else {
            Some((begin, limit))
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let use_kmer_filter = !args.no_kmer_filter;
    let finder =
        TALEFinder::new(&args.hmm_dir, use_kmer_filter, args.min_kmers, args.threshold)?;

    let mut tsv_writer = match &args.output {
        Some(path) => {
            let file = File::create(path)
                .with_context(|| format!("Failed to create output file {}", path))?;
            let mut writer = BufWriter::new(file);
            writeln!(writer, "contig\tstrand\tstart\tend\tscore\tpseudo\trvds")?;
            Some(writer)
        }
        None => None,
    };

    println!(
        "Scanning {} for TALE effectors{} (k-mer filter: {}, peak threshold: {:.1} bits)...",
        args.input,
        if args.metagenome { " [metagenomic]" } else { "" },
        if use_kmer_filter { "enabled" } else { "disabled" },
        finder.peak_threshold_bits
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
            total_tales += process_batch(&batch, &finder, args.metagenome, tsv_writer.as_mut())?;
            batch.clear();
        }
    }

    if !batch.is_empty() {
        total_scanned += batch.len();
        total_tales += process_batch(&batch, &finder, args.metagenome, tsv_writer.as_mut())?;
    }

    println!("Total sequences scanned: {}", total_scanned);
    println!("Total potential TALE effectors found: {}", total_tales);
    Ok(())
}

/// Results sink: `Some` writes machine-readable TSV rows, `None` prints
/// human-readable tables to stdout.
type TsvWriter<'a> = Option<&'a mut BufWriter<File>>;

fn process_batch(
    batch: &[SeqRecord],
    finder: &TALEFinder,
    metagenome: bool,
    tsv_out: TsvWriter<'_>,
) -> Result<usize> {
    let results: Vec<(String, Vec<TALERegion>)> = batch
        .par_iter()
        .map(|record| {
            let id = record.id.clone();
            (id.clone(), finder.scan_sequence(&id, &record.seq))
        })
        .collect();

    let mut tales_found = 0;
    match tsv_out {
        Some(writer) => {
            for (id, matches) in results {
                for region in matches {
                    tales_found += 1;
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                        id,
                        region.strand,
                        region.start,
                        region.end,
                        region.score,
                        region.is_pseudo,
                        region.rvds
                    )?;
                }
            }
        }
        None => {
            for (id, matches) in results {
                if !matches.is_empty() {
                    tales_found += matches.len();
                    println!("\nFound {} potential TALE effectors in {}", matches.len(), id);
                    print_table(&matches);
                } else if !metagenome {
                    println!("No TALE effectors found in {}", id);
                }
            }
        }
    }
    Ok(tales_found)
}

fn print_table(matches: &[TALERegion]) {
    println!(
        "{:<12} {:<6} {:<12} {:<12} {:<10} {:<8}",
        "Strand", "Start", "End", "Score", "Pseudo", "RVDs"
    );
    println!("{}", "-".repeat(80));
    for region in matches {
        println!(
            "{:<12} {:<6} {:<12} {:<12} {:<10} {:<8}",
            region.strand, region.start, region.end, region.score, region.is_pseudo, region.rvds
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pick_peaks_finds_two_separated_maxima_and_stops_below_threshold() {
        let mut vals = vec![0.0; 1000];
        vals[100] = 50.0;
        vals[110] = 40.0; // suppressed neighbor of peak 1
        vals[800] = 45.0;
        vals[999] = 5.0; // below threshold, never picked
        let peaks = pick_peaks(&mut vals, 56, 38.6);
        let positions: Vec<usize> = peaks.iter().map(|p| p.0).collect();
        assert_eq!(positions, vec![100, 800]);
        assert_eq!(peaks[0].1, 50.0);
        assert_eq!(peaks[1].1, 45.0);
    }

    #[test]
    fn pick_peaks_returns_empty_when_all_below_threshold() {
        let mut vals = vec![10.0; 100];
        assert!(pick_peaks(&mut vals, 5, 38.6).is_empty());
    }

    #[test]
    fn cluster_peaks_splits_on_gap_and_extends_by_consensus_len() {
        let peaks = [(100, 40.0), (200, 41.0), (900, 42.0)];
        let regions = cluster_peaks(&peaks, 500, 102);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0], (100, 302, 41.0)); // 200+102, best score kept
        assert_eq!(regions[1], (900, 1002, 42.0));
    }

    #[test]
    fn cluster_peaks_unsorted_input_is_sorted_first() {
        let peaks = [(900, 40.0), (100, 41.0)];
        let regions = cluster_peaks(&peaks, 500, 102);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].0, 100);
    }

    #[test]
    fn dedup_nested_keeps_container_drops_contained_same_strand() {
        let mk = |strand: char, start: usize, end: usize| TALERegion {
            strand,
            start,
            end,
            score: 1.0,
            is_pseudo: false,
            rvds: String::new(),
        };
        let regions = vec![
            mk('+', 100, 900),
            mk('+', 200, 400),  // contained -> dropped
            mk('-', 250, 350),  // other strand -> kept
            mk('+', 5000, 6000), // disjoint -> kept
        ];
        let out = dedup_nested(regions);
        assert_eq!(out.len(), 3);
        assert_eq!((out[0].start, out[0].end), (100, 900));
        assert_eq!(out[0].strand, '+');
        assert_eq!(out[1].strand, '-');
        assert_eq!((out[2].start, out[2].end), (5000, 6000));
    }

    #[test]
    fn terminus_scan_range_matches_java_bounds() {
        let w = 112;
        let cons = 102;
        // N-term: begin = start - w + 10 = 1000-112+10 = 898; floor = 1000-112-200 = 688
        assert_eq!(terminus_scan_range(true, 1000, 1200, 5000, w, cons), Some((898, 688)));
        // N-term near contig start: begin < 0 -> None
        assert_eq!(terminus_scan_range(true, 50, 200, 5000, w, cons), None);
        // C-term: begin = end-10; limit = min(end+200, len-w+1)
        assert_eq!(
            terminus_scan_range(false, 1000, 1200, 1500, w, cons),
            Some((1190, 1389))
        );
        // C-term clamped by contig length (len-w+1 = 1139): non-empty but clamped
        assert_eq!(
            terminus_scan_range(false, 1000, 1100, 1250, w, cons),
            Some((1090, 1139))
        );
        // C-term where begin >= clamped limit -> empty loop -> None
        // (Java: positions empty -> getBestTerminus returns null)
        assert_eq!(
            terminus_scan_range(false, 1000, 1200, 1300, w, cons),
            None
        );
    }

    #[test]
    fn contains_fragment_is_exact_byte_match() {
        assert!(contains_fragment(b"AAAACGTCCCGGG", b"ACGTCCC"));
        assert!(!contains_fragment(b"aaaacgtcccggg", b"ACGTCCC"));
    }

    /// One TALE repeat: 34 aa (102 bp), RVD diresidue at residues 12/13
    /// (1-based) = 0-based codon indices 11/12.
    fn tale_repeat(rvd12: &[u8], rvd13: &[u8]) -> Vec<u8> {
        let mut codons: Vec<&[u8]> = vec![b"CTG", b"ACG", b"CCG"]; // L T P
        codons.resize(11, b"GCT"); // filler A at codons 4..=10
        codons.push(rvd12); // residue 12
        codons.push(rvd13); // residue 13
        codons.resize(34, b"GCT"); // filler through residue 34
        codons.concat()
    }

    #[test]
    fn rvds_from_block_reads_diresidue_via_ltp_anchor() {
        let block = tale_repeat(b"CAT", b"GAT"); // H D
        assert_eq!(rvds_from_block(&block), vec!["HD"]);
    }

    #[test]
    fn rvds_from_block_accepts_lpp_and_ltq_anchor_variants() {
        // LPP start (repeat 3 of PthXo1) and LTQ start both encode RVDs.
        let mut lpp: Vec<&[u8]> = vec![b"CTC", b"CCA", b"CCA"]; // L P P
        lpp.resize(11, b"GCT");
        lpp.push(b"CAT");
        lpp.push(b"GAT");
        lpp.resize(34, b"GCT");
        assert_eq!(rvds_from_block(&lpp.concat()), vec!["HD"]);

        let mut ltq: Vec<&[u8]> = vec![b"CTG", b"ACA", b"CAA"]; // L T Q
        ltq.resize(11, b"GCT");
        ltq.push(b"CAT");
        ltq.push(b"GAT");
        ltq.resize(34, b"GCT");
        assert_eq!(rvds_from_block(&ltq.concat()), vec!["HD"]);
    }

    #[test]
    fn rvds_from_block_rejects_blocks_without_start_anchor() {
        // Anchors within a few residues of the block start are tolerated
        // (HMMER boundary wobble); anchors further in are rejected.
        let mut near = tale_repeat(b"CAT", b"GAT");
        near.splice(0..0, vec![b'A'; 6]); // anchor lands at aa index 2 -> accepted
        assert_eq!(rvds_from_block(&near), vec!["HD"]);

        let mut far = tale_repeat(b"CAT", b"GAT");
        far.splice(0..0, vec![b'A'; 36]); // anchor at aa index 12 -> beyond window
        assert!(rvds_from_block(&far).is_empty());
        assert!(rvds_from_block(&far[1..]).is_empty()); // off-by-one frame
    }
}
