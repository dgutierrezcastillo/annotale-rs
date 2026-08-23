# Benchmark Report: AnnoTALE (Java) vs annotale-rs on Xanthomonas oryzae pv. oryzae PXO99ᴬ

**Date**: 2026-08-22 · **Machine**: Apple Silicon Mac, 10 cores, 16 GB RAM, macOS (darwin)

## Setup

| Component | Version | Source |
|---|---|---|
| AnnoTALE CLI | 1.5 (`predict` tool v1.4.2), Zulu JDK 8u504 (ARM64) | jstacs.de / Azul |
| annotale-rs | local clone @ `fc5f414` + P0–P2 fixes, release build (cargo 1.97.1) | github.com/dgutierrezcastillo/annotale-rs |
| Genome | *Xoo* PXO99ᴬ chromosome CP000967.1, 5,240,075 bp | NCBI Entrez |
| HMM profiles | `starts.hmm`, `repeats.hmm`, `ends.hmm` extracted from AnnoTALEcli jar | identical models for both tools |
| Ground truth | 19 TALE genes (incl. 2 pseudogenes) mapped to CP000967.1 | jstacs.de `List_of_TALEs.tsv` |

Note: the Rust port loads only `repeats.hmm`; the Java pipeline uses all three profiles plus internal filtering.

## Protocol

- Both tools ran with default settings on the same input FASTA.
- **Parallel run**: both launched simultaneously under `/usr/bin/time -l` (wall-clock + peak RSS).
- **Sequential run**: back-to-back for uncontended reference timing.
- Correctness: predictions matched against the 19-gene ground truth within ±1,200 bp tolerance; RVD strings derived from each tool's own predicted sequence/protein using a fixed 102-bp stride from the LTP motif (each tool's native approach).

## Performance

| Metric | Java AnnoTALE 1.5 | annotale-rs (release) |
|---|---|---|
| Wall time (sequential) | 82.7 s | **70.7 s** |
| Wall time (parallel pair) | 92.5 s | **78.0 s** |
| User CPU time | 86.3 s | 60.2 s |
| Peak RSS (sequential) | 3.3 GB | **10.7 GB** |
| Peak RSS (parallel) | 2.27 GB | **9.86 GB** |

**Reading**: Rust is ~15% faster in wall-clock and ~30% less CPU, but allocates **~3× more memory** (≈10.7 GB vs 3.3 GB) — pathological for a 5 Mb genome; likely unbanded dynamic-programming allocations inside `hmmer-pure-rs` on chromosome-length sequences.

## Correctness

### Recall & precision vs ground truth

| Metric | Java | Rust |
|---|---|---|
| True TALEs found (of 19) | **19/19** | **19/19** |
| Strand-correct calls | 19/19 | 19/19 |
| False-positive regions | **0** | **76** (+4 degenerate zero-length rows) |
| Pseudo/CDS labeling | 19/19 correct | not directly comparable |

### RVD extraction agreement (19 true loci)

- Mean positional agreement: **55%**
- Identical across ≥8 leading diresidues: 4/19 loci (perfect 100% on the four near-identical *tal7/tal8* paralog loci)
- Java reports systematically more repeats per gene (e.g., 29 vs 22 at tal1) — its predicted CDS extends further downstream.
- **Root cause of divergence**: both implementations extract RVDs by walking a fixed 102-bp stride from the first LTP motif. Natural TALE arrays contain non-canonical repeat lengths, so both walkers drift after indel events; combined with different CDS boundary choices, strings diverge mid-array. This is a shared algorithmic limitation, not evidence that one tool mis-calls biology.

## Update 2026-08-23: Java-faithful NHMMer port (P0+P1+P2 resolved)

annotale-rs @ `31944f1` replaces the genome-wide `repeats.hmm` HMMER scan with a
faithful port of Java's algorithm (see `notes.md` for extracted semantics):

- rolling consensus 10-mer prefilter → per-window HMM scoring (reusable profile)
- greedy peak picking with ±w/2 suppression, 500 bp clustering
- `getBestTerminus` boundary refinement via `starts.hmm`/`ends.hmm`
- `refine()` ORF recovery: longest stop-to-stop ORF in any frame, ATG snap with
  upstream codon walk, downstream stop extension
- nested-region dedup; clusters without resolvable ORF are dropped

| Metric | Java 1.5 | Rust before | **Rust after** |
|---|---|---|---|
| Recall (±1,200 bp) | 19/19 | 19/19 | **19/19** |
| False positives | 0 | 76 | **0** |
| Degenerate rows | 0 | 4 | **0** |
| Wall time (sequential) | 82.7 s | 70.7 s | **52.8 s** |
| Peak RSS (sequential) | 3.3 GB | 10.7 GB | **60.6 MB** |

The memory blowup disappeared structurally: no chromosome-length DP matrices are
allocated any more — only w=112 bp windows pass through the pipeline.

Known deviations from Java (documented in `notes.md`):
- window score = HMMER domain bitscore (EV-corrected) instead of raw log-likelihood ratio;
  peak threshold auto-derived as `consensus_len·ln(1.3)` in bits (`--threshold` overrides);
- terminus extent from HMMER alignment bounds instead of Viterbi match-state run;
- pseudo-gene flag not yet triggered on the two pseudogenic loci Java labels.

Re-run:

```bash
../annotale-rs/target/release/annotale-rs --input data/PXO99A_v1.fna --hmm-dir hmm/ -O results/rust_fixed.tsv
python3 bin/compare_truth.py results/rust_fixed.tsv
```

Artifacts added: `results/rust_fixed.tsv`, `bin/compare_truth.py`, `notes.md`,
`java_classes/NHMMer.java` (+ TALEPredictionTool source, javap dumps).

## Update 2026-08-23 (later): P3 per-repeat HMM-aligned RVD extraction

annotale-rs @ `d9e671d` replaces the fixed-stride RVD walker (shared flaw with
Java) with per-repeat extraction: CDS is re-scanned with `repeats.hmm`, each
domain block is trimmed to the reading frame (phase locked by majority vote),
anchored on an `L?{P,Q}` motif near its start, and residues 12/13 are read.

**Off-by-one discovered en route**: the legacy readers in *both* tools read
0-based anchor+12/+13, i.e. residues 13/14 — every previously published-style
RVD string from either tool was shifted one residue right (canonical `NI/NG`
came out as `IG/GG`). Fixed in `dna_to_rvds`, the scanner, and `analyze`.

Validation against references:
- **Tal2b/PthXo1** (TalBX1): core array matches reviewed UniProt B2SU53 exactly
  through repeat 12, including its hard `LPP` third repeat; remaining divergence
  confined to the known-degenerate C-terminal repeats.
- **TalAE1**: matches the published TalAE-class alignment (Plant Disease 2020,
  LN4 paper) at ~10/13 positions, diverging only where strain orthologs vary.
- vs Java (both sides frame-corrected): leading-8 identity **10/19** (was 7/19);
  mean positional agreement 53%. The low-agreement loci are ones where Java's
  stride demonstrably walks into non-repeat sequence (e.g. tal6a: Java reports
  27 "repeats" ending in `KQ-KQ-RP-DL-GL`; annotale-rs emits 18 anchored ones).

Detection unchanged: recall 19/19, FP 0, degenerate rows 0.

## Conclusions

1. **Detection parity**: the Rust port finds every annotated TALE on PXO99ᴬ — recall is fully on par with the reference implementation.
2. **Specificity gap is the main weakness**: 80% of Rust's calls are false positives (76 extra regions genome-wide, scores dominated by ~4,100-bit and ~21,500-bit summed clusters). Root cause: it scans only `repeats.hmm`, ignoring `starts.hmm`/`ends.hmm` and Java's N/C-terminal validation. Adding start/end profile checks would likely eliminate most FPs.
3. **Performance win is modest, memory is a liability**: ~15% faster wall-clock does not justify 3× memory today; fixing HMMER memory banding matters more than speed.
4. **Robustness bugs observed**: degenerate zero-length predictions (start==end) when no ORF is refined; output rows include them unfiltered.
5. **Practical recommendation**: usable as a fast *screening* pre-filter (all true sites recovered), but not yet as a drop-in annotation replacement until specificity and RVD walking are fixed.

## Reproduction

```bash
# data (see data/, hmm/, bin/)
java -jar bin/AnnoTALEcli-1.5.jar predict g=data/PXO99A_v1.fna s=PXO99A outdir=results/java_full   # JDK 8 required
annotale-rs --input data/PXO99A_v1.fna --hmm-dir hmm/ -O results/rust_full.tsv
python3 results/correctness_summary.txt   # regenerated comparison
```

Artifacts: `results/{java_full,rust_full.tsv,correctness_summary.txt}`, `logs/*.{log}`, this report.
