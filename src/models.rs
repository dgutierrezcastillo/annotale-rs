use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct TALERegion {
    pub strand: char,
    pub start: usize,
    pub end: usize,
    pub score: f32,
    pub is_pseudo: bool,
    pub rvds: String,
}

#[derive(Clone, Debug)]
pub struct SeqRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

enum Format {
    Fasta,
    Fastq,
    Unknown,
}

pub struct SeqReader<R: BufRead> {
    reader: R,
    format: Format,
    buffer: String,
}

impl<R: BufRead> SeqReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            format: Format::Unknown,
            buffer: String::new(),
        }
    }

    pub fn next_record(&mut self) -> Result<Option<SeqRecord>> {
        if let Format::Unknown = self.format {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(None);
            }
            if buf[0] == b'>' {
                self.format = Format::Fasta;
            } else if buf[0] == b'@' {
                self.format = Format::Fastq;
            } else {
                return Err(anyhow::anyhow!(
                    "Invalid sequence format. First character must be '>' or '@', got '{}'",
                    buf[0] as char
                ));
            }
        }

        match self.format {
            Format::Fasta => {
                let mut header = String::new();
                if !self.buffer.is_empty() {
                    header = self.buffer.clone();
                    self.buffer.clear();
                } else {
                    let bytes_read = self.reader.read_line(&mut header)?;
                    if bytes_read == 0 {
                        return Ok(None);
                    }
                }

                let trimmed_header = header.trim();
                if !trimmed_header.starts_with('>') {
                    return Err(anyhow::anyhow!("FASTA header must start with '>'"));
                }
                let id = trimmed_header[1..].split_whitespace().next().unwrap_or("").to_string();

                let mut seq_bytes = Vec::new();
                let mut line = String::new();
                loop {
                    line.clear();
                    let bytes = self.reader.read_line(&mut line)?;
                    if bytes == 0 {
                        break;
                    }
                    if line.starts_with('>') {
                        self.buffer = line.clone();
                        break;
                    }
                    let trimmed = line.trim();
                    seq_bytes.extend_from_slice(trimmed.as_bytes());
                }

                Ok(Some(SeqRecord {
                    id,
                    seq: seq_bytes,
                    qual: None,
                }))
            }
            Format::Fastq => {
                let mut header = String::new();
                let bytes_read = self.reader.read_line(&mut header)?;
                if bytes_read == 0 {
                    return Ok(None);
                }
                let trimmed_header = header.trim();
                if !trimmed_header.starts_with('@') {
                    return Err(anyhow::anyhow!("FASTQ header must start with '@'"));
                }
                let id = trimmed_header[1..].split_whitespace().next().unwrap_or("").to_string();

                let mut seq_line = String::new();
                if self.reader.read_line(&mut seq_line)? == 0 {
                    return Err(anyhow::anyhow!("Truncated FASTQ record sequence for {}", id));
                }
                let seq = seq_line.trim().as_bytes().to_vec();

                let mut plus_line = String::new();
                if self.reader.read_line(&mut plus_line)? == 0 {
                    return Err(anyhow::anyhow!("Truncated FASTQ '+' line for {}", id));
                }

                let mut qual_line = String::new();
                if self.reader.read_line(&mut qual_line)? == 0 {
                    return Err(anyhow::anyhow!("Truncated FASTQ quality scores for {}", id));
                }
                let qual = qual_line.trim().as_bytes().to_vec();

                Ok(Some(SeqRecord {
                    id,
                    seq,
                    qual: Some(qual),
                }))
            }
            Format::Unknown => unreachable!(),
        }
    }
}

pub fn open_sequence_reader<P: AsRef<Path>>(path: P) -> Result<SeqReader<BufReader<Box<dyn Read>>>> {
    let p = path.as_ref();
    let file = File::open(p).with_context(|| format!("Failed to open {}", p.display()))?;
    let reader: Box<dyn Read> = if p.extension().is_some_and(|ext| ext == "gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(SeqReader::new(BufReader::new(reader)))
}

/// Extract consensus nucleotide sequence from an HMM profile file
pub fn extract_consensus<P: AsRef<Path>>(hmm_path: P) -> Result<String> {
    let file = File::open(hmm_path.as_ref())?;
    let reader = BufReader::new(file);
    let mut consensus = String::new();
    let mut in_hmm = false;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("HMM ") {
            in_hmm = true;
            continue;
        }
        if !in_hmm {
            continue;
        }

        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }

        if let Ok(_state_idx) = tokens[0].parse::<usize>() {
            if tokens.len() > 6 {
                consensus.push_str(tokens[6]);
            }
        }
    }

    Ok(consensus)
}

/// Reverse complement of a DNA sequence (full IUPAC; unknown bytes -> N).
pub fn revcomp(dna: &[u8]) -> Vec<u8> {
    fn complement(b: u8) -> u8 {
        match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            _ => b'N',
        }
    }
    dna.iter().rev().map(|&b| complement(b)).collect()
}

/// Split consensus sequence into non-overlapping k-mer fragments for fast pre-filtering
pub fn extract_kmers(consensus: &str, k: usize) -> Vec<Vec<u8>> {
    let mut kmers = Vec::new();
    if k == 0 || consensus.len() < k {
        return kmers;
    }
    let upper = consensus.to_ascii_uppercase();
    let bytes = upper.as_bytes();
    for chunk in bytes.chunks_exact(k) {
        kmers.push(chunk.to_vec());
    }
    kmers
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_reader_handles_multiline_and_lookahead_buffer() {
        let input: &[u8] = b">rec1 some description\nACGT\nGTCA\n>rec2\nTTTT\n";
        let mut reader = SeqReader::new(input);

        let r1 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.id, "rec1");
        assert_eq!(r1.seq, b"ACGTGTCA");
        assert!(r1.qual.is_none());

        // Second record exercises the buffered-header path (line stored back
        // when a '>' is hit mid-sequence).
        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r2.id, "rec2");
        assert_eq!(r2.seq, b"TTTT");

        assert!(reader.next_record().unwrap().is_none());
        assert!(reader.next_record().unwrap().is_none()); // stays exhausted
    }

    #[test]
    fn fastq_reader_roundtrip_four_line_records() {
        let input: &[u8] = b"@read1\nACGT\n+\nIIII\n@read2\nGG\n+\nHH\n";
        let mut reader = SeqReader::new(input);

        let r1 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.id, "read1");
        assert_eq!(r1.seq, b"ACGT");
        assert_eq!(r1.qual.unwrap(), b"IIII");

        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r2.id, "read2");
        assert_eq!(r2.qual.unwrap(), b"HH");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn fastq_truncation_is_reported_as_error() {
        let cases: [&[u8]; 3] = [
            b"@read1\n",            // missing seq
            b"@read1\nACGT\n",      // missing '+'
            b"@read1\nACGT\n+\n",   // missing qual
        ];
        for input in cases {
            let mut reader = SeqReader::new(input);
            assert!(reader.next_record().is_err(), "input {:?}", String::from_utf8_lossy(input));
        }
    }

    #[test]
    fn invalid_first_byte_is_rejected() {
        let mut reader = SeqReader::new(&b"not-a-sequence-file\n"[..]);
        let err = reader.next_record().unwrap_err().to_string();
        assert!(err.contains("Invalid sequence format"), "{}", err);
    }

    #[test]
    fn extract_kmers_splits_non_overlapping_uppercase() {
        assert_eq!(
            extract_kmers("abcdefghijKl", 5),
            vec![b"ABCDE".to_vec(), b"FGHIJ".to_vec()]
        );
        // trailing remainder shorter than k is dropped
        assert_eq!(extract_kmers("ABCDEFGH", 5), vec![b"ABCDE".to_vec()]);
        assert!(extract_kmers("AC", 5).is_empty());
        assert!(extract_kmers("ACGT", 0).is_empty());
    }

    #[test]
    fn extract_consensus_reads_state_column_from_hmm_body() {
        let hmm: &[u8] = b"# STOCKHOLM 1.0\n\
                          NAME repeats\n\
                          LENG 2\n\
                          HMM A C G T\n\
                          1 99 4 59 36 0 A 1 2 3\n\
                          2 98 4 59 36 0 C 1 2 3\n\
                          //\n";
        let dir = std::env::temp_dir().join(format!("annotale-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("repeats.hmm");
        std::fs::write(&path, hmm).unwrap();

        let consensus = extract_consensus(&path).unwrap();
        assert_eq!(consensus, "AC");

        std::fs::remove_dir_all(&dir).ok();
    }
}

