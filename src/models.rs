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
    let reader: Box<dyn Read> = if p.extension().map_or(false, |ext| ext == "gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(SeqReader::new(BufReader::new(reader)))
}
