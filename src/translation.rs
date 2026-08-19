/// Translate a codon triple to an amino acid
#[inline]
pub fn translate_codon(codon: &[u8]) -> u8 {
    let b0 = codon[0].to_ascii_uppercase();
    let b1 = codon[1].to_ascii_uppercase();
    let b2 = codon[2].to_ascii_uppercase();

    match (b0, b1, b2) {
        (b'G', b'C', _) => b'A',
        (b'T', b'G', b'T' | b'C') => b'C',
        (b'G', b'A', b'T' | b'C') => b'D',
        (b'G', b'A', b'A' | b'G') => b'E',
        (b'T', b'T', b'T' | b'C') => b'F',
        (b'G', b'G', _) => b'G',
        (b'C', b'A', b'T' | b'C') => b'H',
        (b'A', b'T', b'T' | b'C' | b'A') => b'I',
        (b'A', b'A', b'A' | b'G') => b'K',
        (b'T', b'T', b'A' | b'G') | (b'C', b'T', _) => b'L',
        (b'A', b'T', b'G') => b'M',
        (b'A', b'A', b'T' | b'C') => b'N',
        (b'C', b'C', _) => b'P',
        (b'C', b'A', b'A' | b'G') => b'Q',
        (b'C', b'G', _) | (b'A', b'G', b'A' | b'G') => b'R',
        (b'T', b'C', _) | (b'A', b'G', b'T' | b'C') => b'S',
        (b'A', b'C', _) => b'T',
        (b'G', b'T', _) => b'V',
        (b'T', b'G', b'G') => b'W',
        (b'T', b'A', b'T' | b'C') => b'Y',
        (b'T', b'A', b'A' | b'G') | (b'T', b'G', b'A') => b'*',
        _ => b'X',
    }
}

/// Translate a DNA sequence into an amino acid sequence
pub fn translate(dna: &[u8]) -> Vec<u8> {
    let mut aa = Vec::with_capacity(dna.len() / 3);
    for chunk in dna.chunks_exact(3) {
        aa.push(translate_codon(chunk));
    }
    aa
}

/// Zero-cost Translator wrapper
#[derive(Default, Clone, Copy, Debug)]
pub struct Translator;

impl Translator {
    pub fn new() -> Self {
        Self
    }

    #[inline]
    pub fn translate(&self, seq: &[u8]) -> Vec<u8> {
        translate(seq)
    }
}

/// Check if byte slice represents DNA
pub fn is_dna(seq: &[u8]) -> bool {
    if seq.contains(&b'-') || seq.contains(&b',') {
        return false;
    }
    seq.iter().all(|&b| match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' | b'N' | b'U' => true,
        _ => false,
    })
}

/// Convert DNA sequence to RVD list
pub fn dna_to_rvds(seq: &[u8]) -> Vec<String> {
    let mut rvds = Vec::new();
    let aa_seq = translate(seq);
    let mut start_idx = None;
    for i in 0..aa_seq.len().saturating_sub(3) {
        if aa_seq[i] == b'L' && aa_seq[i + 1] == b'T' && aa_seq[i + 2] == b'P' {
            start_idx = Some(i);
            break;
        }
    }
    if let Some(start) = start_idx {
        let mut curr = start * 3;
        while curr + 102 <= seq.len() {
            let repeat_dna = &seq[curr..curr + 102];
            let repeat_aa = translate(repeat_dna);
            if repeat_aa.len() >= 14 {
                let rvd = format!("{}{}", repeat_aa[12] as char, repeat_aa[13] as char);
                rvds.push(rvd);
            }
            curr += 102;
        }
    }
    rvds
}

/// Parse RVD sequence from string representation
pub fn parse_rvd_sequence(seq_str: &str) -> Vec<String> {
    let cleaned: String = seq_str.chars().filter(|c| !c.is_whitespace()).collect();
    if cleaned.contains('-') {
        cleaned.split('-').map(|s| s.to_string()).collect()
    } else if cleaned.contains(',') {
        cleaned.split(',').map(|s| s.to_string()).collect()
    } else {
        let chars: Vec<char> = cleaned.chars().collect();
        let mut rvds = Vec::new();
        let mut i = 0;
        while i < chars.len() {
            if i + 1 < chars.len() {
                rvds.push(format!("{}{}", chars[i], chars[i + 1]));
                i += 2;
            } else {
                rvds.push(chars[i].to_string());
                i += 1;
            }
        }
        rvds
    }
}
