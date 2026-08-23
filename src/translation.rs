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

/// Index of the first LTP motif in an amino acid sequence, if any.
pub fn find_first_ltp(aa: &[u8]) -> Option<usize> {
    (0..aa.len().saturating_sub(2)).find(|&i| aa[i] == b'L' && aa[i + 1] == b'T' && aa[i + 2] == b'P')
}

/// Check if byte slice represents DNA
pub fn is_dna(seq: &[u8]) -> bool {
    if seq.contains(&b'-') || seq.contains(&b',') {
        return false;
    }
    seq.iter()
        .all(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N' | b'U'))
}

/// Convert DNA sequence to RVD list
///
/// Anchors at the first LTP motif and walks a fixed 102-bp stride. The RVD
/// is residues 12/13 of each repeat (1-based), i.e. anchor+11/anchor+12
/// in 0-based indices relative to the L.
pub fn dna_to_rvds(seq: &[u8]) -> Vec<String> {
    let mut rvds = Vec::new();
    let aa_seq = translate(seq);
    let start_idx = find_first_ltp(&aa_seq);
    if let Some(start) = start_idx {
        let mut curr = start * 3;
        while curr + 102 <= seq.len() {
            let repeat_dna = &seq[curr..curr + 102];
            let repeat_aa = translate(repeat_dna);
            if repeat_aa.len() >= 13 {
                let rvd = format!("{}{}", repeat_aa[11] as char, repeat_aa[12] as char);
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Golden standard genetic code table, codon bases ordered T, C, A, G
    /// (index = b0*16 + b1*4 + b2). '*' marks a stop codon.
    const GENETIC_CODE: &str =
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    #[test]
    fn translate_codon_matches_standard_genetic_code() {
        let bases = b"TCAG";
        let mut idx = 0;
        for &b0 in bases {
            for &b1 in bases {
                for &b2 in bases {
                    let expected = GENETIC_CODE.as_bytes()[idx];
                    assert_eq!(
                        translate_codon(&[b0, b1, b2]),
                        expected,
                        "codon {}{}{}",
                        b0 as char,
                        b1 as char,
                        b2 as char
                    );
                    idx += 1;
                }
            }
        }
    }

    #[test]
    fn translate_ignores_trailing_partial_codon_and_uppercases() {
        assert_eq!(translate(b"ATGGCT"), b"MA");
        assert_eq!(translate(b"atggct"), b"MA");
        assert_eq!(translate(b"ATGGCTA"), b"MA"); // trailing base dropped
        assert!(translate(b"").is_empty());
    }

    /// One TALE repeat: 34 aa (102 bp) with the RVD diresidue at aa positions
    /// 12/13. Build three repeats encoding HD and expect ["HD", "HD", "HD"].
    fn tale_repeat(rvd12: &[u8], rvd13: &[u8]) -> Vec<u8> {
        let mut codons: Vec<&[u8]> = vec![b"CTG", b"ACG", b"CCG"]; // L T P
        codons.resize(11, b"GCT"); // filler A at residues 4..=10
        codons.push(rvd12); // residue 12
        codons.push(rvd13); // residue 13
        codons.resize(34, b"GCT"); // filler through residue 34
        codons.concat()
    }

    #[test]
    fn dna_to_rvds_extracts_diresidues_from_three_repeats() {
        let mut seq = Vec::new();
        for _ in 0..3 {
            seq.extend(tale_repeat(b"CAT", b"GAT")); // H D
        }
        assert_eq!(seq.len(), 306);
        assert_eq!(dna_to_rvds(&seq), vec!["HD", "HD", "HD"]);
    }

    #[test]
    fn dna_to_rvds_returns_empty_without_ltp_motif() {
        let seq = vec![b'A'; 306];
        assert!(dna_to_rvds(&seq).is_empty());
    }

    #[test]
    fn parse_rvd_sequence_handles_all_delimiters() {
        assert_eq!(parse_rvd_sequence("NI-NG-NI"), vec!["NI", "NG", "NI"]);
        assert_eq!(parse_rvd_sequence("NI,NG,NI"), vec!["NI", "NG", "NI"]);
        assert_eq!(parse_rvd_sequence("NINGNI"), vec!["NI", "NG", "NI"]);
        assert_eq!(parse_rvd_sequence(" NI NG "), vec!["NI", "NG"]);
        assert_eq!(parse_rvd_sequence("NIG"), vec!["NI", "G"]); // odd trailing char
        assert!(parse_rvd_sequence("").is_empty());
    }

    #[test]
    fn is_dna_rejects_protein_and_gap_characters() {
        assert!(is_dna(b"ACGTNacgtn"));
        assert!(is_dna(b"")); // empty treated as DNA by contract
        assert!(!is_dna(b"LTPQVVAIAS"));
        assert!(!is_dna(b"AC-GT"));
        assert!(!is_dna(b"AC,GT"));
    }
}
