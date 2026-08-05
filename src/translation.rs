use std::collections::HashMap;

pub struct Translator {
    table: HashMap<[u8; 3], u8>,
}

impl Default for Translator {
    fn default() -> Self {
        Self::new()
    }
}

impl Translator {
    pub fn new() -> Self {
        let mut table = HashMap::new();
        let codes = [
            (*b"TTT", b'F'), (*b"TTC", b'F'), (*b"TTA", b'L'), (*b"TTG", b'L'),
            (*b"TCT", b'S'), (*b"TCC", b'S'), (*b"TCA", b'S'), (*b"TCG", b'S'),
            (*b"TAT", b'Y'), (*b"TAC", b'Y'), (*b"TAA", b'*'), (*b"TAG", b'*'),
            (*b"TGT", b'C'), (*b"TGC", b'C'), (*b"TGA", b'*'), (*b"TGG", b'W'),
            (*b"CTT", b'L'), (*b"CTC", b'L'), (*b"CTA", b'L'), (*b"CTG", b'L'),
            (*b"CCT", b'P'), (*b"CCC", b'P'), (*b"CCA", b'P'), (*b"CCG", b'P'),
            (*b"CAT", b'H'), (*b"CAC", b'H'), (*b"CAA", b'Q'), (*b"CAG", b'Q'),
            (*b"CGT", b'R'), (*b"CGC", b'R'), (*b"CGA", b'R'), (*b"CGG", b'R'),
            (*b"ATT", b'I'), (*b"ATC", b'I'), (*b"ATA", b'I'), (*b"ATG", b'M'),
            (*b"ACT", b'T'), (*b"ACC", b'T'), (*b"ACA", b'T'), (*b"ACG", b'T'),
            (*b"AAT", b'N'), (*b"AAC", b'N'), (*b"AAA", b'K'), (*b"AAG", b'K'),
            (*b"AGT", b'S'), (*b"AGC", b'S'), (*b"AGA", b'R'), (*b"AGG", b'R'),
            (*b"GTT", b'V'), (*b"GTC", b'V'), (*b"GTA", b'V'), (*b"GTG", b'V'),
            (*b"GCT", b'A'), (*b"GCC", b'A'), (*b"GCA", b'A'), (*b"GCG", b'A'),
            (*b"GAT", b'D'), (*b"GAC", b'D'), (*b"GAA", b'E'), (*b"GAG", b'E'),
            (*b"GGT", b'G'), (*b"GGC", b'G'), (*b"GGA", b'G'), (*b"GGG", b'G'),
        ];
        for (codon, aa) in codes {
            table.insert(codon, aa);
        }
        Translator { table }
    }

    pub fn translate(&self, seq: &[u8]) -> Vec<u8> {
        let mut protein = Vec::with_capacity(seq.len() / 3);
        for chunk in seq.chunks_exact(3) {
            let upper_chunk = [
                chunk[0].to_ascii_uppercase(),
                chunk[1].to_ascii_uppercase(),
                chunk[2].to_ascii_uppercase(),
            ];
            protein.push(*self.table.get(&upper_chunk).unwrap_or(&b'X'));
        }
        protein
    }
}
