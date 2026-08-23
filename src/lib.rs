pub mod models;
pub mod translation;

pub use models::{
    extract_consensus, extract_kmers, open_sequence_reader, revcomp, SeqRecord, SeqReader,
    TALERegion,
};
pub use translation::{
    dna_to_rvds, find_first_ltp, is_dna, parse_rvd_sequence, translate, translate_codon,
};
