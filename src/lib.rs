pub mod models;
pub mod translation;

pub use models::{open_sequence_reader, SeqRecord, SeqReader, TALERegion};
pub use translation::{
    dna_to_rvds, is_dna, parse_rvd_sequence, translate, translate_codon, Translator,
};
