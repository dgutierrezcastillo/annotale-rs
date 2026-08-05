#[derive(Debug, Clone)]
pub struct TALERegion {
    pub strand: char,
    pub start: usize,
    pub end: usize,
    pub score: f32,
    pub is_pseudo: bool,
    pub rvds: String,
}
