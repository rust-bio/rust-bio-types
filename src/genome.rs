use std::ops::Range;

pub type Position = u32;
pub type Length = u32;

pub trait AbstractInterval {
    /// Identifier for a genomic contig, e.g., a chromosome
    fn contig(&self) -> &str;
    /// Interval on the contig
    fn range(&self) -> Range<Position>;
}

#[derive(Debug, PartialEq, Eq, Clone, Hash, Serialize, Deserialize)]
pub struct Interval {
    contig: String,
    range: Range<Position>,
}

impl AbstractInterval for Interval {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn range(&self) -> Range<Position> {
        self.range.clone()
    }
}

pub trait AbstractLocus {
    /// Identifier for a genomic contig, e.g., a chromosome
    fn contig(&self) -> &str;
    /// Position on the contig
    fn pos(&self) -> Position;
}

#[derive(Debug, PartialEq, Eq, Clone, Hash, Serialize, Deserialize)]
pub struct Locus {
    contig: String,
    pos: Position,
}

impl AbstractLocus for Locus {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn pos(&self) -> Position {
        self.pos
    }
}
