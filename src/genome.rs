use std::ops::Range;

pub type Position = u32;
pub type Length = u32;

pub trait AbstractInterval {
    /// Identifier for a genomic contig, e.g., a chromosome
    fn contig(&self) -> &str;
    /// Interval on the contig
    fn range(&self) -> Range<Position>;
}

#[derive(new, Debug, PartialEq, Eq, Clone, Hash, Serialize, Deserialize)]
pub struct Interval {
    contig: String,
    range: Range<Position>,
}

impl Interval {
    /// Mutable reference to interval on the contig
    pub fn range_mut(&mut self) -> &mut Range<Position> {
        &mut self.range
    }
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

#[derive(new, Debug, PartialEq, Eq, Clone, Hash, Serialize, Deserialize)]
pub struct Locus {
    contig: String,
    pos: Position,
}

impl Locus {
    /// Mutable reference to position.
    fn pos_mut(&mut self) -> &mut Position {
        &mut self.pos
    }
}

impl AbstractLocus for Locus {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn pos(&self) -> Position {
        self.pos
    }
}
