use genome;
use sequence::{Base, Sequence};

/// A trait for providing variant information. This can e.g. be implemented by file readers.
pub trait AbstractVariant: genome::AbstractLocus {
    fn kind(&self) -> &Kind;
}

/// Possible genomic variants.
#[derive(Debug, PartialEq, Eq, Clone, Hash, Serialize, Deserialize)]
pub enum Kind {
    SNV(Base),
    MNV(Sequence),
    Insertion(Sequence),
    Deletion(genome::Length),
    Duplication(genome::Length),
    Inversion(genome::Length),
    None,
}

impl Kind {
    /// Return variant length.
    pub fn len(&self) -> genome::Length {
        match self {
            &Kind::SNV(_) => 1,
            &Kind::MNV(ref s) => s.len() as u64,
            &Kind::Insertion(ref s) => s.len() as u64,
            &Kind::Deletion(l) => l,
            &Kind::Duplication(l) => l,
            &Kind::Inversion(l) => l,
            &Kind::None => 1,
        }
    }
}
