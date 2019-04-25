/// A DNA base
pub type Base = u8;
/// An amino acid
pub type AminoAcid = u8;
/// A biological sequence
pub type Sequence = Vec<u8>;


pub trait SequenceRead {
    fn name(&self) -> &[u8];
    fn base(&self, i: usize) -> u8;
    fn base_qual(&self, i: usize) -> u8;
    fn len(&self) -> usize;
}
