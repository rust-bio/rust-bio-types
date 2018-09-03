// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Data types for strand information on annotations.

use std::fmt::{self, Display, Formatter};
use std::marker;

/// Strand information.
#[derive(Debug, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    /// Returns a `Strand` enum representing the given char.
    ///
    /// The mapping is as follows:
    ///     * '+', 'f', or 'F' becomes `Strand::Forward`
    ///     * '-', 'r', or 'R' becomes `Strand::Reverse`
    ///     * '.', '?' becomes `Strand::Unknown`
    ///     * Any other inputs will return an `Err(StrandError::InvalidChar)`
    pub fn from_char(strand_char: &char) -> Result<Strand, StrandError> {
        match *strand_char {
            '+' | 'f' | 'F' => Ok(Strand::Forward),
            '-' | 'r' | 'R' => Ok(Strand::Reverse),
            '.' | '?' => Ok(Strand::Unknown),
            invalid => Err(StrandError::InvalidChar(invalid)),
        }
    }

    pub fn is_unknown(&self) -> bool {
        if let Strand::Unknown = *self {
            true
        } else {
            false
        }
    }
}

impl Strandedness for Strand {
    fn reverse(&self) -> Self {
        match *self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
            Strand::Unknown => Strand::Unknown,
        }
    }

    fn try_req_strand(&self) -> Option<ReqStrand> {
        match *self {
            Strand::Forward => Some(ReqStrand::Forward),
            Strand::Reverse => Some(ReqStrand::Reverse),
            Strand::Unknown => None,
        }
    }

    fn strand_symbol(&self) -> &str {
        match *self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        }
    }

    fn break_pos_strand(posstr: &str) -> Option<(&str, Self)> {
        if let Some((pos, reqstr)) = ReqStrand::break_pos_strand(posstr) {
            Some((pos, Self::from(reqstr)))
        } else {
            Some((posstr, Strand::Unknown))
        }
    }
}

impl PartialEq for Strand {
    /// Returns true if both are `Forward` or both are `Reverse`, otherwise returns false.
    fn eq(&self, other: &Strand) -> bool {
        match (self, other) {
            (&Strand::Forward, &Strand::Forward) => true,
            (&Strand::Reverse, &Strand::Reverse) => true,
            _ => false,
        }
    }
}

impl Same for Strand {
    fn same(&self, s1: &Self) -> bool {
        match (*self, *s1) {
            (Strand::Forward, Strand::Forward) => true,
            (Strand::Reverse, Strand::Reverse) => true,
            (Strand::Unknown, Strand::Unknown) => true,
            _ => false,
        }
    }
}

impl Display for Strand {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match *self {
            Strand::Unknown => Ok(()),
            _ => write!(f, "({})", self.strand_symbol()),
        }
    }
}

impl From<ReqStrand> for Strand {
    fn from(rstr: ReqStrand) -> Self {
        match rstr {
            ReqStrand::Forward => Strand::Forward,
            ReqStrand::Reverse => Strand::Reverse,
        }
    }
}

/// Strand information for annotations that require a strand.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Ord, PartialOrd, Copy)]
pub enum ReqStrand {
    Forward,
    Reverse,
}

impl Same for ReqStrand {
    fn same(&self, s1: &Self) -> bool {
        self == s1
    }
}

impl Display for ReqStrand {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "({})", self.strand_symbol())
    }
}

impl Strandedness for ReqStrand {
    fn reverse(&self) -> Self {
        match *self {
            ReqStrand::Forward => ReqStrand::Reverse,
            ReqStrand::Reverse => ReqStrand::Forward,
        }
    }

    fn try_req_strand(&self) -> Option<ReqStrand> {
        Some(*self)
    }

    fn strand_symbol(&self) -> &str {
        match *self {
            ReqStrand::Forward => "+",
            ReqStrand::Reverse => "-",
        }
    }

    fn break_pos_strand(posstr: &str) -> Option<(&str, Self)> {
        if posstr.ends_with("(+)") {
            if let Some(pos) = posstr.get(0..(posstr.len() - 3)) {
                return Some((pos, ReqStrand::Forward));
            }
        } else if posstr.ends_with("(-)") {
            if let Some(pos) = posstr.get(0..(posstr.len() - 3)) {
                return Some((pos, ReqStrand::Reverse));
            }
        }

        None
    }
}

impl KnownStrandedness for ReqStrand {
    fn req_strand(&self) -> ReqStrand {
        *self
    }
}

/// Strand information for annotations that definitively have no
/// strand information.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Ord, PartialOrd, Copy)]
pub enum NoStrand {
    Unknown,
}

impl Same for NoStrand {
    fn same(&self, _s1: &Self) -> bool {
        true
    }
}

impl Display for NoStrand {
    fn fmt(&self, _f: &mut Formatter) -> fmt::Result {
        Ok(())
    }
}

impl Strandedness for NoStrand {
    fn reverse(&self) -> Self {
        NoStrand::Unknown
    }
    fn try_req_strand(&self) -> Option<ReqStrand> {
        None
    }
    fn strand_symbol(&self) -> &str {
        "."
    }
    fn break_pos_strand(posstr: &str) -> Option<(&str, Self)> {
        Some((posstr, NoStrand::Unknown))
    }
}

/// Equality-like operator for comparing strand information. Unknown
/// strands are not equal, but they are the "same" as other unknown
/// strands.
pub trait Same {
    /// Indicate when two strands are the "same" -- two
    /// unknown/unspecified strands are the "same" but are not equal.
    fn same(&self, &Self) -> bool;
}

impl<T> Same for Option<T>
where
    T: Same,
{
    fn same(&self, s1: &Self) -> bool {
        match (self, s1) {
            (&Option::None, &Option::None) => true,
            (&Option::Some(ref x), &Option::Some(ref x1)) => x.same(x1),
            (_, _) => false,
        }
    }
}

/// Strand information for annotations on sequences
pub trait Strandedness {
    /// Reversed strand
    fn reverse(&self) -> Self;

    /// Convert into a `ReqStrand`, indicating a specific, known
    /// strand. The return is wrapped in an `Option`, so unknown
    /// strands will return `None`.
    fn try_req_strand(&self) -> Option<ReqStrand>;

    /// Symbol denoting the strand. By convention, in BED and GFF
    /// files, the forward strand is `+`, the reverse strand is `-`,
    /// and unknown or unspecified strands are `.`.
    fn strand_symbol(&self) -> &str;

    /// Parse a strand designator from the end of a location display
    /// string.
    ///
    /// Strand information can be represented by a `(+)` or `(-)` at
    /// the end of a display string. This function extracts that
    /// strand information if it's present, and returns a pair
    /// comprising the rest of the string and the parsed strand.
    ///
    /// When strand information is required but it isn't found in the
    /// display string, then `None` is returned.
    ///
    /// When strand information is optional and it isn't found, then
    /// the entire input string is returned along with the unknown
    /// strand designation.
    fn break_pos_strand(&str) -> Option<(&str, Self)>
    where
        Self: marker::Sized;
}

/// Strand information when the strand is guaranteed to be know.
pub trait KnownStrandedness: Strandedness {
    /// Convert into a `ReqStrand`, indicating a specific, known
    /// strand. As strandedness is guaranteed, this cannot fail.
    fn req_strand(&self) -> ReqStrand;

    fn product<S>(&self, s: S) -> S
    where
        S: Strandedness,
    {
        match self.req_strand() {
            ReqStrand::Forward => s,
            ReqStrand::Reverse => s.reverse(),
        }
    }
}

quick_error! {
    #[derive(Debug, Clone)]
    pub enum StrandError {
        InvalidChar(invalid_char: char) {
            description("invalid character for strand conversion")
            display("character {:?} can not be converted to a Strand", invalid_char)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand() {
        assert_eq!(Strand::from_char(&'+').unwrap(), Strand::Forward);
        assert_eq!(Strand::from_char(&'-').unwrap(), Strand::Reverse);
        assert!(Strand::from_char(&'.').unwrap().is_unknown());
        assert!(Strand::from_char(&'o').is_err());
    }
}
