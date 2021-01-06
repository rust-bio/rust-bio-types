// Copyright 2021 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! This module provides a representation for read pair orientation, which is
//! a readout of read mapping.

use strum_macros::{AsRefStr, Display};

/// Representation of read pair orientation
/// (e.g. F1R2 means that the forward read comes first on the reference contig,
/// followed by the reverse read, on the same contig).
///
/// This enum can be pretty-printed into a readable string repesentation:
///
/// ```rust
/// use bio_types::sequencing::read_pair_orientation::ReadPairOrientation;
///
/// // format into string
/// println!("{}", ReadPairOrientation::F1R2);
/// // obtain string via `AsRef<&'static str>`
/// assert_eq!(ReadPairOrientation::R1F2.as_ref(), "R1F2");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, AsRefStr, Display)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum ReadPairOrientation {
    F1R2,
    F2R1,
    R1F2,
    R2F1,
    F1F2,
    R1R2,
    F2F1,
    R2R1,
    None,
}
