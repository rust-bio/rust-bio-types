//! Data types for positions and regions on named sequences
//! (e.g. chromosomes), useful for annotating features in a genome.
//! For example, these data types let you represent that _TMA22_ is on
//! chromosome X, positions 461,829-462,426, on the forward strand. They
//! also allow coordinate math on these annotations, e.g., that
//! position chrX:461,839 is +10 within _TMA22_ and vice versa.
//!
//! This module provides three concrete data types to represent a
//! single position ([`Pos`](pos/Pos.t.html)), a contiguous region
//! ([`Contig`](contig/Contig.t.html)), or a "spliced" region
//! ([`Spliced`](spliced/Spliced.t.html)) consisting of one or more
//! exons separated by introns.  All three data types implement a
//! location trait [`Loc`](loc/Loc.t.html).
//!
//! These data types are generic over the data type used to "name" the
//! annotated reference sequence (e.g., the chromosome name). It's
//! possible to use an owned `String`, an interned `Rc<String>`, or an
//! integer sequence identifier like the "target id" field in a BAM
//! file.
//!
//! These data types are also generic over the kind of strand
//! information in the annotation. This allows annotations with
//! _required_ strand annotation
//! ([`ReqStrand`](../strand/enum.ReqStrand.html)), _optional_ strand
//! annotation ([`Strand`](../strand/enum.Strand.html)), or _no_
//! strand annotation ([`NoStrand`](../strand/enum.NoStrand.html)).
//!
//! The example below shows how to create the _TMA22_ annotation and
//! find where chrX:461,839 falls within this gene.
//! ```
//! # use bio_types::annot::ParseAnnotError;
//! # fn try_main() -> Result<(), ParseAnnotError> {
//! use bio_types::annot::contig::Contig;
//! use bio_types::annot::loc::Loc;
//! use bio_types::annot::pos::Pos;
//! use bio_types::strand::{ReqStrand,NoStrand};
//! let tma22 = "chrX:461829-462426(+)".parse::<Contig<String, ReqStrand>>()?;
//! let p0 = "chrX:461839".parse::<Pos<String, NoStrand>>()?;
//! let p0_into = tma22.pos_into(&p0).unwrap_or_else(|| panic!("p0 not within TMA22"));
//! assert!(p0_into.pos() == 10);
//! # Ok(())
//! # }
//! # fn main() { try_main().unwrap(); }
//! ```

use strand;

pub mod contig;
pub mod loc;
pub mod pos;
pub mod refids;
pub mod spliced;

/// Errors that arise in parsing annotations.
quick_error! {
    #[derive(Debug, Clone)]
    pub enum ParseAnnotError {
        NoColon {
            description("Missing ':'")
        }
        NoRefid {
            description("Missing reference sequence ID")
        }
        NoRange {
            description("Missing position range")
        }
        NoPosition {
            description("Missing position")
        }
        NoStrand {
            description("Missing strand")
        }
        ParseInt(err: ::std::num::ParseIntError) {
            description("Integer parsing error")
        }
        ParseStrand(err: strand::StrandError) {
            description("Strand parsing error")
        }
        Splicing(err: spliced::SplicingError) {
            description("Bad splicing structure")
        }
        EndBeforeStart {
            description("Ending position < starting position")
        }
    }
}

// Break a position display string into a reference name part and the
// "rest"
fn break_refid(s: &str) -> Result<(&str, &str), ParseAnnotError> {
    let breakpt = s.find(':').ok_or_else(|| ParseAnnotError::NoColon)?;
    let refid = s.get(..breakpt).ok_or_else(|| ParseAnnotError::NoRefid)?;
    let rest = s
        .get((breakpt + 1)..)
        .ok_or_else(|| ParseAnnotError::NoPosition)?;
    Ok((refid, rest))
}

// Break a range into a start and an end, and parse each as isize
// individually.  We can unambiguously parse negative starting and
// ending positions by breaking on the first '-' occurring after the
// first character of the range.
fn break_start_end(s: &str) -> Result<(isize, isize, &str), ParseAnnotError> {
    let mut breakpt = s
        .get(1..s.len())
        .ok_or_else(|| ParseAnnotError::NoPosition)?
        .find('-')
        .ok_or_else(|| ParseAnnotError::NoRange)?;
    breakpt += 1;
    let startstr = s
        .get(0..breakpt)
        .ok_or_else(|| ParseAnnotError::NoPosition)?;
    let start = startstr.parse().map_err(|e| ParseAnnotError::ParseInt(e))?;
    let rest = s
        .get((breakpt + 1)..(s.len()))
        .ok_or_else(|| ParseAnnotError::NoPosition)?;
    let endbreak = rest.find(|c: char| !c.is_numeric()).unwrap_or(rest.len());
    let (endstr, leftover) = rest.split_at(endbreak);
    let end = endstr.parse().map_err(|e| ParseAnnotError::ParseInt(e))?;
    Ok((start, end, leftover))
}

/// Errors that arise in maniuplating annotations
quick_error! {
    #[derive(Debug, Clone)]
    pub enum AnnotError {
        NoStrand {
            description("No strand information")
        }
        BadSplicing {
            description("Invalid splicing structure")
        }
    }
}
