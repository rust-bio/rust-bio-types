//! Data types for positions and regions on named sequences
//! (e.g. chromosomes), useful for annotating features in a genome.

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
fn break_start_end(s: &str) -> Result<(isize, isize), ParseAnnotError> {
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
    let endstr = s
        .get((breakpt + 1)..(s.len()))
        .ok_or_else(|| ParseAnnotError::NoPosition)?;
    let end = endstr.parse().map_err(|e| ParseAnnotError::ParseInt(e))?;
    Ok((start, end))
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
