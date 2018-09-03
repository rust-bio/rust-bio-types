// Copyright 2017 Nicholas Ingolia
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Spliced region on a named sequence, e.g., the reverse strand of
//! chromosome V, exon #1 at 166,885 through 166,875 and exon #2 at
//! 166,771 through 166,237.

use std::cmp::{max, min};
use std::fmt::{self, Display, Formatter};

use annot::contig::Contig;
use annot::loc::Loc;
use annot::pos::Pos;
use annot::*;
use strand::*;

// The spliced location representation inherently cannot represent
// "bad" splicing structures. Locations comprise a first exon length,
// with a vector of intron-length / exon-length pairs, all type usize
// and hence non-negative.
//
// The InEx data type used to represent intron-exon pairs after the
// first exon enforces an additional strict positivity constraint
// (i.e., >0) on the lengths. This further eliminates degeneracy in
// representations: only one unique set of positive-length exons with
// interleaved positive-length introns can represent a splicing
// structure.
mod inex {
    use std::slice::Iter;

    use super::SplicingError;

    #[derive(Debug, Clone, Hash, PartialEq, Eq)]
    pub struct InEx {
        intron_length: usize,
        exon_length: usize,
    }

    impl InEx {
        pub fn new(intron_length: usize, exon_length: usize) -> Result<Self, SplicingError> {
            if intron_length < 1 {
                Err(SplicingError::IntronLength)
            } else if exon_length < 1 {
                Err(SplicingError::ExonLength)
            } else {
                Ok(InEx {
                    intron_length: intron_length,
                    exon_length: exon_length,
                })
            }
        }
        pub fn intron_length(&self) -> usize {
            self.intron_length
        }
        pub fn exon_length(&self) -> usize {
            self.exon_length
        }
        pub fn length(&self) -> usize {
            self.intron_length + self.exon_length
        }
    }

    // Represent just the start (relative to the start of the location) and length of exons
    // Useful for internal coordinate math
    #[derive(Debug, Clone, Hash, PartialEq, Eq)]
    pub struct Ex {
        start: usize,
        length: usize,
    }

    impl Ex {
        pub fn start(&self) -> usize {
            self.start
        }
        pub fn length(&self) -> usize {
            self.length
        }
        pub fn end(&self) -> usize {
            self.start + self.length
        }
    }

    // Iterator over the Ex exons from lowest to highest coordinate
    pub struct Exes<'a> {
        state: ExesState,
        curr_start: usize,
        rest: Iter<'a, InEx>,
    }

    enum ExesState {
        FirstExon(usize),
        LaterExon,
    }

    impl<'a> Exes<'a> {
        pub fn new(exon_0_length: usize, inexes: &'a Vec<InEx>) -> Self {
            Exes {
                state: ExesState::FirstExon(exon_0_length),
                curr_start: 0,
                rest: inexes.iter(),
            }
        }
    }

    impl<'a> Iterator for Exes<'a> {
        type Item = Ex;

        fn next(&mut self) -> Option<Ex> {
            match self.state {
                ExesState::FirstExon(len) => {
                    let ex = Ex {
                        start: self.curr_start,
                        length: len,
                    };
                    self.curr_start += len;
                    self.state = ExesState::LaterExon;
                    Some(ex)
                }
                ExesState::LaterExon => match self.rest.next() {
                    Some(inex) => {
                        let ex = Ex {
                            start: self.curr_start + inex.intron_length(),
                            length: inex.exon_length(),
                        };
                        self.curr_start += inex.length();
                        Some(ex)
                    }
                    None => None,
                },
            }
        }
    }

    // Represent just the start (relative to the start of the location) and length of introns
    // Useful for internal coordinate math
    #[derive(Debug, Clone, Hash, PartialEq, Eq)]
    pub struct In {
        start: usize,
        length: usize,
    }

    #[allow(dead_code)]
    impl In {
        pub fn start(&self) -> usize {
            self.start
        }
        pub fn length(&self) -> usize {
            self.length
        }
        pub fn end(&self) -> usize {
            self.start + self.length
        }
    }

    // Iterator over the Ex introns from lowest to highest coordinate
    pub struct Ins<'a> {
        curr_start: usize,
        rest: Iter<'a, InEx>,
    }

    impl<'a> Ins<'a> {
        #[allow(dead_code)]
        pub fn new(exon_0_length: usize, inexes: &'a Vec<InEx>) -> Self {
            Ins {
                curr_start: exon_0_length,
                rest: inexes.iter(),
            }
        }
    }

    impl<'a> Iterator for Ins<'a> {
        type Item = In;

        fn next(&mut self) -> Option<In> {
            match self.rest.next() {
                Some(inex) => {
                    let intr = In {
                        start: self.curr_start,
                        length: inex.intron_length(),
                    };
                    self.curr_start += inex.length();
                    Some(intr)
                }
                None => None,
            }
        }
    }
}

/// Spliced sequence annotation on a particular, named sequence
/// (e.g. a chromosome).
///
/// Parameterized over the type of the reference sequence identifier
/// and over the strandedness of the position.
///
/// The display format for a `Spliced` is
/// _chr:start_0-end_0;start_1-end_1;...;start_N-end_N(+/-/.)_. The
/// boundaries for each individual exon are given as a half-open
/// 0-based interval, like the Rust `Range` and BED format.
///
/// ```
/// # use bio_types::strand::ReqStrand;
/// # use bio_types::annot::AnnotError;
/// # use bio_types::annot::spliced::{Spliced,SplicingError};
/// # fn try_main() -> Result<(), Box<SplicingError>> {
/// let tad3 = Spliced::with_lengths_starts("chrXII".to_owned(), 765265,
///                                         &vec![808,52,109], &vec![0,864,984],
///                                         ReqStrand::Reverse)?;
/// assert_eq!(tad3.to_string(), "chrXII:765265-766073;766129-766181;766249-766358(-)");
/// let tad3_exons = tad3.exon_contigs();
/// assert_eq!(tad3_exons.len(), 3);
/// assert_eq!(tad3_exons[0].to_string(), "chrXII:766249-766358(-)");
/// assert_eq!(tad3_exons[1].to_string(), "chrXII:766129-766181(-)");
/// assert_eq!(tad3_exons[2].to_string(), "chrXII:765265-766073(-)");
/// # Ok(())
/// # }
/// # fn main() { try_main().unwrap(); }
/// ```
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct Spliced<R, S> {
    refid: R,
    start: isize,
    exon_0_length: usize,
    inexes: Vec<inex::InEx>,
    strand: S,
}

impl<R, S> Spliced<R, S> {
    /// Construct a new, single-exon "spliced" location
    ///
    /// ```
    /// use std::rc::Rc;
    /// use bio_types::annot::spliced::Spliced;
    /// use bio_types::strand::ReqStrand;
    /// let chr = Rc::new("chrX".to_owned());
    /// let tma22 = Spliced::new(chr, 461829, 462426 - 461829, ReqStrand::Forward);
    /// ```
    pub fn new(refid: R, start: isize, exon_0_length: usize, strand: S) -> Self {
        Spliced {
            refid: refid,
            start: start,
            exon_0_length: exon_0_length,
            inexes: Vec::new(),
            strand: strand,
        }
    }

    /// Construct a multi-exon "spliced" location using BED-style exon
    /// starts and lengths.
    pub fn with_lengths_starts(
        refid: R,
        start: isize,
        exon_lengths: &[usize],
        exon_starts: &[usize],
        strand: S,
    ) -> Result<Self, SplicingError> {
        if exon_starts.len() < 1 {
            return Err(SplicingError::NoExons);
        } else if exon_starts[0] != 0 {
            return Err(SplicingError::BlockStart);
        } else if exon_starts.len() != exon_lengths.len() {
            return Err(SplicingError::BlockMismatch);
        }

        let exon_0_length = exon_lengths[0];
        let mut intron_start = exon_0_length;
        let mut inexes = Vec::new();
        for exno in 1..exon_starts.len() {
            let exon_start = exon_starts[exno];
            if intron_start >= exon_start {
                return Err(SplicingError::BlockOverlap);
            }
            let intron_length = exon_start - intron_start;
            let exon_length = exon_lengths[exno];
            if exon_length == 0 {
                return Err(SplicingError::ExonLength);
            }
            inexes.push(inex::InEx::new(intron_length, exon_length)?);
            intron_start = exon_start + exon_length;
        }

        Ok(Spliced {
            refid: refid,
            start: start,
            exon_0_length: exon_0_length,
            inexes: inexes,
            strand: strand,
        })
    }

    /// Number of exons
    pub fn exon_count(&self) -> usize {
        self.inexes.len() + 1
    }

    // This is the unique representation for zero-length Spliced
    // locations, because InEx pairs have positive lengths for both
    // introns and exons.
    #[allow(dead_code)]
    fn is_zero_length(&self) -> bool {
        self.exon_0_length == 0 && self.inexes.is_empty()
    }

    fn exes<'a>(&'a self) -> inex::Exes<'a> {
        inex::Exes::new(self.exon_0_length, &self.inexes)
    }

    /// Vector of exon starting positions, relative to the start of
    /// the location overall.
    ///
    /// These positions run from left to right on the reference
    /// sequence, regardless of the location's strand.
    pub fn exon_starts(&self) -> Vec<usize> {
        let mut starts = vec![0];
        let mut intron_start = self.exon_0_length;
        for inex in self.inexes.iter() {
            starts.push(intron_start + inex.intron_length());
            intron_start += inex.length();
        }
        starts
    }

    /// Vector of exon lengths.
    ///
    /// Exon lengths are given from left to right on the reference
    /// sequence, regardless of the location's strand.
    pub fn exon_lengths(&self) -> Vec<usize> {
        let mut lengths = vec![self.exon_0_length];
        for inex in self.inexes.iter() {
            lengths.push(inex.exon_length());
        }
        lengths
    }

    /// Total length of exons only.
    ///
    /// The `length` method from the `Loc` trait returns the total
    /// length spanned by the annotation, including both introns and
    /// exons.
    pub fn exon_total_length(&self) -> usize {
        self.exes().map(|e| e.length()).sum()
    }

    /// Convert into an unstranded sequence location
    pub fn into_unstranded(self) -> Spliced<R, NoStrand> {
        Spliced {
            refid: self.refid,
            start: self.start,
            exon_0_length: self.exon_0_length,
            inexes: self.inexes,
            strand: NoStrand::Unknown,
        }
    }

    /// Convert into a stranded sequence location on the specified strand
    pub fn into_stranded(self, strand: ReqStrand) -> Spliced<R, ReqStrand> {
        Spliced {
            refid: self.refid,
            start: self.start,
            exon_0_length: self.exon_0_length,
            inexes: self.inexes,
            strand: strand,
        }
    }

    pub fn contig_cover(self) -> Contig<R, S> {
        let length = self.length();
        Contig::new(self.refid, self.start, length, self.strand)
    }

    // Get exons in reference sequence order.
    fn exon_contigs_vec(&self) -> Vec<Contig<R, S>>
    where
        R: Clone,
        S: Copy,
    {
        let mut exons = Vec::new();

        for ex in self.exes() {
            exons.push(Contig::new(
                self.refid().clone(),
                self.start + ex.start() as isize,
                ex.length(),
                self.strand,
            ));
        }

        exons
    }
}

impl<R> Spliced<R, ReqStrand> {
    /// Extend the annotation by `dist` in the upstream direction on the
    /// annotated strand.
    ///
    /// # Arguments
    ///
    /// * `dist` specifies the offset for sliding the position. The
    /// left, 5'-most end of the spliced will expand for forward-strand
    /// annotations and the right, 3'-most end will expand for
    /// reverse-strand annotations.
    // pub fn extend_upstream(&mut self, dist: usize) {
    //     self.length += dist;
    //     if self.strand == ReqStrand::Forward {
    //         self.start -= dist as isize;
    //     }
    // }

    /// Extend the annotation by `dist` in the downstream direction on the
    /// annotated strand.
    ///
    /// # Arguments
    ///
    /// * `dist` specifies the offset for sliding the position. The
    /// right, 3'-most end of the spliced will expand for
    /// forward-strand annotations and the left, 5'-most end will
    /// expand for reverse-strand annotations.
    // pub fn extend_downstream(&mut self, dist: usize) {
    //     self.length += dist;
    //     if self.strand == ReqStrand::Reverse {
    //         self.start -= dist as isize;
    //     }
    // }

    pub fn try_from<S>(x: Spliced<R, S>) -> Result<Self, AnnotError>
    where
        S: Strandedness,
    {
        match x.strand.try_req_strand() {
            None => Err(AnnotError::NoStrand),
            Some(strand) => Ok(Spliced {
                refid: x.refid,
                start: x.start,
                exon_0_length: x.exon_0_length,
                inexes: x.inexes,
                strand: strand,
            }),
        }
    }

    pub fn exon_contigs(&self) -> Vec<Contig<R, ReqStrand>>
    where
        R: Clone,
    {
        let mut exons = self.exon_contigs_vec();
        if self.strand == ReqStrand::Reverse {
            exons.reverse()
        }
        exons
    }
}

impl<R, S> Loc for Spliced<R, S> {
    type RefID = R;
    type Strand = S;
    fn refid(&self) -> &R {
        &self.refid
    }
    fn start(&self) -> isize {
        self.start
    }
    fn length(&self) -> usize {
        let mut len = self.exon_0_length;
        for inex in self.inexes.iter() {
            len += inex.length()
        }
        len
    }
    fn strand(&self) -> S
    where
        S: Copy,
    {
        self.strand
    }

    fn pos_into<T>(&self, pos: &Pos<Self::RefID, T>) -> Option<Pos<(), T>>
    where
        Self::RefID: Eq,
        Self::Strand: KnownStrandedness + Copy,
        T: Strandedness + Copy,
    {
        if self.refid != *pos.refid() {
            return None;
        } else if pos.pos() < self.start {
            return None;
        }

        let pos_offset = (pos.pos() - self.start) as usize;

        let mut offset_before = 0;
        for ex in self.exes() {
            if pos_offset >= ex.start() && pos_offset < ex.end() {
                let offset = (offset_before + pos_offset - ex.start()) as isize;
                let into = match self.strand().req_strand() {
                    ReqStrand::Forward => Pos::new((), offset, self.strand().product(pos.strand())),
                    ReqStrand::Reverse => Pos::new(
                        (),
                        self.exon_total_length() as isize - (offset + 1),
                        self.strand().product(pos.strand()),
                    ),
                };
                return Some(into);
            }
            offset_before += ex.length();
        }

        None
    }

    fn pos_outof<Q, T>(&self, pos: &Pos<Q, T>) -> Option<Pos<Self::RefID, T>>
    where
        Self::RefID: Clone,
        Self::Strand: KnownStrandedness + Copy,
        T: Strandedness + Copy,
    {
        let mut offset = match self.strand().req_strand() {
            ReqStrand::Forward => pos.pos(),
            ReqStrand::Reverse => self.exon_total_length() as isize - (pos.pos() + 1),
        };

        if offset < 0 {
            return None;
        }

        for ex in self.exes() {
            if offset < ex.length() as isize {
                return Some(Pos::new(
                    self.refid.clone(),
                    self.start + ex.start() as isize + offset,
                    self.strand.product(pos.strand()),
                ));
            }
            offset -= ex.length() as isize;
        }

        None
    }

    fn contig_intersection<T>(&self, contig: &Contig<Self::RefID, T>) -> Option<Self>
    where
        Self::RefID: PartialEq + Clone,
        Self::Strand: Copy,
    {
        if self.refid() != contig.refid() {
            return None;
        }

        let contig_rel_start = if contig.start() < self.start {
            0
        } else {
            (contig.start() - self.start) as usize
        };
        let contig_end = contig.start() + contig.length() as isize;
        let contig_rel_end = if contig_end < self.start {
            0
        } else {
            (contig_end - self.start) as usize
        };

        let mut exon_lengths = Vec::new();
        let mut exon_starts = Vec::new();

        for ex in self.exes() {
            let start = max(contig_rel_start, ex.start());
            let end = min(contig_rel_end, ex.end());

            if start < end {
                exon_starts.push(start - contig_rel_start);
                exon_lengths.push(end - start);
            }
        }

        if exon_starts.len() > 0 {
            let first_start = exon_starts[0];
            for start in exon_starts.iter_mut() {
                *start -= first_start;
            }
            let ixn = Self::with_lengths_starts(
                self.refid.clone(),
                max(self.start, contig.start()) + first_start as isize,
                &exon_lengths,
                &exon_starts,
                self.strand,
            ).unwrap_or_else(|e| {
                panic!(
                    "Creating intersection spliced: {:?} for {:?} {:?}",
                    e, exon_lengths, exon_starts
                )
            });

            Some(ixn)
        } else {
            None
        }
    }
}

impl<R, S> Display for Spliced<R, S>
where
    R: Display,
    S: Display,
{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(f, "{}:", self.refid)?;

        let mut sep = false;
        for ex in self.exes() {
            let ex_start = self.start + ex.start() as isize;
            write!(
                f,
                "{}{}-{}",
                if sep { ";" } else { "" },
                ex_start,
                ex_start + ex.length() as isize
            )?;
            sep = true;
        }
        write!(f, "{}", self.strand)
    }
}

// impl <R,S> FromStr for Contig<R,S>
//     where R: From<String>, S: Strandedness
// {
//     type Err = ParseError;

//     fn from_str(s: &str) -> Result<Self, Self::Err> {
//         let (refidstr, rest) = break_refid(s)?;
//         let (posstr, strand) = S::break_pos_strand(rest)
//             .ok_or_else(|| ParseError::BadStrand(s.to_owned()))?;
//         let (start, end) = break_start_end(posstr)?;
//         if start <= end {
//             Ok( Contig::new(R::from(refidstr.to_owned()), start, (end - start) as usize, strand) )
//         } else {
//             Err( ParseError::EndBeforeStart(s.to_owned()) )
//         }
//     }
// }

impl<R> From<Spliced<R, ReqStrand>> for Spliced<R, Strand> {
    fn from(x: Spliced<R, ReqStrand>) -> Self {
        Spliced {
            refid: x.refid,
            start: x.start,
            exon_0_length: x.exon_0_length,
            inexes: x.inexes,
            strand: match x.strand {
                ReqStrand::Forward => Strand::Forward,
                ReqStrand::Reverse => Strand::Reverse,
            },
        }
    }
}

impl<R> From<Spliced<R, NoStrand>> for Spliced<R, Strand> {
    fn from(x: Spliced<R, NoStrand>) -> Self {
        Spliced {
            refid: x.refid,
            start: x.start,
            exon_0_length: x.exon_0_length,
            inexes: x.inexes,
            strand: Strand::Unknown,
        }
    }
}

/// Default stranded sequence position on a reference sequence named
/// by a `String`.
pub type SeqSplicedStranded = Spliced<String, ReqStrand>;

/// Default unstranded sequence position on a reference sequence named
/// by a `String`
pub type SeqSplicedUnstranded = Spliced<String, NoStrand>;

quick_error! {
    #[derive(Debug, Clone)]
    pub enum SplicingError {
        IntronLength {
            description("Invalid (non-positive) intron length")
        }
        ExonLength {
            description("Invalid (non-positive) exon length")
        }
        NoExons {
            description("No exons")
        }
        BlockStart {
            description("Exons do not start at position 0")
        }
        BlockMismatch {
            description("Number of exon starts != number of exon lengths")
        }
        BlockOverlap {
            description("Exon blocks overlap")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length_start_to_contig() {
        //chrV    166236  166885  YER007C-A       0       -       166236  166885  0       2       535,11, 0,638,
        let tma20 = Spliced::with_lengths_starts(
            "chrV".to_owned(),
            166236,
            &vec![535, 11],
            &vec![0, 638],
            ReqStrand::Reverse,
        ).unwrap();
        assert_eq!(tma20.exon_starts(), vec![0, 638]);
        assert_eq!(tma20.exon_lengths(), vec![535, 11]);
        assert_eq!(tma20.to_string(), "chrV:166236-166771;166874-166885(-)");
        let tma20_exons = tma20.exon_contigs();
        assert_eq!(tma20_exons.len(), 2);
        assert_eq!(tma20_exons[0].to_string(), "chrV:166874-166885(-)");
        assert_eq!(tma20_exons[1].to_string(), "chrV:166236-166771(-)");

        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        ).unwrap();
        assert_eq!(
            rpl7b.to_string(),
            "chrXVI:173151-173162;173571-173665;174072-174702(+)"
        );
        let rpl7b_exons = rpl7b.exon_contigs();
        assert_eq!(rpl7b_exons.len(), 3);
        assert_eq!(rpl7b_exons[0].to_string(), "chrXVI:173151-173162(+)");
        assert_eq!(rpl7b_exons[1].to_string(), "chrXVI:173571-173665(+)");
        assert_eq!(rpl7b_exons[2].to_string(), "chrXVI:174072-174702(+)");

        //chrXII  765265  766358  YLR316C 0       -       765265  766358  0       3       808,52,109,     0,864,984,
        let tad3 = Spliced::with_lengths_starts(
            "chrXII".to_owned(),
            765265,
            &vec![808, 52, 109],
            &vec![0, 864, 984],
            ReqStrand::Reverse,
        ).unwrap();
        assert_eq!(
            tad3.to_string(),
            "chrXII:765265-766073;766129-766181;766249-766358(-)"
        );
        let tad3_exons = tad3.exon_contigs();
        assert_eq!(tad3_exons.len(), 3);
        assert_eq!(tad3_exons[0].to_string(), "chrXII:766249-766358(-)");
        assert_eq!(tad3_exons[1].to_string(), "chrXII:766129-766181(-)");
        assert_eq!(tad3_exons[2].to_string(), "chrXII:765265-766073(-)");
    }

    fn test_into_outof(
        loc: &Spliced<String, ReqStrand>,
        outstr: &str,
        in_offset: isize,
        in_strand: ReqStrand,
    ) -> () {
        let p0 = outstr.parse::<Pos<String, ReqStrand>>().unwrap();
        let p0_into_expected = Pos::new((), in_offset, in_strand);
        let p0_into_actual = loc.pos_into(&p0);
        let p0_back_out_actual = loc.pos_outof(&p0_into_expected);
        print!(
            "{}\t{}\t{:?}\t{:?}\t{:?}\n",
            outstr, p0, p0_into_expected, p0_into_actual, p0_back_out_actual
        );
        assert!(Some(p0_into_expected).same(&p0_into_actual));
        assert!(Some(p0).same(&p0_back_out_actual));
    }

    fn test_no_into(loc: &Spliced<String, ReqStrand>, outstr: &str) -> () {
        let p0 = outstr.parse::<Pos<String, ReqStrand>>().unwrap();
        assert!(None.same(&loc.pos_into(&p0)));
    }

    #[test]
    fn into_outof() {
        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        ).unwrap();
        let p0_into = Pos::new((), -1, ReqStrand::Forward);
        assert!(None.same(&rpl7b.pos_outof(&p0_into)));

        test_no_into(&rpl7b, "chrXVI:173150(+)");
        test_into_outof(&rpl7b, "chrXVI:173151(+)", 0, ReqStrand::Forward);
        test_into_outof(&rpl7b, "chrXVI:173152(-)", 1, ReqStrand::Reverse);
        test_into_outof(&rpl7b, "chrXVI:173161(+)", 10, ReqStrand::Forward);
        test_no_into(&rpl7b, "chrXVI:173162(+)");
        test_no_into(&rpl7b, "chrXVI:173570(+)");
        test_into_outof(&rpl7b, "chrXVI:173571(+)", 11, ReqStrand::Forward);
        test_into_outof(&rpl7b, "chrXVI:173664(+)", 104, ReqStrand::Forward);
        test_no_into(&rpl7b, "chrXVI:173665(+)");
        test_no_into(&rpl7b, "chrXVI:174071(+)");
        test_into_outof(&rpl7b, "chrXVI:174072(+)", 105, ReqStrand::Forward);
        test_into_outof(&rpl7b, "chrXVI:174701(+)", 734, ReqStrand::Forward);
        test_no_into(&rpl7b, "chrXVI:174702(+)");

        let p0_into = Pos::new((), 735, ReqStrand::Forward);
        assert!(None.same(&rpl7b.pos_outof(&p0_into)));

        //chrXII  765265  766358  YLR316C 0       -       765265  766358  0       3       808,52,109,     0,864,984,
        let tad3 = Spliced::with_lengths_starts(
            "chrXII".to_owned(),
            765265,
            &vec![808, 52, 109],
            &vec![0, 864, 984],
            ReqStrand::Reverse,
        ).unwrap();

        let p0_into = Pos::new((), -1, ReqStrand::Forward);
        assert!(None.same(&tad3.pos_outof(&p0_into)));

        test_no_into(&tad3, "chrXII:765264(-)");
        test_into_outof(&tad3, "chrXII:765265(-)", 968, ReqStrand::Forward);
        test_into_outof(&tad3, "chrXII:765266(+)", 967, ReqStrand::Reverse);
        test_into_outof(&tad3, "chrXII:766072(-)", 161, ReqStrand::Forward);
        test_no_into(&tad3, "chrXII:766073(-)");

        test_no_into(&tad3, "chrXII:766128(-)");
        test_into_outof(&tad3, "chrXII:766129(-)", 160, ReqStrand::Forward);
        test_into_outof(&tad3, "chrXII:766180(-)", 109, ReqStrand::Forward);
        test_no_into(&tad3, "chrXII:766181(-)");

        test_no_into(&tad3, "chrXII:766248(-)");
        test_into_outof(&tad3, "chrXII:766249(-)", 108, ReqStrand::Forward);
        test_into_outof(&tad3, "chrXII:766357(-)", 0, ReqStrand::Forward);
        test_no_into(&tad3, "chrXII:766358(-)");

        let p0_into = Pos::new((), 969, ReqStrand::Forward);
        assert!(None.same(&tad3.pos_outof(&p0_into)));
    }

    fn test_contig_ixn(
        spl: &Spliced<String, ReqStrand>,
        cb_str: &str,
        cab_str: Option<String>,
    ) -> () {
        let cb = cb_str.parse::<Contig<String, ReqStrand>>().unwrap();
        match spl.contig_intersection(&cb) {
            None => assert_eq!(None, cab_str),
            Some(cab) => assert_eq!(Some(cab.to_string()), cab_str),
        };
    }

    #[test]
    fn intersection() {
        //chrXVI  173151  174702  YPL198W 0       +       173151  174702  0       3       11,94,630,      0,420,921,
        let rpl7b = Spliced::with_lengths_starts(
            "chrXVI".to_owned(),
            173151,
            &vec![11, 94, 630],
            &vec![0, 420, 921],
            ReqStrand::Forward,
        ).unwrap();

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-175000(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-175000(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173151-175000(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173152-175000(+)",
            Some("chrXVI:173152-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-175000(+)",
            Some("chrXVI:173155-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173161-175000(+)",
            Some("chrXVI:173161-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173162-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173500-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173570-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173571-175000(+)",
            Some("chrXVI:173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173572-175000(+)",
            Some("chrXVI:173572-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173600-175000(+)",
            Some("chrXVI:173600-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173664-175000(+)",
            Some("chrXVI:173664-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173665-175000(+)",
            Some("chrXVI:174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:174100-175000(+)",
            Some("chrXVI:174100-174702(+)".to_owned()),
        );
        test_contig_ixn(&rpl7b, "chrXVI:174800-175000(+)", None);

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-174703(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-174702(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174702(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173150-174701(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174701(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-174500(+)",
            Some("chrXVI:173151-173162;173571-173665;174072-174500(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-174072(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173800(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173666(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173665(+)",
            Some("chrXVI:173151-173162;173571-173665(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173664(+)",
            Some("chrXVI:173151-173162;173571-173664(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173600(+)",
            Some("chrXVI:173151-173162;173571-173600(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173571(+)",
            Some("chrXVI:173151-173162(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173300(+)",
            Some("chrXVI:173151-173162(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173000-173155(+)",
            Some("chrXVI:173151-173155(+)".to_owned()),
        );
        test_contig_ixn(&rpl7b, "chrXVI:173000-173100(+)", None);

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-174500(+)",
            Some("chrXVI:173155-173162;173571-173665;174072-174500(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173600-174500(+)",
            Some("chrXVI:173600-173665;174072-174500(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-173600(+)",
            Some("chrXVI:173155-173162;173571-173600(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:173590-173610(+)",
            Some("chrXVI:173590-173610(+)".to_owned()),
        );

        test_contig_ixn(
            &rpl7b,
            "chrXVI:173155-173160(+)",
            Some("chrXVI:173155-173160(+)".to_owned()),
        );
        test_contig_ixn(
            &rpl7b,
            "chrXVI:174400-174500(+)",
            Some("chrXVI:174400-174500(+)".to_owned()),
        );

        test_contig_ixn(&rpl7b, "chrXVI:173200-173300(+)", None);
        test_contig_ixn(&rpl7b, "chrXVI:173800-174000(+)", None);
    }

}
