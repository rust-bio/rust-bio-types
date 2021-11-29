//! A distance matrix is stored as a list of taxon names, together with a square
//! or triangular matrix representing all pairwise distances.
//!
//! # Example
//!
//! ```rust
//! use bio_types::distancematrix::DistanceMatrix;
//!
//! DistanceMatrix::new(vec!["a".to_string(), "b".to_string()],
//! vec![
//!  vec![0.0, 1.0],
//!  vec![1.0, 0.0],
//! ]);
//! ```

use std::{
    cmp::{max, min},
    ops::{Index, IndexMut},
};

use thiserror::Error;

/// The type of the matrix. Either square, lower triangular (excluding the
/// diagonal), or upper triangular (excluding the diagonal).
/// Square matrices are assumed to be symmetric, but this is not enforced nor checked.
#[derive(Debug, PartialEq, Eq)]
pub enum MatrixType {
    Square,
    Lower,
    Upper,
}

/// A distance matrix containing a list of taxon names and a matrix of pairwise distances.
#[derive(Debug)]
pub struct DistanceMatrix {
    pub names: Vec<String>,
    pub distances: Vec<Vec<f32>>,
    pub matrix_type: MatrixType,
}

#[derive(Error, Debug)]
pub enum DistanceMatrixError {
    #[error("names and matrix do not have matching length")]
    LengthError,
    #[error("matrix has unrecognized shape")]
    ShapeError,
}

impl DistanceMatrix {
    /// The number of sequences (rows) in the matrix.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Create a new DistanceMatrix.
    ///
    /// `distances` must be a matrix with the same number of rows as `names`, and follow one of the accepted shapes.
    pub fn new(names: Vec<String>, distances: Vec<Vec<f32>>) -> Result<Self, DistanceMatrixError> {
        if names.len() != distances.len() {
            Err(DistanceMatrixError::LengthError)
        } else {
            let n = names.len();
            let square = distances.iter().all(|x| x.len() == n);
            let lower = distances.iter().enumerate().all(|(i, x)| x.len() == i);
            let upper = distances
                .iter()
                .enumerate()
                .all(|(i, x)| x.len() == n - 1 - i);
            let matrix_type = match (square, lower, upper) {
                (true, _, _) => MatrixType::Square,
                (_, true, _) => MatrixType::Lower,
                (_, _, true) => MatrixType::Upper,
                _ => Err(DistanceMatrixError::ShapeError)?,
            };
            Ok(Self {
                names,
                distances,
                matrix_type,
            })
        }
    }

    // Converts an index (i,j) to the right index to use given the shape of the matrix.
    fn get_index(&self, (i, j): (usize, usize)) -> (usize, usize) {
        match self.matrix_type {
            MatrixType::Square => (i, j),
            MatrixType::Lower => {
                assert!(i != j);
                (max(i, j), min(i, j))
            }
            MatrixType::Upper => {
                assert!(i != j);
                (min(i, j), max(i, j) - min(i, j) - 1)
            }
        }
    }
}

/// Index access into the DistanceMatrix, taking into account the shape.
///
/// Indices should be used as if the matrix was square. For triangular matrices,
/// `i` and `j` must be distinct, and `(i,j)` and `(j,i)` represent the same
/// element.
///
/// # Example
///
/// ```
/// use bio_types::distancematrix::DistanceMatrix;
///
/// let mut t = DistanceMatrix::new(vec!["a".to_string(), "b".to_string()],
/// vec![
///  vec![],
///  vec![1.0],
/// ]).unwrap();
/// assert_eq!(t[(1, 0)], 1.0);
/// t[(0, 1)] = 2.0;
/// assert_eq!(t[(1, 0)], 2.0);
/// ```
impl Index<(usize, usize)> for DistanceMatrix {
    type Output = f32;
    fn index(&self, t: (usize, usize)) -> &Self::Output {
        let (i, j) = self.get_index(t);
        &self.distances[i][j]
    }
}

/// Mutable index access into the DistanceMatrix, taking into account the shape, like `Index`.
impl IndexMut<(usize, usize)> for DistanceMatrix {
    fn index_mut(&mut self, t: (usize, usize)) -> &mut Self::Output {
        let (i, j) = self.get_index(t);
        &mut self.distances[i][j]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexing() {
        let mut m = DistanceMatrix::new(
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
            vec![vec![0., 1., 2.], vec![3., 4., 5.], vec![6., 7., 8.]],
        )
        .unwrap();
        assert_eq!(m.matrix_type, MatrixType::Square);
        assert_eq!(m[(0, 0)], 0.);
        assert_eq!(m[(1, 2)], 5.);
        assert_eq!(m[(2, 1)], 7.);
        m[(1, 2)] = -1.;
        m[(2, 1)] = -2.;
        assert_eq!(m[(1, 2)], -1.);
        assert_eq!(m[(2, 1)], -2.);
    }

    #[test]
    fn lower_indexing() {
        let mut m = DistanceMatrix::new(
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
            vec![vec![], vec![1.], vec![2., 3.]],
        )
        .unwrap();
        assert_eq!(m.matrix_type, MatrixType::Lower);
        assert_eq!(m[(2, 0)], 2.);
        assert_eq!(m[(0, 2)], 2.);
        assert_eq!(m[(1, 2)], 3.);
        assert_eq!(m[(2, 1)], 3.);
        m[(0, 2)] = -1.;
        m[(2, 1)] = -2.;
        m[(1, 2)] = -3.;
        assert_eq!(m[(0, 2)], -1.);
        assert_eq!(m[(2, 0)], -1.);
        assert_eq!(m[(1, 2)], -3.);
        assert_eq!(m[(2, 1)], -3.);
    }

    #[test]
    fn upper_indexing() {
        let mut m = DistanceMatrix::new(
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
            vec![vec![1., 2.], vec![3.], vec![]],
        )
        .unwrap();
        assert_eq!(m.matrix_type, MatrixType::Upper);
        assert_eq!(m[(0, 2)], 2.);
        assert_eq!(m[(1, 2)], 3.);
        assert_eq!(m[(2, 1)], 3.);
        m[(0, 2)] = -1.;
        m[(2, 1)] = -2.;
        m[(1, 2)] = -3.;
        assert_eq!(m[(0, 2)], -1.);
        assert_eq!(m[(2, 0)], -1.);
        assert_eq!(m[(1, 2)], -3.);
        assert_eq!(m[(2, 1)], -3.);
    }

    #[test]
    fn bad_length() {
        let m = DistanceMatrix::new(
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
            vec![vec![1., 2.], vec![3.]],
        );
        match m {
            Err(DistanceMatrixError::LengthError) => (),
            _ => assert!(false),
        };
    }

    #[test]
    fn bad_shape() {
        let m = DistanceMatrix::new(
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
            vec![vec![1., 2.], vec![3.], vec![4.]],
        );
        match m {
            Err(DistanceMatrixError::ShapeError) => (),
            _ => assert!(false),
        };
    }
}
