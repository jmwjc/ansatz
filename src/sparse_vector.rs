//! Sparse vector assembly: `Doublet` entries and a constructor, analogous to the `Triplet` of faer's sparse matrices.

use std::fmt;

use faer::prelude::*;

/// An assembly entry for a sparse vector: a `(index, value)` pair, analogous to the `Triplet` of faer's sparse matrices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Doublet<I, T> {
    pub index: I,
    pub value: T,
}

impl<I, T> Doublet<I, T> {
    pub fn new(index: I, value: T) -> Self {
        Doublet { index, value }
    }
}

/// Errors produced while assembling a sparse vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DoubletError {
    /// A doublet references an index outside `0..len`.
    IndexOutOfRange { index: usize, len: usize },
}

impl fmt::Display for DoubletError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IndexOutOfRange { index, len } => {
                write!(f, "doublet index {index} out of range for length {len}")
            }
        }
    }
}

impl std::error::Error for DoubletError {}

/// Construct a dense vector from a list of `Doublet`s, accumulating values with duplicate indices.
///
/// `n` is the vector length (total number of DOFs). Accumulation is naturally conflict-free,
/// so this returns a dense vector directly (analogous to faer's
/// [`faer::sparse::SparseColMat::try_new_from_triplets`]).
///
/// # Errors
///
/// Returns [`DoubletError::IndexOutOfRange`] if a doublet index is not smaller
/// than `n`.
pub fn try_new_from_doublet(
    n: usize,
    doublets: &[Doublet<usize, f64>],
) -> Result<Col<f64>, DoubletError> {
    let mut v = Col::<f64>::zeros(n);
    for d in doublets {
        if d.index >= n {
            return Err(DoubletError::IndexOutOfRange {
                index: d.index,
                len: n,
            });
        }
        v[d.index] += d.value;
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accumulates_duplicate_indices() {
        let doublets = [
            Doublet::new(0, 1.0),
            Doublet::new(1, 2.0),
            Doublet::new(0, 3.0),
        ];
        let v = try_new_from_doublet(2, &doublets).unwrap();
        assert_eq!(v[0], 4.0);
        assert_eq!(v[1], 2.0);
    }

    #[test]
    fn rejects_out_of_range_index() {
        let err = try_new_from_doublet(2, &[Doublet::new(2, 1.0)]).unwrap_err();
        assert_eq!(err, DoubletError::IndexOutOfRange { index: 2, len: 2 });
    }
}
