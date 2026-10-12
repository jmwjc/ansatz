//! Sparse vector assembly: `Doublet` entries and a constructor, analogous to the `Triplet` of faer's sparse matrices.

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

/// Construct a dense vector from a list of `Doublet`s, accumulating values with duplicate indices.
///
/// `n` is the vector length (total number of DOFs). Accumulation is naturally conflict-free,
/// so this returns a dense vector directly (analogous to faer's
/// [`faer::sparse::SparseColMat::try_new_from_triplets`]).
pub fn try_new_from_doublet(n: usize, doublets: &[Doublet<usize, f64>]) -> Col<f64> {
    let mut v = Col::<f64>::zeros(n);
    for d in doublets {
        v[d.index] += d.value;
    }
    v
}
