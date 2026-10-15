pub mod bar;
pub mod heat;

use crate::element::Element;
use crate::sparse_vector::Doublet;
use faer::sparse::Triplet;

/// Sparse-matrix assembly triplets (entries of `K`).
pub type Triplets = Vec<Triplet<usize, usize, f64>>;

/// Sparse-vector assembly doublets (entries of `F`).
pub type Doublets = Vec<Doublet<usize, f64>>;

/// A stiffness-triplets / load-doublets pair, returned by operators that
/// assemble `K` and `F` together.
pub type TripletsAndDoublets = (Triplets, Doublets);

pub trait Variable<const D: usize, const P: usize, const N: usize> {
    type Item: Element<D, P, N>;
    fn iter(&self) -> impl Iterator<Item = &Self::Item>;
}
impl<T, const D: usize, const P: usize, const N: usize> Variable<D, P, N> for &[T]
where
    T: Element<D, P, N>,
{
    type Item = T;

    fn iter(&self) -> impl Iterator<Item = &T> {
        <[T]>::iter(self)
    }
}
impl<T, const D: usize, const P: usize, const N: usize> Variable<D, P, N> for Vec<T>
where
    T: Element<D, P, N>,
{
    type Item = T;
    fn iter(&self) -> impl Iterator<Item = &T> {
        self.as_slice().iter()
    }
}
