pub mod bar;
pub mod heat;

use crate::element::Element;

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
