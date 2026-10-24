#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Node<const D: usize> {
    pub id: usize,
    pub coordinates: [f64; D],
}

impl<const D: usize> Node<D> {
    pub fn new(id: usize, coordinates: [f64; D]) -> Self {
        Node { id, coordinates }
    }
}

