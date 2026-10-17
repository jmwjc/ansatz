pub mod bar;
pub mod elasticity;
pub mod heat;

use std::collections::BTreeSet;

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

/// Number of independent components of a fully symmetric tensor of order `O`
/// in `D` spatial dimensions:
/// `C(D + O - 1, O) = (D + O - 1)! / ((D - 1)! · O!)`.
///
/// This yields `1` for `O = 0`, `D` for `O = 1`, and `D(D+1)/2` for `O = 2`,
/// and generalizes to any tensor order. It is evaluated as the binomial
/// product `C(n, k)` with `n = D + O - 1` and `k = min(O, D - 1)`, keeping the
/// intermediate values bounded instead of computing the factorials first.
const fn symmetric_components<const D: usize>(order: usize) -> usize {
    let n = D + order - 1;
    let k = if order < D - 1 { order } else { D - 1 };
    let mut c = 1usize;
    let mut i = 0usize;
    while i < k {
        c = c * (n - k + 1 + i) / (i + 1);
        i += 1;
    }
    c
}

/// A collection of elements over which an unknown (variable) field is defined.
///
/// The const generic parameters are:
/// - `D`: spatial dimension of the element;
/// - `O`: tensor order of the variable (0 = scalar, 1 = vector, 2 = second-order
///   tensor, …);
/// - `P`: parametric dimension of the element;
/// - `N`: number of nodes per element.
pub trait Variable<const D: usize, const O: usize, const P: usize, const N: usize> {
    type Item: Element<D, P, N>;
    fn iter(&self) -> impl Iterator<Item = &Self::Item>;

    /// Number of global scalar degrees of freedom (DOFs) over which this
    /// variable is defined: the number of unique node ids appearing as element
    /// vertices, multiplied by the number of independent components per node of
    /// a fully symmetric tensor of order `O`
    /// ([`symmetric_components`] = `C(D + O - 1, O)`).
    fn num_dofs(&self) -> usize {
        self.iter()
            .flat_map(|elm| elm.id())
            .collect::<BTreeSet<_>>()
            .len()
            * const { symmetric_components::<D>(O) }
    }

    /// Evaluates the field represented by this variable at the physical point
    /// `x`, by locating the element containing `x` ([`Element::contains`]) and
    /// interpolating its nodal coefficients through the shape functions.
    ///
    /// `coefficients` is the global coefficient vector, of length
    /// `self.num_dofs()`; node `n` owns `m` scalar coefficients
    /// `coefficients[n * m .. (n + 1) * m]` (node ids are expected to be
    /// 0-based DOF indices), where
    /// `m = symmetric_components::<D>(O) = C(D + O - 1, O)` is the number of
    /// independent components of a fully symmetric tensor of order `O`. The
    /// result is a `Vec<f64>` of length `m`.
    ///
    /// Returns `None` when `x` lies outside every element; on a shared boundary
    /// the first element reporting containment wins (the interpolation is
    /// continuous there, so the result is independent of that choice).
    fn value(&self, x: [f64; D], coefficients: &[f64]) -> Option<Vec<f64>> {
        for elm in self.iter() {
            if let Some(ξ) = elm.local_coordinates(x) {
                let shape = elm.shape(ξ);
                let id = elm.id();
                let m = const { symmetric_components::<D>(O) };
                let mut acc = vec![0.0; m];
                for i in 0..N {
                    let start = id[i] * m;
                    for k in 0..m {
                        acc[k] += shape[i] * coefficients[start + k];
                    }
                }
                return Some(acc);
            }
        }
        None
    }
}
impl<T, const D: usize, const O: usize, const P: usize, const N: usize> Variable<D, O, P, N> for &[T]
where
    T: Element<D, P, N>,
{
    type Item = T;

    fn iter(&self) -> impl Iterator<Item = &T> {
        <[T]>::iter(self)
    }
}
impl<T, const D: usize, const O: usize, const P: usize, const N: usize> Variable<D, O, P, N> for Vec<T>
where
    T: Element<D, P, N>,
{
    type Item = T;
    fn iter(&self) -> impl Iterator<Item = &T> {
        self.as_slice().iter()
    }
}

#[cfg(test)]
mod tests {
    use super::Variable;
    use crate::element::seg2::Seg2;
    use crate::node::Node;

    #[test]
    fn num_dofs_counts_unique_node_ids() {
        // Three Seg2 elements sharing interior nodes: node ids 0,1,2,3 appear,
        // with 1 and 2 repeated across elements. Unique DOFs = 4.
        let bars: Vec<Seg2<1>> = vec![
            Seg2::new(Node::new(0, [0.0]), Node::new(1, [1.0])),
            Seg2::new(Node::new(1, [1.0]), Node::new(2, [2.0])),
            Seg2::new(Node::new(2, [2.0]), Node::new(3, [3.0])),
        ];

        // Vec<T> and &[T] blanket impls must agree.
        assert_eq!(Variable::<1, 0, 1, 2>::num_dofs(&bars), 4);
        assert_eq!(Variable::<1, 0, 1, 2>::num_dofs(&bars.as_slice()), 4);
    }

    #[test]
    fn num_dofs_empty_variable_is_zero() {
        let bars: Vec<Seg2<1>> = Vec::new();
        assert_eq!(Variable::<1, 0, 1, 2>::num_dofs(&bars), 0);
    }

    #[test]
    fn num_dofs_scales_with_tensor_order() {
        // One Seg2<2> element: two unique nodes.
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
        )];
        assert_eq!(Variable::<2, 0, 1, 2>::num_dofs(&bars), 2); // scalar: 2 * 1
        assert_eq!(Variable::<2, 1, 1, 2>::num_dofs(&bars), 4); // vector: 2 * 2
        assert_eq!(Variable::<2, 2, 1, 2>::num_dofs(&bars), 6); // symmetric 2nd-order: 2 * 3
    }

    #[test]
    fn value_interpolates_at_physical_point() {
        // Two Seg2 elements over [0, 4]; nodal coefficients 0, 2, 4 give u(x) = x.
        let bars: Vec<Seg2<1>> = vec![
            Seg2::new(Node::new(0, [0.0]), Node::new(1, [2.0])),
            Seg2::new(Node::new(1, [2.0]), Node::new(2, [4.0])),
        ];
        let coefficients = [0.0, 2.0, 4.0];

        assert_eq!(
            Variable::<1, 0, 1, 2>::value(&bars, [1.0], &coefficients),
            Some(vec![1.0])
        );
        assert_eq!(
            Variable::<1, 0, 1, 2>::value(&bars, [3.0], &coefficients),
            Some(vec![3.0])
        );
        assert_eq!(
            Variable::<1, 0, 1, 2>::value(&bars, [5.0], &coefficients),
            None
        );
    }

    #[test]
    fn value_returns_vector_for_order_one() {
        // One Seg2<2> element; node 0 coefficient vector (1, 2), node 1 (3, 0).
        // At the midpoint the interpolated vector is the average (2, 1).
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [2.0, 0.0]),
        )];
        let coefficients = [1.0, 2.0, 3.0, 0.0];

        assert_eq!(
            Variable::<2, 1, 1, 2>::value(&bars, [1.0, 0.0], &coefficients),
            Some(vec![2.0, 1.0])
        );
    }

    #[test]
    fn value_returns_second_order_tensor() {
        // One Seg2<2> element; node 0 Voigt tensor [σ11, σ12, σ22] = [1, 2, 3],
        // node 1 = [4, 5, 6]. At the midpoint this averages to [2.5, 3.5, 4.5].
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [2.0, 0.0]),
        )];
        let coefficients = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

        assert_eq!(
            Variable::<2, 2, 1, 2>::value(&bars, [1.0, 0.0], &coefficients),
            Some(vec![2.5, 3.5, 4.5])
        );
    }
}
