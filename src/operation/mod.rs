pub mod bar;
pub mod elasticity;
pub mod heat;

use std::collections::BTreeSet;
use std::fmt;

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

fn node_coefficients(
    coefficients: &[f64],
    node_id: usize,
    number_of_components: usize,
) -> Result<&[f64], VariableError> {
    let start = node_id
        .checked_mul(number_of_components)
        .ok_or(VariableError::DofOverflow)?;
    let end = start
        .checked_add(number_of_components)
        .ok_or(VariableError::DofOverflow)?;
    coefficients
        .get(start..end)
        .ok_or(VariableError::CoefficientOutOfRange {
            node_id,
            start,
            len: coefficients.len(),
            number_of_components,
        })
}

/// Errors produced while evaluating a variable field, or while assembling an
/// operator over two element sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariableError {
    /// The coefficient vector is too short for a node's degree of freedom range.
    CoefficientOutOfRange {
        node_id: usize,
        start: usize,
        len: usize,
        number_of_components: usize,
    },
    /// Computing a node's coefficient offset overflowed `usize`.
    DofOverflow,
    /// Two element sets that must be paired element-by-element differ in length.
    ElementCountMismatch {
        /// Number of elements in the first set (e.g. the temperature elements).
        first: usize,
        /// Number of elements in the second set (e.g. the multiplier elements).
        second: usize,
    },
}

impl fmt::Display for VariableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CoefficientOutOfRange {
                node_id,
                start,
                len,
                number_of_components,
            } => write!(
                f,
                "node {node_id} needs coefficients [{start}..{}], but the vector has length {len}",
                start.saturating_add(*number_of_components)
            ),
            Self::DofOverflow => write!(f, "degree-of-freedom offset overflowed usize"),
            Self::ElementCountMismatch { first, second } => write!(
                f,
                "cannot pair {first} elements with {second} elements: both element sets must have the same length"
            ),
        }
    }
}

impl std::error::Error for VariableError {}

/// A collection of elements over which an unknown (variable) field is defined.
///
/// The const generic parameters are:
/// - `D`: spatial dimension of the element;
/// - `C`: number of scalar degrees of freedom per node;
/// - `P`: parametric dimension of the element;
/// - `N`: number of nodes per element.
pub trait Variable<const D: usize, const C: usize, const P: usize, const N: usize> {
    type Item: Element<D, P, N>;
    fn iter(&self) -> impl Iterator<Item = &Self::Item>;

    /// Number of global scalar degrees of freedom (DOFs): number of unique
    /// node ids multiplied by the number of components `C` per node.
    fn num_dofs(&self) -> usize {
        self.iter()
            .flat_map(|elm| elm.id())
            .collect::<BTreeSet<_>>()
            .len()
            * C
    }

    /// Locates the element containing the physical point `x`.
    ///
    /// This is a linear scan; it avoids storing a spatial index for the
    /// low-frequency point queries this API targets.
    fn locate(&self, x: [f64; D]) -> Option<(&Self::Item, [f64; P])> {
        self.iter()
            .find_map(|elm| elm.local_coordinates(x).map(|xi| (elm, xi)))
    }

    /// Evaluates the field at the physical point `x`.
    ///
    /// `coefficients` uses the global layout `node * C + component`. Returns
    /// `Ok(None)` if `x` is outside every element.
    ///
    /// # Errors
    ///
    /// Returns [`VariableError::CoefficientOutOfRange`] or
    /// [`VariableError::DofOverflow`] for an invalid coefficient layout.
    fn value(
        &self,
        x: [f64; D],
        coefficients: &[f64],
    ) -> Result<Option<[f64; C]>, VariableError> {
        let Some((elm, xi)) = self.locate(x) else {
            return Ok(None);
        };

        let shape = elm.shape(xi);
        let nodes_id = elm.id();
        let mut acc = [0.0; C];

        for a in 0..N {
            let coefficients_at_node = node_coefficients(coefficients, nodes_id[a], C)?;
            for k in 0..C {
                acc[k] += shape[a] * coefficients_at_node[k];
            }
        }

        Ok(Some(acc))
    }

    /// Derivative of each independent field component at physical point `x`.
    ///
    /// The result is indexed as `derivative[i][k] = ∂u_k / ∂x_i`, with shape
    /// `[[f64; C]; D]`. Returns `Ok(None)` if `x` is outside every element.
    fn derivative_value(
        &self,
        x: [f64; D],
        coefficients: &[f64],
    ) -> Result<Option<[[f64; C]; D]>, VariableError> {
        let Some((elm, xi)) = self.locate(x) else {
            return Ok(None);
        };

        let dshape = elm.derivative_shape(xi);
        let node_ids = elm.id();
        let mut derivative = [[0.0_f64; C]; D];

        for a in 0..N {
            let coefficients_at_node = node_coefficients(coefficients, node_ids[a], C)?;
            for (i, row) in derivative.iter_mut().enumerate() {
                for (k, value) in row.iter_mut().enumerate() {
                    *value += dshape[a][i] * coefficients_at_node[k];
                }
            }
        }

        Ok(Some(derivative))
    }
}
impl<T, const D: usize, const C: usize, const P: usize, const N: usize>
    Variable<D, C, P, N> for &[T]
where
    T: Element<D, P, N>,
{
    type Item = T;

    fn iter(&self) -> impl Iterator<Item = &T> {
        <[T]>::iter(self)
    }
}
impl<T, const D: usize, const C: usize, const P: usize, const N: usize>
    Variable<D, C, P, N> for Vec<T>
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
        assert_eq!(Variable::<1, 1, 1, 2>::num_dofs(&bars), 4);
        assert_eq!(Variable::<1, 1, 1, 2>::num_dofs(&bars.as_slice()), 4);
    }

    #[test]
    fn num_dofs_empty_variable_is_zero() {
        let bars: Vec<Seg2<1>> = Vec::new();
        assert_eq!(Variable::<1, 1, 1, 2>::num_dofs(&bars), 0);
    }

    #[test]
    fn num_dofs_scales_with_tensor_order() {
        // One Seg2<2> element: two unique nodes.
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
        )];
        assert_eq!(Variable::<2, 1, 1, 2>::num_dofs(&bars), 2); // scalar: 2 * 1
        assert_eq!(Variable::<2, 2, 1, 2>::num_dofs(&bars), 4); // vector: 2 * 2
        assert_eq!(Variable::<2, 3, 1, 2>::num_dofs(&bars), 6); // symmetric 2nd-order: 2 * 3
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
            Variable::<1, 1, 1, 2>::value(&bars, [1.0], &coefficients).unwrap(),
            Some([1.0])
        );
        assert_eq!(
            Variable::<1, 1, 1, 2>::value(&bars, [3.0], &coefficients).unwrap(),
            Some([3.0])
        );
        assert_eq!(
            Variable::<1, 1, 1, 2>::value(&bars, [5.0], &coefficients).unwrap(),
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
            Variable::<2, 2, 1, 2>::value(&bars, [1.0, 0.0], &coefficients).unwrap(),
            Some([2.0, 1.0])
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
            Variable::<2, 3, 1, 2>::value(&bars, [1.0, 0.0], &coefficients).unwrap(),
            Some([2.5, 3.5, 4.5])
        );
    }

    #[test]
    fn value_rejects_short_coefficients() {
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [1.0]))];

        let err = Variable::<1, 1, 1, 2>::value(&bars, [0.5], &[0.0]).unwrap_err();
        assert!(matches!(
            err,
            super::VariableError::CoefficientOutOfRange { .. }
        ));
    }

    #[test]
    fn derivative_value_of_linear_vector_field() {
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
        )];
        // u = (1 + 2x, 2 - 2x)
        let coefficients = [1.0, 2.0, 3.0, 0.0];

        let derivative = Variable::<2, 2, 1, 2>::derivative_value(&bars, [0.5, 0.0], &coefficients)
            .unwrap()
            .unwrap();
        assert_eq!(derivative, [[2.0, -2.0], [0.0, 0.0]]);
    }

}
