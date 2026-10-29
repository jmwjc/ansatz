//! Piecewise-constant (P0) triangle: [`P0Tri3`].
//!
//! `P0Tri3` is the scalar element whose local space is the constants on a
//! triangle: one degree of freedom per cell, shape function `1`, zero gradient.
//! It is the natural partner of the lowest-order Raviart–Thomas flux
//! ([`crate::element::rttri3::RTTri3`]) in the mixed (Hellinger–Reissner) heat
//! discretization, where the temperature is an `L²` multiplier of the balance
//! equation; see `book/src/theory-heat.md` and `book/src/theory-rt.md`.
//!
//! The element is named after its polynomial degree rather than its dof count
//! (the crate-wide `<shape><dofs>` convention would give `Tri1`): the single
//! dof belongs to the **cell**, not to a vertex, so `nodes`/`vertices` have no
//! geometric meaning for the approximation. The geometry (Jacobian, inverse
//! map, containment) is delegated to [`crate::element::tri3::Tri3`] through the
//! three corner vertices, while the single node stored at the centroid only
//! carries the global cell-dof id.

use crate::element::tri3::Tri3;
use crate::element::{
    Dunavant, Element, ElementBase, FullIntegration, IntegrationScheme, ReducedIntegration,
};
use crate::node::Node;

/// Piecewise-constant (P0) triangle with one cell degree of freedom.
///
/// The local space is `span{1}`, so `shape = [1]` and `derivative_shape = [0]`;
/// `jacobe` is `2 × area` and `local_coordinates` is the inverse map of the
/// underlying triangle. The single node is the centroid of the triangle and its
/// [`Node::id`] is the global cell dof.
#[derive(Clone, Debug, PartialEq)]
pub struct P0Tri3<const D: usize> {
    vertices: [Node<D>; 3],
    centroid: Node<D>,
}

impl<const D: usize> P0Tri3<D> {
    /// Builds the element from its three corner vertices and the global id of
    /// the cell dof.
    pub fn new(vertices: [Node<D>; 3], dof_id: usize) -> Self {
        let coordinates: [f64; D] = std::array::from_fn(|d| {
            vertices
                .iter()
                .map(|vertex| vertex.coordinates[d])
                .sum::<f64>()
                / 3.0
        });
        P0Tri3 {
            vertices,
            centroid: Node::new(dof_id, coordinates),
        }
    }

    /// Builds the P0 element of the triangle `tri`, with `dof_id` as the global
    /// id of the cell dof.
    pub fn from_tri3(tri: &Tri3<D>, dof_id: usize) -> Self {
        let id = tri.id();
        let coordinates = tri.vertices_coordinates();
        Self::new(
            [
                Node::new(id[0], coordinates[0]),
                Node::new(id[1], coordinates[1]),
                Node::new(id[2], coordinates[2]),
            ],
            dof_id,
        )
    }

    /// The three corner vertices, in the local ordering of the underlying
    /// [`Tri3`].
    pub fn vertices(&self) -> &[Node<D>; 3] {
        &self.vertices
    }

    /// The centroid (the single node), which carries the cell dof id.
    pub fn centroid(&self) -> &Node<D> {
        &self.centroid
    }

    /// The underlying linear triangle, used for the geometry (Jacobian,
    /// inverse map, containment).
    fn as_tri3(&self) -> Tri3<D> {
        Tri3::new(self.vertices[0], self.vertices[1], self.vertices[2])
    }
}

impl<const D: usize> IntegrationScheme<2, 1> for P0Tri3<D> {
    fn points() -> [[f64; 2]; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<2, 3> for P0Tri3<D> {
    fn points() -> [[f64; 2]; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::points()
    }
    fn weights() -> [f64; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::weights()
    }
}

impl<const D: usize> FullIntegration<2, 3> for P0Tri3<D> {}
impl<const D: usize> ReducedIntegration<2, 1> for P0Tri3<D> {}

impl<const D: usize> ElementBase<D, 2, 1> for P0Tri3<D> {
    fn id(&self) -> [usize; 1] {
        [self.centroid.id]
    }

    /// The single node is the centroid; it carries the cell dof id but is not a
    /// vertex of the approximation.
    fn vertices_coordinates(&self) -> [[f64; D]; 1] {
        [self.centroid.coordinates]
    }

    fn coordinates(&self, parametric_coordinates: [f64; 2]) -> [f64; D] {
        self.as_tri3().coordinates(parametric_coordinates)
    }

    fn jacobe(&self, parametric_coordinates: [f64; 2]) -> f64 {
        self.as_tri3().jacobe(parametric_coordinates)
    }

    fn jacobe_mat(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 2] {
        self.as_tri3().jacobe_mat(parametric_coordinates)
    }

    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 2]> {
        self.as_tri3().local_coordinates(x)
    }
}

impl<const D: usize> Element<D, 2, 1> for P0Tri3<D> {
    /// The constant shape function `1`.
    fn shape(&self, _parametric_coordinates: [f64; 2]) -> [f64; 1] {
        [1.0]
    }

    /// The gradient of a constant is zero.
    fn derivative_shape(&self, _parametric_coordinates: [f64; 2]) -> [[f64; D]; 1] {
        [[0.0; D]; 1]
    }
}

#[cfg(test)]
mod tests {
    use super::P0Tri3;
    use crate::element::{Element, ElementBase, IntegrationScheme};
    use crate::node::Node;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    /// Right triangle `(0,0)`, `(3,0)`, `(0,4)`: area 6, centroid `(1, 4/3)`.
    fn right_triangle() -> P0Tri3<2> {
        P0Tri3::new(
            [
                Node::new(0, [0.0, 0.0]),
                Node::new(1, [3.0, 0.0]),
                Node::new(2, [0.0, 4.0]),
            ],
            7,
        )
    }

    #[test]
    fn p0tri3_shape_gradient_and_geometry() {
        let p0 = right_triangle();
        assert_eq!(p0.id(), [7]);
        assert_eq!(p0.shape([0.25, 0.5]), [1.0]);
        assert_eq!(p0.derivative_shape([0.25, 0.5]), [[0.0, 0.0]]);

        assert!(close(p0.jacobe([0.0, 0.0]), 12.0));
        assert!(close(p0.coordinates([1.0 / 3.0, 1.0 / 3.0])[0], 1.0));
        assert!(close(p0.coordinates([1.0 / 3.0, 1.0 / 3.0])[1], 4.0 / 3.0));
        assert_eq!(p0.vertices_coordinates(), [[1.0, 4.0 / 3.0]]);
        assert!(p0.contains([1.0, 1.0]));
        assert!(!p0.contains([3.0, 4.0]));
    }

    #[test]
    fn p0tri3_integrates_constants_exactly() {
        // ∫_T 1 dΩ = area = 6, with both integration rules (the integrand is constant).
        let p0 = right_triangle();
        let integrate = |points: &[[f64; 2]], weights: &[f64]| -> f64 {
            points
                .iter()
                .zip(weights)
                .map(|(xi, w)| p0.shape(*xi)[0] * p0.jacobe(*xi) * w)
                .sum()
        };
        let one = <P0Tri3<2> as IntegrationScheme<2, 1>>::points();
        let one_w = <P0Tri3<2> as IntegrationScheme<2, 1>>::weights();
        let full = <P0Tri3<2> as IntegrationScheme<2, 3>>::points();
        let full_w = <P0Tri3<2> as IntegrationScheme<2, 3>>::weights();
        assert!(close(integrate(&one, &one_w), 6.0));
        assert!(close(integrate(&full, &full_w), 6.0));
    }

    #[test]
    fn p0tri3_from_tri3_uses_the_cell_dof() {
        use crate::element::tri3::Tri3;

        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        let p0 = P0Tri3::from_tri3(&tri, 11);
        assert_eq!(p0.id(), [11]);
        assert_eq!(p0.vertices(), tri.nodes());
        assert!(close(p0.jacobe([0.0, 0.0]), 12.0));
    }
}
