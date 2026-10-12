use crate::element::{
    Element, FullIntegration, IntegrationScheme, ReducedIntegration, TensorProductQuad,
};
use crate::node::Node;

pub struct Quad4<const D: usize> {
    nodes: [Node<D>; 4],
}

impl<const D: usize> Quad4<D> {
    pub fn new(n1: Node<D>, n2: Node<D>, n3: Node<D>, n4: Node<D>) -> Self {
        Quad4 {
            nodes: [n1, n2, n3, n4],
        }
    }
}

// Quadrilateral full integration: 2×2 tensor product (4 points, degree 3); reduced integration: 1 point (centroid).
impl<const D: usize> IntegrationScheme<2, 1> for Quad4<D> {
    fn points() -> [[f64; 2]; 1] {
        <TensorProductQuad<1> as IntegrationScheme<2, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <TensorProductQuad<1> as IntegrationScheme<2, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<2, 4> for Quad4<D> {
    fn points() -> [[f64; 2]; 4] {
        <TensorProductQuad<2> as IntegrationScheme<2, 4>>::points()
    }
    fn weights() -> [f64; 4] {
        <TensorProductQuad<2> as IntegrationScheme<2, 4>>::weights()
    }
}

impl<const D: usize> FullIntegration<2, 4> for Quad4<D> {}
impl<const D: usize> ReducedIntegration<2, 1> for Quad4<D> {}

impl<const D: usize> Element<D, 2, 4> for Quad4<D> {
    fn id(&self) -> [usize; 4] {
        [
            self.nodes[0].id,
            self.nodes[1].id,
            self.nodes[2].id,
            self.nodes[3].id,
        ]
    }

    /// Bilinear shape functions: nodes 1, 2, 3, 4 correspond respectively to the reference
    /// quadrilateral vertices `(-1, -1)`, `(1, -1)`, `(1, 1)`, `(-1, 1)`.
    fn shape(&self, parametric_coordinates: [f64; 2]) -> [f64; 4] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        [
            0.25 * (1.0 - ξ) * (1.0 - η),
            0.25 * (1.0 + ξ) * (1.0 - η),
            0.25 * (1.0 + ξ) * (1.0 + η),
            0.25 * (1.0 - ξ) * (1.0 + η),
        ]
    }

    /// Partial derivative with respect to the first parametric coordinate: `∂N/∂ξ`.
    ///
    /// Note: `Element`'s method returns `[f64; N]`, which cannot hold the full gradient
    /// of a 2D element; the derivative with respect to `η` and the physical gradient
    /// must be obtained separately, e.g. via [`Self::jacobe_mat`].
    fn derivative_shape(&self, parametric_coordinates: [f64; 2]) -> [f64; 4] {
        let η = parametric_coordinates[1];
        [
            -0.25 * (1.0 - η),
            0.25 * (1.0 - η),
            0.25 * (1.0 + η),
            -0.25 * (1.0 + η),
        ]
    }

    fn vertices_coordinates(&self) -> [[f64; D]; 4] {
        [
            self.nodes[0].coordinates,
            self.nodes[1].coordinates,
            self.nodes[2].coordinates,
            self.nodes[3].coordinates,
        ]
    }

    fn coordinates(&self, parametric_coordinates: [f64; 2]) -> [f64; D] {
        let shape = self.shape(parametric_coordinates);
        let vertices = self.vertices_coordinates();
        std::array::from_fn(|i| (0..4).map(|k| vertices[k][i] * shape[k]).sum())
    }

    /// Area ratio between the reference and physical quadrilaterals: the square root of
    /// the determinant of `JᵀJ` (i.e. `|det J|` when `D = 2`; for `D = 3` it is the norm
    /// of the cross product of the two tangent vectors).
    fn jacobe(&self, parametric_coordinates: [f64; 2]) -> f64 {
        let j = self.jacobe_mat(parametric_coordinates);
        let aa: f64 = (0..D).map(|i| j[i][0] * j[i][0]).sum();
        let bb: f64 = (0..D).map(|i| j[i][1] * j[i][1]).sum();
        let ab: f64 = (0..D).map(|i| j[i][0] * j[i][1]).sum();
        (aa * bb - ab * ab).sqrt()
    }

    /// `J[i][0] = ∂x_i/∂ξ`, `J[i][1] = ∂x_i/∂η` (bilinear mapping, varies with `(ξ, η)`).
    fn jacobe_mat(&self, parametric_coordinates: [f64; 2]) -> [[f64; 2]; D] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        let dn_dξ = [
            -0.25 * (1.0 - η),
            0.25 * (1.0 - η),
            0.25 * (1.0 + η),
            -0.25 * (1.0 + η),
        ];
        let dn_dη = [
            -0.25 * (1.0 - ξ),
            -0.25 * (1.0 + ξ),
            0.25 * (1.0 + ξ),
            0.25 * (1.0 - ξ),
        ];
        std::array::from_fn(|i| {
            [
                (0..4).map(|k| self.nodes[k].coordinates[i] * dn_dξ[k]).sum(),
                (0..4).map(|k| self.nodes[k].coordinates[i] * dn_dη[k]).sum(),
            ]
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Quad4;
    use crate::element::Element;
    use crate::node::Node;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn quad4_shape_and_jacobe() {
        // Unit square (0,0), (1,0), (1,1), (0,1): area 1, jacobe(0,0) = 1/4.
        let q = Quad4::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
            Node::new(2, [1.0, 1.0]),
            Node::new(3, [0.0, 1.0]),
        );
        assert_eq!(q.id(), [0, 1, 2, 3]);
        // Shape function values at the vertices
        assert_eq!(q.shape([-1.0, -1.0]), [1.0, 0.0, 0.0, 0.0]);
        assert_eq!(q.shape([1.0, -1.0]), [0.0, 1.0, 0.0, 0.0]);
        assert_eq!(q.shape([1.0, 1.0]), [0.0, 0.0, 1.0, 0.0]);
        assert_eq!(q.shape([-1.0, 1.0]), [0.0, 0.0, 0.0, 1.0]);
        // The centroid (ξ, η) = (0, 0) maps to the physical centroid (0.5, 0.5)
        let c = q.coordinates([0.0, 0.0]);
        assert!(close(c[0], 0.5));
        assert!(close(c[1], 0.5));
        // jacobe = |det J| = 1/4
        assert!(close(q.jacobe([0.0, 0.0]), 0.25));
    }

    #[test]
    fn quad4_jacobe_3d() {
        // Spatial quadrilateral (0,0,0), (2,0,0), (2,3,0), (0,3,0): area 6, jacobe = 6/4 = 1.5.
        let q = Quad4::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [2.0, 3.0, 0.0]),
            Node::new(3, [0.0, 3.0, 0.0]),
        );
        assert!(close(q.jacobe([0.0, 0.0]), 1.5));
    }
}
