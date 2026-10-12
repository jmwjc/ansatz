use crate::element::{Dunavant, Element, FullIntegration, IntegrationScheme, ReducedIntegration};
use crate::node::Node;

pub struct Tri3<const D: usize> {
    nodes: [Node<D>; 3],
}

impl<const D: usize> Tri3<D> {
    pub fn new(n1: Node<D>, n2: Node<D>, n3: Node<D>) -> Self {
        Tri3 {
            nodes: [n1, n2, n3],
        }
    }
}

// Full integration for the triangle: 3-point Dunavant rule (algebraic degree 2);
// reduced integration: 1 point (centroid).
impl<const D: usize> IntegrationScheme<2, 1> for Tri3<D> {
    fn points() -> [[f64; 2]; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<2, 3> for Tri3<D> {
    fn points() -> [[f64; 2]; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::points()
    }
    fn weights() -> [f64; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::weights()
    }
}

impl<const D: usize> FullIntegration<2, 3> for Tri3<D> {}
impl<const D: usize> ReducedIntegration<2, 1> for Tri3<D> {}

impl<const D: usize> Element<D, 2, 3> for Tri3<D> {
    fn id(&self) -> [usize; 3] {
        [self.nodes[0].id, self.nodes[1].id, self.nodes[2].id]
    }

    /// Area-coordinate shape functions: nodes 1, 2, 3 correspond to the reference
    /// triangle vertices `(0, 0)`, `(1, 0)`, `(0, 1)`.
    fn shape(&self, parametric_coordinates: [f64; 2]) -> [f64; 3] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        [1.0 - ξ - η, ξ, η]
    }

    /// Partial derivative with respect to the first parametric coordinate `ξ`, `∂N/∂ξ`.
    ///
    /// Note: the `Element` trait method returns `[f64; N]`, which cannot hold the full
    /// gradient of a 2D element; the derivative with respect to `η` and the physical
    /// gradient must be obtained separately, combined with [`Self::jacobe_mat`].
    fn derivative_shape(&self, _parametric_coordinates: [f64; 2]) -> [f64; 3] {
        [-1.0, 1.0, 0.0]
    }

    fn vertices_coordinates(&self) -> [[f64; D]; 3] {
        [
            self.nodes[0].coordinates,
            self.nodes[1].coordinates,
            self.nodes[2].coordinates,
        ]
    }

    fn coordinates(&self, parametric_coordinates: [f64; 2]) -> [f64; D] {
        let shape = self.shape(parametric_coordinates);
        let vertices = self.vertices_coordinates();
        std::array::from_fn(|i| {
            vertices[0][i] * shape[0] + vertices[1][i] * shape[1] + vertices[2][i] * shape[2]
        })
    }

    /// Area ratio between the reference triangle and the physical triangle: the square root
    /// of the determinant of `JᵀJ` (for `D = 2` this is `|det J|`; for `D = 3` it is the
    /// norm of the cross product of the two tangent vectors).
    fn jacobe(&self, parametric_coordinates: [f64; 2]) -> f64 {
        let j = self.jacobe_mat(parametric_coordinates);
        let aa: f64 = (0..D).map(|i| j[i][0] * j[i][0]).sum();
        let bb: f64 = (0..D).map(|i| j[i][1] * j[i][1]).sum();
        let ab: f64 = (0..D).map(|i| j[i][0] * j[i][1]).sum();
        (aa * bb - ab * ab).sqrt()
    }

    /// `J[i][0] = ∂x_i/∂ξ`, `J[i][1] = ∂x_i/∂η` (constant for a linear triangle).
    fn jacobe_mat(&self, _parametric_coordinates: [f64; 2]) -> [[f64; 2]; D] {
        std::array::from_fn(|i| {
            [
                self.nodes[1].coordinates[i] - self.nodes[0].coordinates[i],
                self.nodes[2].coordinates[i] - self.nodes[0].coordinates[i],
            ]
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Tri3;
    use crate::element::Element;
    use crate::node::Node;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn tri3_shape_and_area() {
        // Right triangle (0,0), (3,0), (0,4): area 6, jacobe should be 2×area = 12.
        let t = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        assert_eq!(t.id(), [0, 1, 2]);
        // Shape function values at the vertices
        assert_eq!(t.shape([0.0, 0.0]), [1.0, 0.0, 0.0]);
        assert_eq!(t.shape([1.0, 0.0]), [0.0, 1.0, 0.0]);
        assert_eq!(t.shape([0.0, 1.0]), [0.0, 0.0, 1.0]);
        // Centroid (ξ, η) = (1/3, 1/3) maps to the physical centroid (1, 4/3)
        let c = t.coordinates([1.0 / 3.0, 1.0 / 3.0]);
        assert!(close(c[0], 1.0));
        assert!(close(c[1], 4.0 / 3.0));
        // jacobe = |det J| = 2×area = 12
        assert!(close(t.jacobe([0.0, 0.0]), 12.0));
    }

    #[test]
    fn tri3_jacobe_3d() {
        // Spatial triangle (0,0,0), (2,0,0), (0,3,0): area 3, jacobe = 2×area = 6.
        let t = Tri3::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [0.0, 3.0, 0.0]),
        );
        assert!(close(t.jacobe([0.0, 0.0]), 6.0));
    }
}
