use crate::element::{
    Dunavant, Element, ElementBase, FullIntegration, IntegrationScheme, ReducedIntegration,
};
use crate::node::Node;

#[derive(Clone, Debug, PartialEq)]
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
impl<const D: usize> Tri3<D> {
    /// The element's nodes, in local ordering.
    pub fn nodes(&self) -> &[Node<D>; 3] {
        &self.nodes
    }

    /// Node at local index, or `None` if out of bounds.
    pub fn node(&self, index: usize) -> Option<&Node<D>> {
        self.nodes.get(index)
    }
}

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

impl<const D: usize> ElementBase<D, 2, 3> for Tri3<D> {
    fn id(&self) -> [usize; 3] {
        [self.nodes[0].id, self.nodes[1].id, self.nodes[2].id]
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
        let aa: f64 = j[0].iter().map(|v| v * v).sum();
        let bb: f64 = j[1].iter().map(|v| v * v).sum();
        let ab: f64 = j[0].iter().zip(j[1]).map(|(a, b)| a * b).sum();
        (aa * bb - ab * ab).sqrt()
    }

    /// `J[0][i] = ∂x_i/∂ξ`, `J[1][i] = ∂x_i/∂η` (constant for a linear triangle).
    fn jacobe_mat(&self, _parametric_coordinates: [f64; 2]) -> [[f64; D]; 2] {
        [
            std::array::from_fn(|i| self.nodes[1].coordinates[i] - self.nodes[0].coordinates[i]),
            std::array::from_fn(|i| self.nodes[2].coordinates[i] - self.nodes[0].coordinates[i]),
        ]
    }

    /// Solves the least-squares projection of `x` onto the triangle plane for
    /// the barycentric coordinates `(λ₂, λ₃)` (with `λ₁ = 1 − λ₂ − λ₃`, the
    /// reference triangle being `(0,0)`, `(1,0)`, `(0,1)`). Returns `Some` when
    /// all three barycentric coordinates are non-negative and the out-of-plane
    /// residual (for `D = 3`) is within tolerance.
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 2]> {
        let tol = super::CONTAINS_TOL;
        let a = self.nodes[0].coordinates;
        let b = self.nodes[1].coordinates;
        let c = self.nodes[2].coordinates;
        let u: [f64; D] = std::array::from_fn(|i| b[i] - a[i]);
        let v: [f64; D] = std::array::from_fn(|i| c[i] - a[i]);
        let w: [f64; D] = std::array::from_fn(|i| x[i] - a[i]);

        let uu: f64 = (0..D).map(|i| u[i] * u[i]).sum();
        let uv: f64 = (0..D).map(|i| u[i] * v[i]).sum();
        let vv: f64 = (0..D).map(|i| v[i] * v[i]).sum();
        let wu: f64 = (0..D).map(|i| w[i] * u[i]).sum();
        let wv: f64 = (0..D).map(|i| w[i] * v[i]).sum();
        let det = uu * vv - uv * uv;
        if det.abs() <= 1e-16 * (uu * vv) {
            return None;
        }

        let l2 = (wu * vv - wv * uv) / det;
        let l3 = (uu * wv - uv * wu) / det;
        let l1 = 1.0 - l2 - l3;
        if l1 < -tol || l2 < -tol || l3 < -tol {
            return None;
        }

        let r2: f64 = (0..D).map(|i| (w[i] - l2 * u[i] - l3 * v[i]).powi(2)).sum();
        let scale = uu.max(vv).max(1.0);
        if r2 > tol * tol * scale {
            return None;
        }
        Some([l2, l3])
    }
}

impl<const D: usize> Element<D, 2, 3> for Tri3<D> {
    /// Area-coordinate shape functions: nodes 1, 2, 3 correspond to the reference
    /// triangle vertices `(0, 0)`, `(1, 0)`, `(0, 1)`.
    fn shape(&self, parametric_coordinates: [f64; 2]) -> [f64; 3] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        [1.0 - ξ - η, ξ, η]
    }

    /// Physical gradient of the area-coordinate shape functions, computed
    /// directly via the metric tensor: `∇ₓN = j_α g^{αβ} ∂_β N`. Valid for
    /// both `D = 2` (plane triangle) and `D = 3` (spatial triangle).
    fn derivative_shape(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 3] {
        let j = self.jacobe_mat(parametric_coordinates);
        let j0 = &j[0]; // ∂x/∂ξ
        let j1 = &j[1]; // ∂x/∂η

        // Metric tensor g_{αβ} = j_α · j_β.
        let g00: f64 = j0.iter().map(|v| v * v).sum();
        let g01: f64 = j0.iter().zip(j1).map(|(a, b)| a * b).sum();
        let g11: f64 = j1.iter().map(|v| v * v).sum();
        let det = g00 * g11 - g01 * g01;

        // Inverse metric g^{αβ}.
        let gi00 = g11 / det;
        let gi01 = -g01 / det;
        let gi11 = g00 / det;

        // Parametric gradients ∂_β N_a = (∂N_a/∂ξ, ∂N_a/∂η).
        let dn: [[f64; 2]; 3] = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];

        // ∇ₓN_a = j_ξ (g^{ξξ}∂_ξN + g^{ξη}∂_ηN) + j_η (g^{ηξ}∂_ξN + g^{ηη}∂_ηN).
        std::array::from_fn(|a| {
            let (dξ, dη) = (dn[a][0], dn[a][1]);
            let c0 = gi00 * dξ + gi01 * dη;
            let c1 = gi01 * dξ + gi11 * dη;
            std::array::from_fn(|i| j0[i] * c0 + j1[i] * c1)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Tri3;
    use crate::element::{Element, ElementBase};
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
    fn tri3_contains_and_local_coordinates() {
        // Right triangle (0,0), (3,0), (0,4): centroid barycentric = (1/3, 1/3).
        let t = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        assert!(t.contains([1.0, 4.0 / 3.0]));
        let lc = t.local_coordinates([1.0, 4.0 / 3.0]).unwrap();
        assert!(close(lc[0], 1.0 / 3.0));
        assert!(close(lc[1], 1.0 / 3.0));
        assert!(!t.contains([3.0, 4.0]));

        // Spatial triangle in the xy-plane: off-plane points are rejected.
        let t3 = Tri3::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [0.0, 3.0, 0.0]),
        );
        assert!(t3.contains([1.0, 1.0, 0.0]));
        assert!(!t3.contains([1.0, 1.0, 0.1]));
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

    #[test]
    fn tri3_physical_gradient() {
        // Right triangle (0,0), (3,0), (0,4): physical gradients
        // ∇N1 = (-1/3, -1/4), ∇N2 = (1/3, 0), ∇N3 = (0, 1/4).
        let t = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        let g = t.derivative_shape([0.0, 0.0]);
        assert!(close(g[0][0], -1.0 / 3.0));
        assert!(close(g[0][1], -1.0 / 4.0));
        assert!(close(g[1][0], 1.0 / 3.0));
        assert!(close(g[1][1], 0.0));
        assert!(close(g[2][0], 0.0));
        assert!(close(g[2][1], 1.0 / 4.0));
    }

    #[test]
    fn tri3_physical_gradient_3d() {
        // Spatial triangle (0,0,0), (2,0,0), (0,3,0): physical gradients
        // ∇N1 = (-1/2, -1/3, 0), ∇N2 = (1/2, 0, 0), ∇N3 = (0, 1/3, 0).
        let t = Tri3::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [0.0, 3.0, 0.0]),
        );
        let g = t.derivative_shape([0.0, 0.0]);
        assert!(close(g[0][0], -0.5));
        assert!(close(g[0][1], -1.0 / 3.0));
        assert!(close(g[0][2], 0.0));
        assert!(close(g[1][0], 0.5));
        assert!(close(g[1][1], 0.0));
        assert!(close(g[1][2], 0.0));
        assert!(close(g[2][0], 0.0));
        assert!(close(g[2][1], 1.0 / 3.0));
        assert!(close(g[2][2], 0.0));
    }
}
