use crate::element::{
    Element, ElementBase, FullIntegration, IntegrationScheme, ReducedIntegration, TensorProductQuad,
};
use crate::node::Node;

#[derive(Clone, Debug, PartialEq)]
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
impl<const D: usize> Quad4<D> {
    /// The element's nodes, in local ordering.
    pub fn nodes(&self) -> &[Node<D>; 4] {
        &self.nodes
    }

    /// Node at local index, or `None` if out of bounds.
    pub fn node(&self, index: usize) -> Option<&Node<D>> {
        self.nodes.get(index)
    }
}

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

impl<const D: usize> ElementBase<D, 2, 4> for Quad4<D> {
    fn id(&self) -> [usize; 4] {
        [
            self.nodes[0].id,
            self.nodes[1].id,
            self.nodes[2].id,
            self.nodes[3].id,
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
        let aa: f64 = j[0].iter().map(|v| v * v).sum();
        let bb: f64 = j[1].iter().map(|v| v * v).sum();
        let ab: f64 = j[0].iter().zip(j[1]).map(|(a, b)| a * b).sum();
        (aa * bb - ab * ab).sqrt()
    }

    /// `J[0][i] = ∂x_i/∂ξ`, `J[1][i] = ∂x_i/∂η` (bilinear mapping, varies with `(ξ, η)`).
    fn jacobe_mat(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 2] {
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
        [
            std::array::from_fn(|i| {
                (0..4)
                    .map(|k| self.nodes[k].coordinates[i] * dn_dξ[k])
                    .sum()
            }),
            std::array::from_fn(|i| {
                (0..4)
                    .map(|k| self.nodes[k].coordinates[i] * dn_dη[k])
                    .sum()
            }),
        ]
    }

    /// Inverse of the bilinear isoparametric map, solved with a Gauss–Newton
    /// (least-squares) iteration starting from the reference centroid `(0, 0)`.
    /// Returns `Some([ξ, η])` only when the iteration converges inside the
    /// reference square `[-1, 1]²` and the physical residual is within
    /// tolerance (relevant for `D = 3`).
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 2]> {
        let tol = super::CONTAINS_TOL;
        let mut ξ = 0.0_f64;
        let mut η = 0.0_f64;
        for _ in 0..20 {
            let j = self.jacobe_mat([ξ, η]);
            let g00: f64 = (0..D).map(|i| j[0][i] * j[0][i]).sum();
            let g01: f64 = (0..D).map(|i| j[0][i] * j[1][i]).sum();
            let g11: f64 = (0..D).map(|i| j[1][i] * j[1][i]).sum();
            let det = g00 * g11 - g01 * g01;
            if det.abs() <= 1e-20 {
                return None;
            }
            let f = self.coordinates([ξ, η]);
            let r: [f64; D] = std::array::from_fn(|i| x[i] - f[i]);
            let r0: f64 = (0..D).map(|i| j[0][i] * r[i]).sum();
            let r1: f64 = (0..D).map(|i| j[1][i] * r[i]).sum();
            let dξ = (g11 * r0 - g01 * r1) / det;
            let dη = (g00 * r1 - g01 * r0) / det;
            ξ += dξ;
            η += dη;
            if dξ.abs() < 1e-13 && dη.abs() < 1e-13 {
                break;
            }
        }

        if ξ < -1.0 - tol || ξ > 1.0 + tol || η < -1.0 - tol || η > 1.0 + tol {
            return None;
        }

        let f = self.coordinates([ξ, η]);
        let r2: f64 = (0..D).map(|i| (x[i] - f[i]).powi(2)).sum();
        let vertices = self.vertices_coordinates();
        let mut scale = 0.0_f64;
        for k in 0..4 {
            let d2: f64 = (0..D)
                .map(|i| (vertices[k][i] - vertices[0][i]).powi(2))
                .sum();
            scale = scale.max(d2);
        }
        if r2 > tol * tol * scale.max(1.0) {
            return None;
        }
        Some([ξ, η])
    }
}

impl<const D: usize> Element<D, 2, 4> for Quad4<D> {
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

    /// Physical gradient of the bilinear shape functions, computed directly via
    /// the metric tensor: `∇ₓN = j_α g^{αβ} ∂_β N`. Valid for both `D = 2`
    /// (plane quadrilateral) and `D = 3` (spatial quadrilateral).
    fn derivative_shape(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 4] {
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
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        let dn: [[f64; 2]; 4] = [
            [-0.25 * (1.0 - η), -0.25 * (1.0 - ξ)],
            [0.25 * (1.0 - η), -0.25 * (1.0 + ξ)],
            [0.25 * (1.0 + η), 0.25 * (1.0 + ξ)],
            [-0.25 * (1.0 + η), 0.25 * (1.0 - ξ)],
        ];

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
    use super::Quad4;
    use crate::element::{Element, ElementBase};
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

    #[test]
    fn quad4_physical_gradient() {
        // Unit square at its centroid: ∇N = (±1/2, ±1/2).
        let q = Quad4::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
            Node::new(2, [1.0, 1.0]),
            Node::new(3, [0.0, 1.0]),
        );
        let g = q.derivative_shape([0.0, 0.0]);
        assert!(close(g[0][0], -0.5));
        assert!(close(g[0][1], -0.5));
        assert!(close(g[1][0], 0.5));
        assert!(close(g[1][1], -0.5));
        assert!(close(g[2][0], 0.5));
        assert!(close(g[2][1], 0.5));
        assert!(close(g[3][0], -0.5));
        assert!(close(g[3][1], 0.5));
    }

    #[test]
    fn quad4_contains_and_local_coordinates() {
        // Unit square mapped to [0,2]×[0,2]: reference centroid (0,0) maps to (1,1).
        let q = Quad4::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [2.0, 0.0]),
            Node::new(2, [2.0, 2.0]),
            Node::new(3, [0.0, 2.0]),
        );
        assert!(q.contains([1.0, 1.0]));
        let lc = q.local_coordinates([1.0, 1.0]).unwrap();
        assert!(close(lc[0], 0.0));
        assert!(close(lc[1], 0.0));
        assert!(q.contains([2.0, 2.0]));
        assert!(!q.contains([3.0, 1.0]));

        // Distorted (trapezoidal) quadrilateral: reference centroid maps to the
        // average of the vertices.
        let q2 = Quad4::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [2.0, 0.0]),
            Node::new(2, [3.0, 2.0]),
            Node::new(3, [0.0, 2.0]),
        );
        assert!(q2.contains([1.25, 1.0]));
    }

    #[test]
    fn quad4_physical_gradient_3d() {
        // Spatial quadrilateral (0,0,0), (2,0,0), (2,3,0), (0,3,0) at its
        // centroid: ∇N1 = (-1/4, -1/6, 0), ∇N2 = (1/4, -1/6, 0),
        // ∇N3 = (1/4, 1/6, 0), ∇N4 = (-1/4, 1/6, 0).
        let q = Quad4::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [2.0, 3.0, 0.0]),
            Node::new(3, [0.0, 3.0, 0.0]),
        );
        let g = q.derivative_shape([0.0, 0.0]);
        assert!(close(g[0][0], -0.25));
        assert!(close(g[0][1], -1.0 / 6.0));
        assert!(close(g[0][2], 0.0));
        assert!(close(g[1][0], 0.25));
        assert!(close(g[1][1], -1.0 / 6.0));
        assert!(close(g[1][2], 0.0));
        assert!(close(g[2][0], 0.25));
        assert!(close(g[2][1], 1.0 / 6.0));
        assert!(close(g[2][2], 0.0));
        assert!(close(g[3][0], -0.25));
        assert!(close(g[3][1], 1.0 / 6.0));
        assert!(close(g[3][2], 0.0));
    }
}
