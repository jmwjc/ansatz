//! Linear four-node tetrahedron (`P1`, conforming).

use crate::element::inverse_symmetric_3x3;
use crate::element::{Element, FullIntegration, IntegrationScheme, Keast, ReducedIntegration};
use crate::io::gmsh::{FromGmsh, MeshError, NodeTable};
use crate::node::Node;

/// Four-node linear (conforming) tetrahedron, `Element<D, 3, 4>`.
///
/// The reference tetrahedron is `(0,0,0)`, `(1,0,0)`, `(0,1,0)`, `(0,0,1)`, so the
/// shape functions are the barycentric coordinates
///
/// `N₁ = 1 − ξ − η − ζ`, `N₂ = ξ`, `N₃ = η`, `N₄ = ζ`.
///
/// The gradient and the inverse isoparametric map go through the metric tensor
/// `g_{αβ} = J_α · J_β`, exactly as for [`crate::element::tri3::Tri3`], so the element
/// also works when embedded in a higher-dimensional space (`D > 3`).
#[derive(Clone, Debug, PartialEq)]
pub struct Tet4<const D: usize> {
    nodes: [Node<D>; 4],
}

impl<const D: usize> Tet4<D> {
    /// Builds the element from its four nodes.
    pub fn new(n1: Node<D>, n2: Node<D>, n3: Node<D>, n4: Node<D>) -> Self {
        Tet4 {
            nodes: [n1, n2, n3, n4],
        }
    }

    /// The element's nodes, in local ordering.
    pub fn nodes(&self) -> &[Node<D>; 4] {
        &self.nodes
    }

    /// Node at local index, or `None` if out of bounds.
    pub fn node(&self, index: usize) -> Option<&Node<D>> {
        self.nodes.get(index)
    }

    /// The element's volume.
    pub fn volume(&self) -> f64 {
        self.jacobe([0.0; 3]) / 6.0
    }

    /// Tangent vectors of the affine map, `[x₂ − x₁, x₃ − x₁, x₄ − x₁]`.
    fn tangents(&self) -> [[f64; D]; 3] {
        let a = self.vertices_coordinates();
        std::array::from_fn(|alpha| {
            let b = a[alpha + 1];
            std::array::from_fn(|d| b[d] - a[0][d])
        })
    }

    /// Metric tensor `g_{αβ} = J_α · J_β` packed as
    /// `[g00, g01, g02, g11, g12, g22]`.
    fn metric(&self) -> [f64; 6] {
        let j = self.tangents();
        let dot = |p: [f64; D], q: [f64; D]| -> f64 { (0..D).map(|d| p[d] * q[d]).sum() };
        [
            dot(j[0], j[0]),
            dot(j[0], j[1]),
            dot(j[0], j[2]),
            dot(j[1], j[1]),
            dot(j[1], j[2]),
            dot(j[2], j[2]),
        ]
    }
}

impl<const D: usize> IntegrationScheme<3, 1> for Tet4<D> {
    fn points() -> [[f64; 3]; 1] {
        <Keast<1> as IntegrationScheme<3, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Keast<1> as IntegrationScheme<3, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<3, 4> for Tet4<D> {
    fn points() -> [[f64; 3]; 4] {
        <Keast<4> as IntegrationScheme<3, 4>>::points()
    }
    fn weights() -> [f64; 4] {
        <Keast<4> as IntegrationScheme<3, 4>>::weights()
    }
}

// Full integration: 4-point Keast rule (algebraic degree 2, exact for the quadratic
// mass integrand). Reduced integration: 1 point (centroid), exact for the constant
// stiffness integrand of a linear tetrahedron.
impl<const D: usize> FullIntegration<3, 4> for Tet4<D> {}
impl<const D: usize> ReducedIntegration<3, 1> for Tet4<D> {}

impl<const D: usize> Element<D, 3, 4> for Tet4<D> {
    fn id(&self) -> [usize; 4] {
        [
            self.nodes[0].id,
            self.nodes[1].id,
            self.nodes[2].id,
            self.nodes[3].id,
        ]
    }

    /// Barycentric (area/volume coordinate) shape functions on the reference
    /// tetrahedron `(0,0,0)`, `(1,0,0)`, `(0,1,0)`, `(0,0,1)`.
    fn shape(&self, parametric_coordinates: [f64; 3]) -> [f64; 4] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        let ζ = parametric_coordinates[2];
        [1.0 - ξ - η - ζ, ξ, η, ζ]
    }

    /// Physical gradient `∇ₓN = J_α g^{αβ} ∂_βN`, constant on the element.
    fn derivative_shape(&self, _parametric_coordinates: [f64; 3]) -> [[f64; D]; 4] {
        let j = self.tangents();
        let (inverse, _) = inverse_symmetric_3x3(self.metric());

        // Parametric gradients ∂_β N_a.
        let dn = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        std::array::from_fn(|a| {
            let (dξ, dη, dζ) = (dn[a][0], dn[a][1], dn[a][2]);
            let c0 = inverse[0] * dξ + inverse[1] * dη + inverse[2] * dζ;
            let c1 = inverse[1] * dξ + inverse[3] * dη + inverse[4] * dζ;
            let c2 = inverse[2] * dξ + inverse[4] * dη + inverse[5] * dζ;
            std::array::from_fn(|d| j[0][d] * c0 + j[1][d] * c1 + j[2][d] * c2)
        })
    }

    fn vertices_coordinates(&self) -> [[f64; D]; 4] {
        [
            self.nodes[0].coordinates,
            self.nodes[1].coordinates,
            self.nodes[2].coordinates,
            self.nodes[3].coordinates,
        ]
    }

    /// Isoparametric map `x = Σₐ Nₐ xₐ` (affine for a straight tetrahedron).
    fn coordinates(&self, parametric_coordinates: [f64; 3]) -> [f64; D] {
        let shape = self.shape(parametric_coordinates);
        let vertices = self.vertices_coordinates();
        std::array::from_fn(|d| (0..4).map(|a| shape[a] * vertices[a][d]).sum())
    }

    /// Volume scaling `√det(JᵀJ) = 6 × volume`.
    fn jacobe(&self, _parametric_coordinates: [f64; 3]) -> f64 {
        let (_, det) = inverse_symmetric_3x3(self.metric());
        det.max(0.0).sqrt()
    }

    /// `J[α][d] = ∂x_d/∂ξ_α` (constant for a straight tetrahedron).
    fn jacobe_mat(&self, _parametric_coordinates: [f64; 3]) -> [[f64; D]; 3] {
        self.tangents()
    }

    /// Inverse isoparametric map through the affine frame `[x₂−x₁, x₃−x₁, x₄−x₁]`.
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 3]> {
        let tol = super::CONTAINS_TOL;
        let vertices = self.vertices_coordinates();
        let j = self.tangents();
        let d: [f64; D] = std::array::from_fn(|i| x[i] - vertices[0][i]);

        let metric = self.metric();
        let (inverse, det) = inverse_symmetric_3x3(metric);
        let scale = metric[0].max(metric[3]).max(metric[5]).max(1.0);
        if det.abs() <= 1e-16 * scale * scale * scale {
            return None;
        }

        let dot = |p: [f64; D], q: [f64; D]| -> f64 { (0..D).map(|i| p[i] * q[i]).sum() };
        let rhs = [dot(d, j[0]), dot(d, j[1]), dot(d, j[2])];
        let ξ = inverse[0] * rhs[0] + inverse[1] * rhs[1] + inverse[2] * rhs[2];
        let η = inverse[1] * rhs[0] + inverse[3] * rhs[1] + inverse[4] * rhs[2];
        let ζ = inverse[2] * rhs[0] + inverse[4] * rhs[1] + inverse[5] * rhs[2];

        if 1.0 - ξ - η - ζ < -tol || ξ < -tol || η < -tol || ζ < -tol {
            return None;
        }

        // Out-of-volume residual (only relevant when `D > 3`).
        let residual: f64 = (0..D)
            .map(|i| (d[i] - ξ * j[0][i] - η * j[1][i] - ζ * j[2][i]).powi(2))
            .sum();
        if residual > tol * tol * scale {
            return None;
        }

        Some([ξ, η, ζ])
    }
}

impl<const D: usize> FromGmsh<D> for Tet4<D> {
    const GMSH_ELEMENT_TYPE: i32 = 4;
    const GMSH_NODE_COUNT: usize = 4;

    fn from_gmsh(ids: &[usize], nodes: &NodeTable<D>) -> Result<Self, MeshError> {
        let nodes = [
            Node::new(ids[0], nodes.get(ids[0])?),
            Node::new(ids[1], nodes.get(ids[1])?),
            Node::new(ids[2], nodes.get(ids[2])?),
            Node::new(ids[3], nodes.get(ids[3])?),
        ];
        Ok(Tet4::new(nodes[0], nodes[1], nodes[2], nodes[3]))
    }
}

#[cfg(test)]
mod tests {
    use super::Tet4;
    use crate::element::{Element, IntegrationScheme};
    use crate::io::gmsh::FromGmsh;
    use crate::node::Node;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    /// Skewed tetrahedron (0,0,0), (2,0,0), (0,3,0), (0,0,4): volume 4.
    fn skewed_tet() -> Tet4<3> {
        Tet4::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [0.0, 3.0, 0.0]),
            Node::new(3, [0.0, 0.0, 4.0]),
        )
    }

    #[test]
    fn tet4_shape_and_volume() {
        let tet = skewed_tet();
        assert_eq!(tet.id(), [0, 1, 2, 3]);

        // Kronecker delta at the vertices and partition of unity.
        assert_eq!(tet.shape([0.0, 0.0, 0.0]), [1.0, 0.0, 0.0, 0.0]);
        assert_eq!(tet.shape([1.0, 0.0, 0.0]), [0.0, 1.0, 0.0, 0.0]);
        assert_eq!(tet.shape([0.0, 1.0, 0.0]), [0.0, 0.0, 1.0, 0.0]);
        assert_eq!(tet.shape([0.0, 0.0, 1.0]), [0.0, 0.0, 0.0, 1.0]);
        for xi in [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.3, 0.1]] {
            assert!(close(tet.shape(xi).iter().sum::<f64>(), 1.0));
        }

        // jacobe = 6 × volume = 24.
        assert!(close(tet.volume(), 4.0));
        assert!(close(tet.jacobe([0.0, 0.0, 0.0]), 24.0));

        // The centroid maps to the average of the vertices.
        let centroid = tet.coordinates([0.25, 0.25, 0.25]);
        assert!(close(centroid[0], 0.5));
        assert!(close(centroid[1], 0.75));
        assert!(close(centroid[2], 1.0));
    }

    #[test]
    fn tet4_physical_gradient() {
        // ∇N₁ = (−1/2, −1/3, −1/4), ∇N₂ = (1/2, 0, 0), ∇N₃ = (0, 1/3, 0),
        // ∇N₄ = (0, 0, 1/4); they sum to zero.
        let tet = skewed_tet();
        let g = tet.derivative_shape([0.0, 0.0, 0.0]);
        let expected = [
            [-0.5, -1.0 / 3.0, -0.25],
            [0.5, 0.0, 0.0],
            [0.0, 1.0 / 3.0, 0.0],
            [0.0, 0.0, 0.25],
        ];
        for (a, row) in g.iter().enumerate() {
            for (d, value) in row.iter().enumerate() {
                assert!(close(*value, expected[a][d]), "∇N{a}[{d}]");
            }
        }
        for d in 0..3 {
            assert!(close(g.iter().map(|row| row[d]).sum::<f64>(), 0.0));
        }
    }

    #[test]
    fn tet4_stiffness_and_mass_have_the_expected_values() {
        // Reference tetrahedron: |T| = 1/6, jacobe = 1.
        let reference = Tet4::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [1.0, 0.0, 0.0]),
            Node::new(2, [0.0, 1.0, 0.0]),
            Node::new(3, [0.0, 0.0, 1.0]),
        );
        let jacobe = reference.jacobe([0.0, 0.0, 0.0]);
        let g = reference.derivative_shape([0.0, 0.0, 0.0]);
        let weight = <Tet4<3> as IntegrationScheme<3, 1>>::weights()[0];

        // Kᵢⱼ = |T| ∇Nᵢ·∇Nⱼ with ∫ = |T| (one point, weight 1/6).
        let mut k = [[0.0_f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                let dot: f64 = (0..3).map(|d| g[i][d] * g[j][d]).sum();
                k[i][j] = dot * jacobe * weight;
            }
        }
        assert!(close(k[0][0], 0.5));
        assert!(close(k[0][1], -1.0 / 6.0));
        assert!(close(k[1][1], 1.0 / 6.0));
        assert!(close(k[1][2], 0.0));
        for row in &k {
            assert!(close(row.iter().sum::<f64>(), 0.0));
        }

        // Consistent mass: Mᵢᵢ = |T|/10, Mᵢⱼ = |T|/20.
        let points = <Tet4<3> as IntegrationScheme<3, 4>>::points();
        let weights = <Tet4<3> as IntegrationScheme<3, 4>>::weights();
        let mut m = [[0.0_f64; 4]; 4];
        for g in 0..4 {
            let shape = reference.shape(points[g]);
            for i in 0..4 {
                for j in 0..4 {
                    m[i][j] += shape[i] * shape[j] * jacobe * weights[g];
                }
            }
        }
        for (i, row) in m.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let expected = if i == j { 1.0 / 60.0 } else { 1.0 / 120.0 };
                assert!((value - expected).abs() < 1e-12, "M[{i}][{j}] = {value}");
            }
            assert!((row.iter().sum::<f64>() - 1.0 / 24.0).abs() < 1e-12);
        }
    }

    #[test]
    fn tet4_contains_and_local_coordinates() {
        let tet = skewed_tet();
        assert!(tet.contains([0.3, 0.4, 0.5]));
        let xi = tet.local_coordinates([0.3, 0.4, 0.5]).unwrap();
        assert!(close(xi[0], 0.15) && close(xi[1], 0.4 / 3.0) && close(xi[2], 0.125));
        assert!(!tet.contains([2.0, 3.0, 4.0]));
        assert!(!tet.contains([1.0, 1.0, 1.0]));

        // A spatial tetrahedron embedded in 3D stays consistent with `coordinates`.
        for xi in [[0.1, 0.2, 0.3], [0.6, 0.1, 0.2]] {
            let point = tet.coordinates(xi);
            let back = tet.local_coordinates(point).unwrap();
            for d in 0..3 {
                assert!(close(back[d], xi[d]));
            }
        }
    }

    #[test]
    fn tet4_from_gmsh() {
        use crate::io::gmsh::NodeTable;

        let table = NodeTable::new(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);
        let tet = Tet4::from_gmsh(&[0, 1, 2, 3], &table).unwrap();
        assert_eq!(tet.id(), [0, 1, 2, 3]);
        assert!(close(tet.volume(), 1.0 / 6.0));
    }
}
