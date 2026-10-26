//! Crouzeix–Raviart (CR) elements on tetrahedra: [`CRTet4`] and its boundary trace
//! [`CRTri4`].
//!
//! The 3D counterpart of [`crate::element::crtri3`]; see that module for the CR family
//! overview (facet averages, non-interpolatory basis, `H¹` nonconformity, naming) and
//! `book/src/theory-cr.md` for the theory.
//!
//! The CR tetrahedron keeps the piecewise-linear local space of
//! [`crate::element::tet4::Tet4`] but moves the degrees of freedom from the vertices to
//! the **face averages**: the dof attached to face `f` is `(1/|f|) ∫_f v dA`, which for
//! a linear function equals its value at the face centroid. Basis functions are
//! `θₐ = 1 − 3λₐ`, so `∫_{f_b} θₐ dA = |f_b| δ_ab`.
//!
//! Unlike the 2D case the mass matrix is **not** diagonal: `M_aa = 2|T|/5` and
//! `M_ab = −|T|/20` (`a ≠ b`).
//!
//! Global numbering: the dofs belong to **faces**, not to mesh vertices, so a global
//! face table is needed; [`CRTet4::from_tet4`] builds the element from a
//! [`crate::element::tet4::Tet4`] plus a face-id lookup. Boundary conditions use
//! [`CRTri4`], the trace on one face, whose dofs are the parent's four face dofs.

use std::convert::Infallible;

use crate::element::inverse_symmetric_3x3;
use crate::element::tet4::Tet4;
use crate::element::tri3::Tri3;
use crate::element::{
    Dunavant, Element, FullIntegration, IntegrationScheme, Keast, ReducedIntegration,
};
use crate::node::Node;

/// Crouzeix–Raviart nonconforming linear tetrahedron.
///
/// The 3D counterpart of [`crate::element::crtri3::CRTri3`]: the dofs are the **face averages**, which for a
/// linear function are the values at the face centroids. `nodes[a]` is the centroid of
/// the face opposite vertex `a` and carries the global id of that face's dof, exactly
/// as [`crate::element::crtri3::CRTri3`] stores edge midpoints. Basis functions are
///
/// `θₐ = 1 − 3λₐ`,  so `∫_{f_b} θₐ dA = |f_b| δ_ab`,
///
/// i.e. the global space is not contained in `H¹`: neighbouring tetrahedra only share
/// the *face average*. Unlike the 2D case the mass matrix is **not** diagonal
/// (`M_aa = 2|T|/5`, `M_ab = −|T|/20` for `a ≠ b`).
///
/// Since `Σₐ θₐ Cₐ = Σₐ λₐ Aₐ` for the face centroids `Cₐ` (the barycentric identity
/// of the midpoints), the centroids carry the full affine geometry: the isoparametric
/// map built from them is exact, and [`CRTet4::tetrahedron_vertices`] recovers the
/// vertices by `Aₐ = Σ_b C_b − 3Cₐ`.
///
/// Global face numbering is the caller's business (a tetrahedron has four faces, each
/// shared by two elements); [`CRTet4::from_vertices`] builds the element from four
/// corner nodes plus a face-id lookup, and [`crate::element::crtri3::CRSeg3`]/[`CRTri4`] consume the same ids
/// for boundary conditions. See `book/src/theory-cr.md` for the theory.
#[derive(Clone, Debug, PartialEq)]
pub struct CRTet4<const D: usize> {
    nodes: [Node<D>; 4],
}

impl<const D: usize> CRTet4<D> {
    /// Builds the element from its four face centroids.
    ///
    /// `f1` is the centroid of the face opposite vertex 1, etc.; each node's id is
    /// the global degree of freedom of that face.
    pub fn new(f1: Node<D>, f2: Node<D>, f3: Node<D>, f4: Node<D>) -> Self {
        CRTet4 {
            nodes: [f1, f2, f3, f4],
        }
    }

    /// Builds the element from the four corner nodes of the tetrahedron.
    ///
    /// `face_id` receives the global ids of a face's three corners in **ascending**
    /// order (`[min, mid, max]`) and returns the global id of that face's dof; it is
    /// called four times, once per face, in the order "opposite vertex 0, 1, 2, 3".
    pub fn from_vertices(vertices: [Node<D>; 4], face_id: impl FnMut([usize; 3]) -> usize) -> Self {
        let mut face_id = face_id;
        // `Infallible` makes the fallible path total, so this cannot fail.
        match Self::try_from_vertices(vertices, |face| Ok::<usize, Infallible>(face_id(face))) {
            Ok(element) => element,
            Err(never) => match never {},
        }
    }

    /// Fallible variant of [`CRTet4::from_vertices`], for face tables that only know
    /// already-numbered faces (e.g. a numbering created while reading a mesh).
    ///
    /// # Errors
    ///
    /// Returns whatever error `face_id` produces; the element is not built.
    pub fn try_from_vertices<E>(
        vertices: [Node<D>; 4],
        mut face_id: impl FnMut([usize; 3]) -> Result<usize, E>,
    ) -> Result<Self, E> {
        // Face `a` is opposite local vertex `a`, i.e. the other three vertices.
        let local_faces = [(1_usize, 2_usize, 3_usize), (0, 2, 3), (0, 1, 3), (0, 1, 2)];
        let mut nodes = [Node::new(0, [0.0; D]); 4];
        for (a, (i, j, k)) in local_faces.into_iter().enumerate() {
            let mut ids = [vertices[i].id, vertices[j].id, vertices[k].id];
            ids.sort_unstable();
            let coordinates = std::array::from_fn(|d| {
                (vertices[i].coordinates[d]
                    + vertices[j].coordinates[d]
                    + vertices[k].coordinates[d])
                    / 3.0
            });
            nodes[a] = Node::new(face_id(ids)?, coordinates);
        }
        Ok(CRTet4 { nodes })
    }

    /// Builds the CR element of the tetrahedron `tet`.
    ///
    /// `face_id` receives the global ids of a face's three corners in **ascending**
    /// order (`[min, mid, max]`) and returns the global id of that face's dof; it is
    /// called once per face, in the order "opposite vertex 0, 1, 2, 3".
    pub fn from_tet4(tet: &Tet4<D>, face_id: impl FnMut([usize; 3]) -> usize) -> Self {
        Self::from_vertices(tet_corners(tet), face_id)
    }

    /// Fallible variant of [`CRTet4::from_tet4`], for face tables that only know
    /// already-numbered faces.
    ///
    /// # Errors
    ///
    /// Returns whatever error `face_id` produces; the element is not built.
    pub fn try_from_tet4<E>(
        tet: &Tet4<D>,
        face_id: impl FnMut([usize; 3]) -> Result<usize, E>,
    ) -> Result<Self, E> {
        Self::try_from_vertices(tet_corners(tet), face_id)
    }

    /// The element's nodes, i.e. the four face centroids, in local ordering.
    pub fn nodes(&self) -> &[Node<D>; 4] {
        &self.nodes
    }

    /// Node at local index, or `None` if out of bounds.
    pub fn node(&self, index: usize) -> Option<&Node<D>> {
        self.nodes.get(index)
    }

    /// Tetrahedron vertices recovered from the face centroids,
    /// `Aₐ = Σ_b C_b − 3 Cₐ`.
    pub fn tetrahedron_vertices(&self) -> [[f64; D]; 4] {
        let centroids = self.vertices_coordinates();
        let sum: [f64; D] = std::array::from_fn(|d| centroids.iter().map(|c| c[d]).sum());
        std::array::from_fn(|a| std::array::from_fn(|d| sum[d] - 3.0 * centroids[a][d]))
    }

    /// The element's volume.
    pub fn volume(&self) -> f64 {
        self.jacobe([0.0; 3]) / 6.0
    }

    /// Tangent vectors of the affine map, `[A₂ − A₁, A₃ − A₁, A₄ − A₁]`.
    fn tangents(&self) -> [[f64; D]; 3] {
        let a = self.tetrahedron_vertices();
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

/// The four corner nodes of a [`Tet4`], in local order.
fn tet_corners<const D: usize>(tet: &Tet4<D>) -> [Node<D>; 4] {
    let id = tet.id();
    let vertices = tet.vertices_coordinates();
    std::array::from_fn(|a| Node::new(id[a], vertices[a]))
}

impl<const D: usize> IntegrationScheme<3, 1> for CRTet4<D> {
    fn points() -> [[f64; 3]; 1] {
        <Keast<1> as IntegrationScheme<3, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Keast<1> as IntegrationScheme<3, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<3, 4> for CRTet4<D> {
    fn points() -> [[f64; 3]; 4] {
        <Keast<4> as IntegrationScheme<3, 4>>::points()
    }
    fn weights() -> [f64; 4] {
        <Keast<4> as IntegrationScheme<3, 4>>::weights()
    }
}

impl<const D: usize> FullIntegration<3, 4> for CRTet4<D> {}
impl<const D: usize> ReducedIntegration<3, 1> for CRTet4<D> {}

impl<const D: usize> Element<D, 3, 4> for CRTet4<D> {
    fn id(&self) -> [usize; 4] {
        [
            self.nodes[0].id,
            self.nodes[1].id,
            self.nodes[2].id,
            self.nodes[3].id,
        ]
    }

    /// Face-average basis `θₐ = 1 − 3λₐ` on the reference tetrahedron
    /// `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`, with `λ = (1 − ξ − η − ζ, ξ, η, ζ)`.
    fn shape(&self, parametric_coordinates: [f64; 3]) -> [f64; 4] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        let ζ = parametric_coordinates[2];
        let lambda = [1.0 - ξ - η - ζ, ξ, η, ζ];
        std::array::from_fn(|a| 1.0 - 3.0 * lambda[a])
    }

    /// Physical gradient `∇ₓθₐ = −3 ∇ₓλₐ`, constant on the element and computed with
    /// the metric tensor: `∇ₓN = J_α g^{αβ} ∂_βN`.
    fn derivative_shape(&self, _parametric_coordinates: [f64; 3]) -> [[f64; D]; 4] {
        let j = self.tangents();
        let (inverse, _) = inverse_symmetric_3x3(self.metric());

        // Parametric gradients ∂_β N_a of `θₐ = 1 − 3λₐ`.
        let dn = [
            [3.0, 3.0, 3.0],
            [-3.0, 0.0, 0.0],
            [0.0, -3.0, 0.0],
            [0.0, 0.0, -3.0],
        ];
        std::array::from_fn(|a| {
            let (dξ, dη, dζ) = (dn[a][0], dn[a][1], dn[a][2]);
            let c0 = inverse[0] * dξ + inverse[1] * dη + inverse[2] * dζ;
            let c1 = inverse[1] * dξ + inverse[3] * dη + inverse[4] * dζ;
            let c2 = inverse[2] * dξ + inverse[4] * dη + inverse[5] * dζ;
            std::array::from_fn(|d| j[0][d] * c0 + j[1][d] * c1 + j[2][d] * c2)
        })
    }

    /// The four face centroids (the element's nodes).
    fn vertices_coordinates(&self) -> [[f64; D]; 4] {
        [
            self.nodes[0].coordinates,
            self.nodes[1].coordinates,
            self.nodes[2].coordinates,
            self.nodes[3].coordinates,
        ]
    }

    /// Isoparametric map built from the face centroids: `x = Σₐ θₐ Cₐ`, which equals
    /// the affine map `Σₐ λₐ Aₐ` of the tetrahedron.
    fn coordinates(&self, parametric_coordinates: [f64; 3]) -> [f64; D] {
        let shape = self.shape(parametric_coordinates);
        let centroids = self.vertices_coordinates();
        std::array::from_fn(|d| (0..4).map(|a| shape[a] * centroids[a][d]).sum())
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

    /// Inverse isoparametric map through the affine frame `[A₂−A₁, A₃−A₁, A₄−A₁]`.
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 3]> {
        let tol = crate::element::CONTAINS_TOL;
        let vertices = self.tetrahedron_vertices();
        let j = self.tangents();
        let d: [f64; D] = std::array::from_fn(|i| x[i] - vertices[0][i]);

        let dot = |p: [f64; D], q: [f64; D]| -> f64 { (0..D).map(|i| p[i] * q[i]).sum() };
        let metric = self.metric();
        if metric.iter().all(|g| *g == 0.0) {
            return None;
        }
        let (inverse, det) = inverse_symmetric_3x3(metric);
        let scale = metric[0].max(metric[3]).max(metric[5]).max(1.0);
        if det.abs() <= 1e-16 * scale * scale * scale {
            return None;
        }

        let rhs = [dot(d, j[0]), dot(d, j[1]), dot(d, j[2])];
        let ξ = inverse[0] * rhs[0] + inverse[1] * rhs[1] + inverse[2] * rhs[2];
        let η = inverse[1] * rhs[0] + inverse[3] * rhs[1] + inverse[4] * rhs[2];
        let ζ = inverse[2] * rhs[0] + inverse[4] * rhs[1] + inverse[5] * rhs[2];

        if 1.0 - ξ - η - ζ < -tol || ξ < -tol || η < -tol || ζ < -tol {
            return None;
        }

        // Out-of-plane residual (only relevant when `D > 3`).
        let residual: f64 = (0..D)
            .map(|i| (d[i] - ξ * j[0][i] - η * j[1][i] - ζ * j[2][i]).powi(2))
            .sum();
        if residual > tol * tol * scale {
            return None;
        }

        Some([ξ, η, ζ])
    }
}

/// Trace of a CR tetrahedron on one of its four facets: the CR boundary triangle.
///
/// This is the 3D counterpart of [`crate::element::crtri3::CRSeg3`]: the dofs are the **four** facet dofs of
/// the parent tetrahedron (a CR tetrahedron has no vertex dofs), the shape functions
/// are the restrictions `θₐ|face = 1 − 3λₐ`, and `facet = k` again means "the facet
/// opposite parent vertex `k`", with the facet's vertices listed in ascending parent
/// index order.
///
/// The parent CR tetrahedron type is not implemented yet, so this element is
/// currently constructed directly with [`CRTri4::new`]; it exists so that 3D
/// boundary handling can be written (and tested) dimension-agnostically.
#[derive(Clone, Debug, PartialEq)]
pub struct CRTri4<const D: usize> {
    facets: [usize; 4],
    vertices: [Node<D>; 4],
    facet: usize,
}

impl<const D: usize> CRTri4<D> {
    /// Builds a facet element from the parent's four facet dofs, the parent's four
    /// vertices (local order) and the local facet index (`0..4`).
    ///
    /// Returns `None` if `facet` is not in `0..4`.
    pub fn new(facets: [usize; 4], vertices: [Node<D>; 4], facet: usize) -> Option<Self> {
        (facet < 4).then_some(CRTri4 {
            facets,
            vertices,
            facet,
        })
    }

    /// Traces the CR tetrahedron `parent` on its local facet `facet`.
    ///
    /// Returns `None` if `facet` is not in `0..4`.
    pub fn from_crtet4(parent: &CRTet4<D>, facet: usize) -> Option<Self> {
        let vertices = parent.tetrahedron_vertices();
        Self::new(
            parent.id(),
            std::array::from_fn(|a| Node::new(a, vertices[a])),
            facet,
        )
    }

    /// The local facet index: the facet opposite parent vertex `local_facet()`.
    pub fn local_facet(&self) -> usize {
        self.facet
    }

    /// The parent element's vertices, in local order.
    pub fn parent_vertices(&self) -> &[Node<D>; 4] {
        &self.vertices
    }

    /// The facet's three vertices, in ascending parent index order.
    pub fn facet_vertices(&self) -> [Node<D>; 3] {
        let mut out = [self.vertices[0]; 3];
        let mut next = 0;
        for (index, vertex) in self.vertices.iter().enumerate() {
            if index != self.facet {
                out[next] = *vertex;
                next += 1;
            }
        }
        out
    }

    /// The facet geometry as a linear triangle.
    fn as_tri3(&self) -> Tri3<D> {
        let [a, b, c] = self.facet_vertices();
        Tri3::new(a, b, c)
    }
}

impl<const D: usize> Element<D, 2, 4> for CRTri4<D> {
    fn id(&self) -> [usize; 4] {
        self.facets
    }

    /// Trace basis on the facet: `shapeₐ = 1 − 3λₐ` with the parent's `λ` scattered
    /// from the facet's barycentric coordinates (the opposite vertex is zero).
    /// At the facet centroid the facet's own dof is `1` and the other three are `0`.
    fn shape(&self, parametric_coordinates: [f64; 2]) -> [f64; 4] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        let facet_lambda = [1.0 - ξ - η, ξ, η];

        let mut lambda = [0.0_f64; 4];
        let mut next = 0;
        for (index, value) in lambda.iter_mut().enumerate() {
            if index != self.facet {
                *value = facet_lambda[next];
                next += 1;
            }
        }

        std::array::from_fn(|a| 1.0 - 3.0 * lambda[a])
    }

    /// In-facet (tangential) gradient of the trace: `−3 ∇λₐ`, obtained from the
    /// facet's [`Tri3`] gradients. Boundary operators do not use it.
    fn derivative_shape(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 4] {
        let grad = self.as_tri3().derivative_shape(parametric_coordinates);
        let mut out = [[0.0_f64; D]; 4];
        let mut next = 0;
        for (index, row) in out.iter_mut().enumerate() {
            if index != self.facet {
                *row = std::array::from_fn(|d| -3.0 * grad[next][d]);
                next += 1;
            }
        }
        out
    }

    /// The parent's vertices (the element's nodes are facets, not points).
    fn vertices_coordinates(&self) -> [[f64; D]; 4] {
        self.vertices.map(|vertex| vertex.coordinates)
    }

    fn coordinates(&self, parametric_coordinates: [f64; 2]) -> [f64; D] {
        self.as_tri3().coordinates(parametric_coordinates)
    }

    /// Area scaling of the facet: `2 × area`.
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

impl<const D: usize> IntegrationScheme<2, 1> for CRTri4<D> {
    fn points() -> [[f64; 2]; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<2, 3> for CRTri4<D> {
    fn points() -> [[f64; 2]; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::points()
    }
    fn weights() -> [f64; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::weights()
    }
}

impl<const D: usize> FullIntegration<2, 3> for CRTri4<D> {}
impl<const D: usize> ReducedIntegration<2, 1> for CRTri4<D> {}

#[cfg(test)]
mod tests {
    use super::{CRTet4, CRTri4};
    use crate::element::tet4::Tet4;
    use crate::element::{Element, IntegrationScheme};
    use crate::node::Node;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    /// Skewed tetrahedron (0,0,0), (2,0,0), (0,3,0), (0,0,4): volume 4.
    fn skewed_tet4() -> Tet4<3> {
        Tet4::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [2.0, 0.0, 0.0]),
            Node::new(2, [0.0, 3.0, 0.0]),
            Node::new(3, [0.0, 0.0, 4.0]),
        )
    }

    /// The same tetrahedron as a `CRTet4`: `jacobe = 24`; the face dof ids are derived
    /// from the ascending face vertex ids so that they are easy to inspect.
    fn skewed_tet() -> CRTet4<3> {
        CRTet4::from_tet4(&skewed_tet4(), |face| {
            10 + face[0] * 100 + face[1] * 10 + face[2]
        })
    }

    #[test]
    fn crtet4_from_tet4_uses_the_face_centroids() {
        let tet = skewed_tet4();
        let cr = skewed_tet();
        let vertices = tet.vertices_coordinates();
        let local_faces = [[1_usize, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for (a, face) in local_faces.iter().enumerate() {
            let expected: [f64; 3] = std::array::from_fn(|d| {
                (vertices[face[0]][d] + vertices[face[1]][d] + vertices[face[2]][d]) / 3.0
            });
            assert_eq!(cr.nodes()[a].coordinates, expected, "face {a} centroid");
        }
        // The face dof ids come from the ascending vertex triple of the face, which
        // for face 0 (opposite vertex 0) is made of vertices 1, 2 and 3.
        let mut ids = [tet.id()[1], tet.id()[2], tet.id()[3]];
        ids.sort_unstable();
        assert_eq!(cr.id()[0], 10 + ids[0] * 100 + ids[1] * 10 + ids[2]);
        assert_eq!(cr.tetrahedron_vertices(), vertices);
    }

    #[test]
    fn crtet4_shape_is_delta_on_face_centroids() {
        let tet = skewed_tet();

        // The centroids, in parametric coordinates, are the reference face centroids.
        let centroids = [
            [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
            [0.0, 1.0 / 3.0, 1.0 / 3.0],
            [1.0 / 3.0, 0.0, 1.0 / 3.0],
            [1.0 / 3.0, 1.0 / 3.0, 0.0],
        ];
        for (b, xi) in centroids.iter().enumerate() {
            let shape = tet.shape(*xi);
            for (a, value) in shape.iter().enumerate() {
                let expected = if a == b { 1.0 } else { 0.0 };
                assert!(
                    (value - expected).abs() < 1e-12,
                    "θ{a}(C{b}) = {value} != {expected}"
                );
            }
        }

        // At the vertices θₐ = 1 − 3λₐ, and the basis partitions unity.
        assert_eq!(tet.shape([0.0, 0.0, 0.0]), [-2.0, 1.0, 1.0, 1.0]);
        assert_eq!(tet.shape([1.0, 0.0, 0.0]), [1.0, -2.0, 1.0, 1.0]);
        for xi in [[0.0, 0.0, 0.0], [0.5, 0.25, 0.125], [0.25, 0.25, 0.25]] {
            assert!(close(tet.shape(xi).iter().sum::<f64>(), 1.0));
        }
    }

    #[test]
    fn crtet4_geometry_matches_the_tetrahedron() {
        let tet = skewed_tet();
        let vertices = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ];

        assert_eq!(tet.tetrahedron_vertices(), vertices);
        assert!(close(tet.volume(), 4.0));
        assert!(close(tet.jacobe([0.0, 0.0, 0.0]), 24.0));

        // The isoparametric map from the centroids is the affine map.
        for xi in [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.6, 0.2, 0.1]] {
            let point = tet.coordinates(xi);
            let lambda = [1.0 - xi[0] - xi[1] - xi[2], xi[0], xi[1], xi[2]];
            for (d, value) in point.iter().enumerate() {
                let expected: f64 = (0..4).map(|a| lambda[a] * vertices[a][d]).sum();
                assert!(close(*value, expected));
            }
        }

        // ∇θₐ = −3∇λₐ with ∇λ = (−1/2, −1/3, −1/4), (1/2, 0, 0), (0, 1/3, 0), (0, 0, 1/4).
        let g = tet.derivative_shape([0.0, 0.0, 0.0]);
        let expected = [
            [1.5, 1.0, 0.75],
            [-1.5, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -0.75],
        ];
        for (a, row) in g.iter().enumerate() {
            for (d, value) in row.iter().enumerate() {
                assert!(close(*value, expected[a][d]), "∇θ{a}[{d}]");
            }
        }
        for d in 0..3 {
            let sum: f64 = g.iter().map(|row| row[d]).sum();
            assert!(close(sum, 0.0));
        }
    }

    #[test]
    fn crtet4_stiffness_is_the_scaled_barycentric_form() {
        // Reference tetrahedron: |T| = 1/6, ∇λ₁ = (−1,−1,−1) etc., so
        // K = 9|T| ∇λᵢ·∇λⱼ = 1.5 × (row sums vanish).
        let reference = CRTet4::new(
            Node::new(0, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]),
            Node::new(1, [0.0, 1.0 / 3.0, 1.0 / 3.0]),
            Node::new(2, [1.0 / 3.0, 0.0, 1.0 / 3.0]),
            Node::new(3, [1.0 / 3.0, 1.0 / 3.0, 0.0]),
        );
        let jacobe = reference.jacobe([0.0, 0.0, 0.0]);
        assert!(close(jacobe, 1.0), "the reference tetrahedron has 6|T| = 1");

        let g = reference.derivative_shape([0.0, 0.0, 0.0]);
        let mut k = [[0.0_f64; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                let dot: f64 = (0..3).map(|d| g[i][d] * g[j][d]).sum();
                k[i][j] = dot * jacobe * <CRTet4<3> as IntegrationScheme<3, 1>>::weights()[0];
            }
        }
        assert!(close(k[0][0], 4.5));
        assert!(close(k[0][1], -1.5));
        assert!(close(k[0][2], -1.5));
        assert!(close(k[0][3], -1.5));
        assert!(close(k[1][1], 1.5));
        assert!(close(k[1][2], 0.0));
        assert!(close(k[2][2], 1.5));
        assert!(close(k[3][3], 1.5));
        for row in &k {
            assert!(close(row.iter().sum::<f64>(), 0.0));
        }
    }

    #[test]
    fn crtet4_mass_matrix_is_not_diagonal() {
        // ∫θₐ² = 2|T|/5, ∫θₐθ_b = −|T|/20 (a ≠ b); row sum = |T|/4 = 1/24.
        let reference = CRTet4::new(
            Node::new(0, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]),
            Node::new(1, [0.0, 1.0 / 3.0, 1.0 / 3.0]),
            Node::new(2, [1.0 / 3.0, 0.0, 1.0 / 3.0]),
            Node::new(3, [1.0 / 3.0, 1.0 / 3.0, 0.0]),
        );
        let volume = 1.0 / 6.0;
        let points = <CRTet4<3> as IntegrationScheme<3, 4>>::points();
        let weights = <CRTet4<3> as IntegrationScheme<3, 4>>::weights();
        let jacobe = reference.jacobe([0.0, 0.0, 0.0]);

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
                let expected = if i == j {
                    2.0 * volume / 5.0
                } else {
                    -volume / 20.0
                };
                assert!(
                    (value - expected).abs() < 1e-12,
                    "M[{i}][{j}] = {value} != {expected}"
                );
            }
            assert!((row.iter().sum::<f64>() - volume / 4.0).abs() < 1e-12);
        }
    }

    #[test]
    fn crtet4_contains_and_local_coordinates() {
        let tet = skewed_tet();
        assert!(tet.contains([0.3, 0.4, 0.5]));
        let xi = tet.local_coordinates([0.3, 0.4, 0.5]).unwrap();
        assert!(close(xi[0], 0.15) && close(xi[1], 0.4 / 3.0) && close(xi[2], 0.125));
        assert!(!tet.contains([2.0, 3.0, 4.0]));
        assert!(!tet.contains([1.0, 1.0, 1.0]));
    }

    #[test]
    fn crtri4_matches_the_parent_tet_trace() {
        // The face trace's basis equals the parent CRTet4 basis at the corresponding
        // barycentric coordinates.
        let tet = skewed_tet();
        for facet in 0..4 {
            let trace = CRTri4::from_crtet4(&tet, facet).unwrap();
            assert_eq!(trace.id(), tet.id(), "the dofs are the parent's face dofs");
            for xi in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 0.2]] {
                let point = trace.coordinates(xi);
                let lambda = tet.local_coordinates(point).unwrap();
                let parent = tet.shape(lambda);
                let shape = trace.shape(xi);
                for a in 0..4 {
                    assert!(
                        (shape[a] - parent[a]).abs() < 1e-12,
                        "facet {facet}, ξ = {xi:?}, basis {a}: {} != {}",
                        shape[a],
                        parent[a]
                    );
                }
            }
        }
        assert!(CRTri4::from_crtet4(&tet, 4).is_none());
    }

    #[test]
    fn crtri4_trace_basis_and_geometry() {
        // Parent tetrahedron (0,0,0), (1,0,0), (0,1,0), (0,0,1), traced on each of
        // its four facets: `shapeₐ = 1 − 3λₐ` with the parent's λ (the opposite
        // vertex has λ = 0).
        let vertices = [
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [1.0, 0.0, 0.0]),
            Node::new(2, [0.0, 1.0, 0.0]),
            Node::new(3, [0.0, 0.0, 1.0]),
        ];
        let facet_area = |facet: usize| {
            if facet == 0 {
                3.0_f64.sqrt() / 2.0
            } else {
                0.5
            }
        };

        for facet in 0..4 {
            let trace = CRTri4::new([10, 11, 12, 13], vertices, facet).unwrap();
            assert_eq!(
                trace.id(),
                [10, 11, 12, 13],
                "the dofs are the parent's facets"
            );
            assert_eq!(trace.local_facet(), facet);

            // At a facet vertex with parent index `vertex` the parent λ is a unit
            // vector, so `shape = 1 − 3` there and `1` elsewhere (and `1` at the
            // facet's own dof).
            let ends = trace.facet_vertices();
            for end in &ends {
                let vertex = vertices.iter().position(|v| v == end).unwrap();
                let local = trace.local_coordinates(end.coordinates).unwrap();
                let shape = trace.shape(local);
                for (a, value) in shape.iter().enumerate() {
                    let expected = if a == vertex { -2.0 } else { 1.0 };
                    assert!(
                        (value - expected).abs() < 1e-12,
                        "facet {facet}, vertex {vertex}, basis {a}: {value} != {expected}"
                    );
                }
            }

            // At the facet centroid the trace is `1` for the facet's own dof and `0`
            // for the other three (this is where the parent's dof lives).
            let centroid = trace.shape([1.0 / 3.0, 1.0 / 3.0]);
            for (a, value) in centroid.iter().enumerate() {
                let expected = if a == facet { 1.0 } else { 0.0 };
                assert!(
                    (value - expected).abs() < 1e-12,
                    "facet {facet}, centroid, basis {a}: {value} != {expected}"
                );
            }

            // Partition of unity everywhere.
            for xi in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.3, 0.2]] {
                assert!(close(trace.shape(xi).iter().sum::<f64>(), 1.0));
            }

            // jacobe = 2 × facet area.
            let jacobe = trace.jacobe([1.0 / 3.0, 1.0 / 3.0]);
            assert!(
                (jacobe - 2.0 * facet_area(facet)).abs() < 1e-12,
                "facet {facet}: {jacobe} != {}",
                2.0 * facet_area(facet)
            );
        }

        // Out-of-range facet indices are rejected.
        assert!(CRTri4::new([10, 11, 12, 13], vertices, 4).is_none());
    }
}
