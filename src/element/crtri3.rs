//! Crouzeix–Raviart (CR) elements on triangles: [`CRTri3`] and its boundary trace
//! [`CRSeg3`].
//!
//! The CR family keeps the piecewise-linear local space of the conforming simplices
//! ([`crate::element::tri3::Tri3`], [`crate::element::tet4::Tet4`]) but moves the
//! degrees of freedom from the vertices to the **facet averages**: the dof attached to
//! facet `f` is `(1/|f|) ∫_f v ds`, which for a linear function equals its value at the
//! facet's centroid. Basis functions are `θₐ = 1 − D λₐ` (`D` the space dimension, `λₐ`
//! the barycentric coordinates), so `∫_{f_b} θₐ ds = |f_b| δ_ab`.
//!
//! Consequently the elements are **not interpolatory** (`θₐ` takes `−D` at its own
//! vertex and `+1` at the others) and a global CR space is not a subspace of `H¹`:
//! neighbouring elements only share the *facet averages* of the field. On the other
//! hand the 2D mass matrix is exactly diagonal (`M = (|T|/3) I`), and the classic uses
//! are incompressible (Stokes, nearly incompressible) problems. See
//! `book/src/theory-cr.md` for the theory.
//!
//! # Contents
//!
//! An interior element and its boundary trace live in one module, because a trace
//! carries the *whole* parent element (the field on a facet is a linear function of all
//! of the parent's facet dofs, not only of the facet's own dof):
//!
//! - this module: [`CRTri3`] (triangle, 3 edge dofs) with its trace [`CRSeg3`] (a
//!   segment carrying the parent's three edge dofs);
//! - [`crate::element::crtet4`]: [`crate::element::crtet4::CRTet4`] (tetrahedron, 4 face
//!   dofs) with its trace [`crate::element::crtet4::CRTri4`] (a triangle carrying the
//!   parent's four face dofs).
//!
//! # Naming
//!
//! Names follow the crate-wide `<shape><dof-count>` convention (`Seg2`, `Tri3`,
//! `Quad4`), where the shape is the element's **own** shape and the number counts its
//! degrees of freedom. Hence `Tri` appears in the 3D trace as well — a face of a
//! tetrahedron is a triangle — while `Tet` names the 3D *interior* element:
//!
//! | Type | Shape | Dofs | Role |
//! |---|---|---|---|
//! | [`CRTri3`] | triangle | 3 (its edges) | 2D interior element |
//! | [`CRSeg3`] | segment | 3 (inherited from the parent triangle) | 2D boundary trace |
//! | [`CRTet4`](crate::element::crtet4::CRTet4) | tetrahedron | 4 (its faces) | 3D interior element |
//! | [`CRTri4`](crate::element::crtet4::CRTri4) | triangle | 4 (inherited from the parent tetrahedron) | 3D boundary trace |
//!
//! In every CR element the number of dofs equals the number of facets of the parent
//! simplex, which is why the traces reuse `3` (2D) and `4` (3D).
//!
//! A dof belongs to a **facet**, not to a mesh vertex, so the global numbering needs a
//! facet table. `io::gmsh` builds one while reading
//! ([`crate::io::gmsh::GmshOptions`]) and exposes the elements through
//! `crtri3_elements`, `crtri3_boundary`, `crtet4_elements` and `crtet4_boundary`.
//!
//! # 2D specifics
//!
//! The CR triangle keeps the piecewise-linear local space of [`crate::element::tri3::Tri3`]
//! but moves the degrees of freedom from the vertices to the **edge averages**: the dof
//! attached to edge `eᵢ` is `(1/|eᵢ|) ∫_{eᵢ} v ds`, which for a linear function equals
//! its value at the edge midpoint. Basis functions are `θᵢ = 1 − 2λᵢ`, so
//! `∫_{eⱼ} θᵢ ds = |eⱼ| δᵢⱼ`; the mass matrix is `M = (|T|/3) I`.
//!
//! Global numbering: the dofs belong to **edges**. A CR element is therefore built
//! from three edge midpoints, each of which carries the global edge-dof id in its
//! [`Node::id`]; [`CRTri3::from_tri3`] derives those ids (and the midpoints) from a
//! [`crate::element::tri3::Tri3`] plus an edge-id lookup. Boundary conditions use
//! [`CRSeg3`], the trace on one edge, whose dofs are the parent's three edge dofs.

use std::convert::Infallible;

use crate::element::seg2::Seg2;
use crate::element::tri3::Tri3;
use crate::element::{
    Dunavant, Element, ElementBase, FullIntegration, GaussLegendre1D, IntegrationScheme,
    ReducedIntegration,
};
use crate::node::Node;

/// Crouzeix–Raviart linear triangle with one edge-average dof per edge.
///
/// `nodes[i]` is the **midpoint of the edge opposite vertex `i`**, together with
/// the global id of that edge's dof. The reference element is the unit triangle
/// `(0, 0)`, `(1, 0)`, `(0, 1)`, so vertex 1 sits at `ξ = (0, 0)` and the edge
/// opposite it is spanned by the other two vertices.
#[derive(Clone, Debug, PartialEq)]
pub struct CRTri3<const D: usize> {
    nodes: [Node<D>; 3],
}

impl<const D: usize> CRTri3<D> {
    /// Builds the element from its three edge midpoints.
    ///
    /// `m1` is the midpoint of the edge opposite vertex 1, `m2` of the edge
    /// opposite vertex 2, `m3` of the edge opposite vertex 3. Each node's id is
    /// the global degree of freedom of that edge.
    pub fn new(m1: Node<D>, m2: Node<D>, m3: Node<D>) -> Self {
        CRTri3 {
            nodes: [m1, m2, m3],
        }
    }

    /// Builds the CR element of the triangle `tri`.
    ///
    /// `edge_id` maps an edge, given by its two end vertex ids in ascending order,
    /// to the global id of that edge's dof. It is called three times per element
    /// (once for each edge), so callers normally pass a closure that looks the edge
    /// up in — or inserts it into — an edge table.
    pub fn from_tri3(tri: &Tri3<D>, edge_id: impl FnMut(usize, usize) -> usize) -> Self {
        let mut edge_id = edge_id;
        // `Infallible` makes the fallible path total, so this cannot fail.
        let built = Self::try_from_tri3(tri, |a, b| Ok::<usize, Infallible>(edge_id(a, b)));
        match built {
            Ok(cr) => cr,
            Err(never) => match never {},
        }
    }

    /// Fallible variant of [`CRTri3::from_tri3`], for edge tables that only know
    /// already-numbered edges (e.g. a numbering created while reading a mesh).
    ///
    /// # Errors
    ///
    /// Returns whatever error `edge_id` produces; the element is not built.
    pub fn try_from_tri3<E>(
        tri: &Tri3<D>,
        mut edge_id: impl FnMut(usize, usize) -> Result<usize, E>,
    ) -> Result<Self, E> {
        let id = tri.id();
        let vertices = tri.vertices_coordinates();

        // Edge `i` is the edge opposite local vertex `i`, i.e. the pair of the other
        // two local vertices.
        let edges = [(1_usize, 2_usize), (0, 2), (0, 1)];
        let mut nodes = [Node::new(0, [0.0; D]); 3];
        for (i, (a, b)) in edges.into_iter().enumerate() {
            let (lo, hi) = if id[a] <= id[b] {
                (id[a], id[b])
            } else {
                (id[b], id[a])
            };
            let coordinates = std::array::from_fn(|k| 0.5 * (vertices[a][k] + vertices[b][k]));
            nodes[i] = Node::new(edge_id(lo, hi)?, coordinates);
        }

        Ok(CRTri3 { nodes })
    }

    /// The element's nodes, i.e. the three edge midpoints, in local ordering.
    pub fn nodes(&self) -> &[Node<D>; 3] {
        &self.nodes
    }

    /// Node at local index, or `None` if out of bounds.
    pub fn node(&self, index: usize) -> Option<&Node<D>> {
        self.nodes.get(index)
    }

    /// Triangle vertices recovered from the midpoints by reflection,
    /// `Aᵢ = Mⱼ + Mₖ − Mᵢ`.
    ///
    /// Since `Σᵢ θᵢ Mᵢ = Σᵢ λᵢ Aᵢ`, the midpoints carry the full affine geometry;
    /// this method (`× example`: for plotting or mesh-size estimates) exposes the
    /// underlying triangle explicitly.
    pub fn triangle_vertices(&self) -> [[f64; D]; 3] {
        let [m1, m2, m3] = self.vertices_coordinates();
        let reflect =
            |a: [f64; D], b: [f64; D], c: [f64; D]| std::array::from_fn(|i| b[i] + c[i] - a[i]);
        [
            reflect(m1, m2, m3),
            reflect(m2, m3, m1),
            reflect(m3, m1, m2),
        ]
    }

    /// The underlying linear triangle, used for the geometry (Jacobian, inverse
    /// map) and for the barycentric-coordinate gradients.
    fn as_tri3(&self) -> Tri3<D> {
        let [a, b, c] = self.triangle_vertices();
        Tri3::new(Node::new(0, a), Node::new(1, b), Node::new(2, c))
    }
}

// Full integration: 3-point Dunavant rule (algebraic degree 2), needed for the
// quadratic mass integrand θᵢθⱼ. Reduced integration: 1 point (centroid), which
// is already exact for the constant stiffness integrand ∇θᵢ·∇θⱼ.
impl<const D: usize> IntegrationScheme<2, 1> for CRTri3<D> {
    fn points() -> [[f64; 2]; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<2, 3> for CRTri3<D> {
    fn points() -> [[f64; 2]; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::points()
    }
    fn weights() -> [f64; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::weights()
    }
}

impl<const D: usize> FullIntegration<2, 3> for CRTri3<D> {}
impl<const D: usize> ReducedIntegration<2, 1> for CRTri3<D> {}

impl<const D: usize> ElementBase<D, 2, 3> for CRTri3<D> {
    fn id(&self) -> [usize; 3] {
        [self.nodes[0].id, self.nodes[1].id, self.nodes[2].id]
    }

    /// Coordinates of the three edge midpoints (the element's nodes).
    ///
    /// Note that these are *not* the triangle vertices; use
    /// [`CRTri3::triangle_vertices`] for those. The midpoints still reproduce the
    /// affine geometry, because `Σᵢ θᵢ Mᵢ = Σᵢ λᵢ Aᵢ`.
    fn vertices_coordinates(&self) -> [[f64; D]; 3] {
        [
            self.nodes[0].coordinates,
            self.nodes[1].coordinates,
            self.nodes[2].coordinates,
        ]
    }

    /// Isoparametric map built from the midpoints: `x(ξ) = Σᵢ θᵢ(ξ) Mᵢ`.
    ///
    /// This is the affine map of the underlying triangle, so it is exact.
    fn coordinates(&self, parametric_coordinates: [f64; 2]) -> [f64; D] {
        let shape = self.shape(parametric_coordinates);
        let midpoints = self.vertices_coordinates();
        std::array::from_fn(|i| {
            midpoints[0][i] * shape[0] + midpoints[1][i] * shape[1] + midpoints[2][i] * shape[2]
        })
    }

    /// Area scaling `√det(JᵀJ) = 2 × area(T)` (or the cross-product norm for a
    /// spatial triangle, `D = 3`).
    fn jacobe(&self, parametric_coordinates: [f64; 2]) -> f64 {
        self.as_tri3().jacobe(parametric_coordinates)
    }

    /// `J[0][i] = ∂xᵢ/∂ξ`, `J[1][i] = ∂xᵢ/∂η` (constant for a straight triangle).
    fn jacobe_mat(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 2] {
        self.as_tri3().jacobe_mat(parametric_coordinates)
    }

    /// Inverse isoparametric map: parametric coordinates of `x` inside the
    /// triangle, or `None` if it lies outside (see [`ElementBase::local_coordinates`]).
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 2]> {
        self.as_tri3().local_coordinates(x)
    }
}

impl<const D: usize> Element<D, 2, 3> for CRTri3<D> {
    /// CR basis functions `θᵢ = 1 − 2λᵢ`, with `λ = (1 − ξ − η, ξ, η)` the
    /// barycentric coordinates of the reference triangle.
    fn shape(&self, parametric_coordinates: [f64; 2]) -> [f64; 3] {
        let ξ = parametric_coordinates[0];
        let η = parametric_coordinates[1];
        let lambda = [1.0 - ξ - η, ξ, η];
        std::array::from_fn(|i| 1.0 - 2.0 * lambda[i])
    }

    /// Physical gradient `∇ₓθᵢ = −2 ∇ₓλᵢ`; constant on the element.
    fn derivative_shape(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 3] {
        let grad_lambda = self.as_tri3().derivative_shape(parametric_coordinates);
        std::array::from_fn(|i| std::array::from_fn(|k| -2.0 * grad_lambda[i][k]))
    }
}

/// Trace of a [`CRTri3`] on one of its three facets: the CR boundary segment.
///
/// A nonconforming CR triangle has no vertex degrees of freedom, and the field on an
/// edge is a linear function of **all three** edge dofs of the adjacent triangle
/// (only one of them belongs to that edge). A boundary integral — a Neumann flux, a
/// Robin term, a Dirichlet penalty or multiplier — therefore cannot use [`Seg2`],
/// whose dofs are vertex values: it needs this element, whose dofs are the parent's
/// three edge dofs and whose shape functions are the restrictions `θₐ|edge`.
///
/// # Local facet numbering
///
/// `facet = k` means "the facet opposite parent vertex `k`" (0-based), exactly the
/// ordering of [`CRTri3`]'s dofs (`nodes[i]` is the edge opposite vertex `i`).
/// Nothing else is numbered: [`ElementBase::id`] returns the parent's three edge dofs
/// unchanged, so a trace element and its parent share one global numbering. The
/// facet's own parametric direction is fixed by listing its vertices in **ascending
/// parent index** order, which makes `shape` independent of how a mesh reader
/// ordered the boundary element.
///
/// # Geometry
///
/// The geometry is the facet's straight segment, delegated to [`Seg2`]
/// (`jacobe = L/2`); the unknown basis is the CR *trace*, so the element is **not**
/// isoparametric: `Σₐ shapeₐ xₐ` is not the coordinate map
/// (`vertices_coordinates` returns the parent vertices). This is harmless, boundary
/// operators only use `shape`, `jacobe` and `coordinates`.
///
/// See `book/src/theory-cr.md` for the derivation.
#[derive(Clone, Debug, PartialEq)]
pub struct CRSeg3<const D: usize> {
    facets: [usize; 3],
    vertices: [Node<D>; 3],
    facet: usize,
}

impl<const D: usize> CRSeg3<D> {
    /// Builds a facet element from the parent's three edge dofs, the parent's three
    /// vertices (local order) and the local facet index (`0..3`).
    ///
    /// Returns `None` if `facet` is not in `0..3`.
    pub fn new(facets: [usize; 3], vertices: [Node<D>; 3], facet: usize) -> Option<Self> {
        (facet < 3).then_some(CRSeg3 {
            facets,
            vertices,
            facet,
        })
    }

    /// Traces the CR triangle `parent` on its local facet `facet`.
    ///
    /// Returns `None` if `facet` is not in `0..3`.
    pub fn from_crtri3(parent: &CRTri3<D>, facet: usize) -> Option<Self> {
        let [a, b, c] = parent.triangle_vertices();
        Self::new(
            parent.id(),
            [Node::new(0, a), Node::new(1, b), Node::new(2, c)],
            facet,
        )
    }

    /// The local facet index: the facet opposite parent vertex `local_facet()`.
    pub fn local_facet(&self) -> usize {
        self.facet
    }

    /// The parent element's vertices, in local order.
    pub fn parent_vertices(&self) -> &[Node<D>; 3] {
        &self.vertices
    }

    /// The facet's two vertices, in ascending parent index order.
    pub fn facet_vertices(&self) -> [Node<D>; 2] {
        let mut out = [self.vertices[0]; 2];
        let mut next = 0;
        for (index, vertex) in self.vertices.iter().enumerate() {
            if index != self.facet {
                out[next] = *vertex;
                next += 1;
            }
        }
        out
    }

    /// The facet geometry as a linear segment.
    fn as_seg2(&self) -> Seg2<D> {
        let [start, end] = self.facet_vertices();
        Seg2::new(start, end)
    }
}

impl<const D: usize> ElementBase<D, 1, 3> for CRSeg3<D> {
    fn id(&self) -> [usize; 3] {
        self.facets
    }

    /// The parent's vertices (the element's nodes are facets, not points).
    fn vertices_coordinates(&self) -> [[f64; D]; 3] {
        self.vertices.map(|vertex| vertex.coordinates)
    }

    fn coordinates(&self, parametric_coordinates: [f64; 1]) -> [f64; D] {
        self.as_seg2().coordinates(parametric_coordinates)
    }

    /// Length scaling of the facet: `L/2`, so that `∫_facet ds = L` with the
    /// Gauss–Legendre rules of the reference interval `[-1, 1]`.
    fn jacobe(&self, parametric_coordinates: [f64; 1]) -> f64 {
        self.as_seg2().jacobe(parametric_coordinates)
    }

    fn jacobe_mat(&self, parametric_coordinates: [f64; 1]) -> [[f64; D]; 1] {
        self.as_seg2().jacobe_mat(parametric_coordinates)
    }

    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 1]> {
        self.as_seg2().local_coordinates(x)
    }
}

impl<const D: usize> Element<D, 1, 3> for CRSeg3<D> {
    /// Trace basis on the facet: with `λ` the parent's barycentric coordinates
    /// (zero at the opposite vertex, the facet's own coordinates elsewhere),
    /// `shapeₐ = 1 − 2λₐ`. With the facet's vertices ordered by ascending parent
    /// index `lo < hi`, this is `shape[facet] = 1`, `shape[lo] = ξ`,
    /// `shape[hi] = −ξ`, i.e. `facet = 0 → [1, ξ, −ξ]`,
    /// `facet = 1 → [ξ, 1, −ξ]`, `facet = 2 → [ξ, −ξ, 1]`.
    fn shape(&self, parametric_coordinates: [f64; 1]) -> [f64; 3] {
        let ξ = parametric_coordinates[0];
        let facet_lambda = [(1.0 - ξ) / 2.0, (1.0 + ξ) / 2.0];

        let mut lambda = [0.0_f64; 3];
        let mut next = 0;
        for (index, value) in lambda.iter_mut().enumerate() {
            if index != self.facet {
                *value = facet_lambda[next];
                next += 1;
            }
        }

        std::array::from_fn(|a| 1.0 - 2.0 * lambda[a])
    }

    /// Tangential gradient of the trace along the facet, `∓2 ∇Nᵢ` of the facet's
    /// linear basis (zero for the opposite vertex). Boundary operators do not use it.
    fn derivative_shape(&self, parametric_coordinates: [f64; 1]) -> [[f64; D]; 3] {
        let grad = self.as_seg2().derivative_shape(parametric_coordinates);
        let mut out = [[0.0_f64; D]; 3];
        let mut next = 0;
        for (index, row) in out.iter_mut().enumerate() {
            if index != self.facet {
                *row = std::array::from_fn(|d| -2.0 * grad[next][d]);
                next += 1;
            }
        }
        out
    }
}

impl<const D: usize> IntegrationScheme<1, 1> for CRSeg3<D> {
    fn points() -> [[f64; 1]; 1] {
        <GaussLegendre1D<1> as IntegrationScheme<1, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <GaussLegendre1D<1> as IntegrationScheme<1, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<1, 2> for CRSeg3<D> {
    fn points() -> [[f64; 1]; 2] {
        <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::points()
    }
    fn weights() -> [f64; 2] {
        <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::weights()
    }
}

impl<const D: usize> FullIntegration<1, 2> for CRSeg3<D> {}
impl<const D: usize> ReducedIntegration<1, 1> for CRSeg3<D> {}

#[cfg(test)]
mod tests {
    use super::{CRSeg3, CRTri3};
    use crate::element::tri3::Tri3;
    use crate::element::{Element, ElementBase, IntegrationScheme};
    use crate::node::Node;
    use std::collections::HashMap;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    /// Right triangle (0,0), (3,0), (0,4) as a `CRTri3`: vertices 1, 2, 3 and edge
    /// midpoints M₁ = (1.5, 2), M₂ = (0, 2), M₃ = (1.5, 0); area 6.
    fn right_triangle_cr() -> CRTri3<2> {
        CRTri3::new(
            Node::new(0, [1.5, 2.0]),
            Node::new(1, [0.0, 2.0]),
            Node::new(2, [1.5, 0.0]),
        )
    }

    #[test]
    fn crtri3_shape_is_delta_on_edge_midpoints_and_partitions_unity() {
        let cr = right_triangle_cr();

        // The midpoints in parametric coordinates: M₁ = (1/2, 1/2), M₂ = (0, 1/2),
        // M₃ = (1/2, 0).
        assert_eq!(cr.shape([0.5, 0.5]), [1.0, 0.0, 0.0]);
        assert_eq!(cr.shape([0.0, 0.5]), [0.0, 1.0, 0.0]);
        assert_eq!(cr.shape([0.5, 0.0]), [0.0, 0.0, 1.0]);

        // Not interpolatory: at the vertices θᵢ takes −1 at its own vertex and +1
        // at the other two.
        assert_eq!(cr.shape([0.0, 0.0]), [-1.0, 1.0, 1.0]);
        assert_eq!(cr.shape([1.0, 0.0]), [1.0, -1.0, 1.0]);
        assert_eq!(cr.shape([0.0, 1.0]), [1.0, 1.0, -1.0]);

        // Partition of unity everywhere → constant fields are reproduced.
        for xi in [[0.0, 0.0], [0.1, 0.7], [0.4, 0.4], [1.0 / 3.0, 1.0 / 3.0]] {
            let sum: f64 = cr.shape(xi).iter().sum();
            assert!(close(sum, 1.0));
        }
    }

    #[test]
    fn crtri3_geometry_matches_the_underlying_triangle() {
        let cr = right_triangle_cr();
        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );

        // Midpoints are recovered from the vertices by reflection.
        assert_eq!(cr.triangle_vertices(), tri.vertices_coordinates());
        assert_eq!(cr.id(), [0, 1, 2]);

        // jacobe = 2 × area = 12, and the isoparametric map is the affine map.
        assert!(close(cr.jacobe([0.0, 0.0]), 12.0));
        for xi in [[0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0], [0.6, 0.3], [0.2, 0.7]] {
            let a = cr.coordinates(xi);
            let b = tri.coordinates(xi);
            assert!(close(a[0], b[0]) && close(a[1], b[1]));
        }

        // Physical gradients are −2 ∇λᵢ: (2/3, 1/2), (−2/3, 0), (0, −1/2), and they
        // sum to zero (constant fields are in the kernel).
        let g = cr.derivative_shape([0.25, 0.25]);
        assert!(close(g[0][0], 2.0 / 3.0) && close(g[0][1], 0.5));
        assert!(close(g[1][0], -2.0 / 3.0) && close(g[1][1], 0.0));
        assert!(close(g[2][0], 0.0) && close(g[2][1], -0.5));
        for i in 0..2 {
            let sum: f64 = g.iter().map(|ga| ga[i]).sum();
            assert!(close(sum, 0.0));
        }
    }

    #[test]
    fn crtri3_edge_averages_are_kronecker_delta() {
        // The defining property: (1/|eⱼ|) ∫_{eⱼ} θᵢ ds = δᵢⱼ. Integrate along
        // every edge with the 2-point Gauss rule (exact: θᵢ is linear on the edge).
        let cr = right_triangle_cr();
        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        let vertices = tri.vertices_coordinates();
        let edges = [[1_usize, 2], [0, 2], [0, 1]];
        let t = [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()];

        for (j, [a, b]) in edges.into_iter().enumerate() {
            let length: f64 = (0..2)
                .map(|k| (vertices[b][k] - vertices[a][k]).powi(2))
                .sum::<f64>()
                .sqrt();
            let mut integral = [0.0_f64; 3];
            for parameter in t {
                let point: [f64; 2] = std::array::from_fn(|k| {
                    vertices[a][k] + parameter * (vertices[b][k] - vertices[a][k])
                });
                let xi = cr.local_coordinates(point).unwrap();
                let shape = cr.shape(xi);
                for (i, value) in integral.iter_mut().enumerate() {
                    *value += shape[i] * length / 2.0;
                }
            }
            for (i, value) in integral.iter().enumerate() {
                let expected = if i == j { length } else { 0.0 };
                assert!(
                    (value - expected).abs() < 1e-12,
                    "∫_{{e{j}}} θ{i} = {value} != {expected}"
                );
            }
        }
    }

    #[test]
    fn crtri3_mass_matrix_is_diagonal_and_needs_three_points() {
        let cr = right_triangle_cr();
        let area = 6.0;

        // Three-point rule: the quadratic integrand is integrated exactly and the
        // 2D CR mass matrix is (|T|/3) I.
        let xi = <CRTri3<2> as IntegrationScheme<2, 3>>::points();
        let w = <CRTri3<2> as IntegrationScheme<2, 3>>::weights();
        let jacobe = cr.jacobe([0.0, 0.0]);
        let mut m = [[0.0_f64; 3]; 3];
        for g in 0..3 {
            let shape = cr.shape(xi[g]);
            for i in 0..3 {
                for j in 0..3 {
                    m[i][j] += shape[i] * shape[j] * jacobe * w[g];
                }
            }
        }
        for (i, row) in m.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let expected = if i == j { area / 3.0 } else { 0.0 };
                assert!(
                    (value - expected).abs() < 1e-10,
                    "M[{i}][{j}] = {value} != {expected}"
                );
            }
        }

        // One-point rule (centroid): θᵢ = 1/3, so all entries become |T|/9 — the
        // mass matrix must therefore not be integrated with the reduced rule.
        let one_point = <CRTri3<2> as IntegrationScheme<2, 1>>::points()[0];
        let one_weight = <CRTri3<2> as IntegrationScheme<2, 1>>::weights()[0];
        let shape = cr.shape(one_point);
        for value in shape {
            assert!(close(value, 1.0 / 3.0));
        }
        let m1 = |i: usize, j: usize| shape[i] * shape[j] * jacobe * one_weight;
        assert!(close(m1(0, 0), area / 9.0));
        assert!(close(m1(0, 1), area / 9.0));
    }

    #[test]
    fn crtri3_stiffness_is_the_scaled_barycentric_form() {
        // Kᵢⱼ = 4 |T| ∇λᵢ·∇λⱼ, and the constant integrand makes the one-point rule
        // exact. For ∇λ₁ = (−1/3, −1/4), ∇λ₂ = (1/3, 0), ∇λ₃ = (0, 1/4):
        // K₁₁ = 25/6, K₁₂ = −8/3, K₁₃ = −3/2, K₂₂ = 8/3, K₂₃ = 0, K₃₃ = 3/2.
        let cr = right_triangle_cr();
        let jacobe = cr.jacobe([0.0, 0.0]);
        let one_point = <CRTri3<2> as IntegrationScheme<2, 1>>::points()[0];
        let one_weight = <CRTri3<2> as IntegrationScheme<2, 1>>::weights()[0];
        let g = cr.derivative_shape(one_point);
        let mut k = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..2).map(|d| g[i][d] * g[j][d]).sum();
                k[i][j] = dot * jacobe * one_weight;
            }
        }

        assert!(close(k[0][0], 25.0 / 6.0));
        assert!(close(k[0][1], -8.0 / 3.0));
        assert!(close(k[0][2], -1.5));
        assert!(close(k[1][1], 8.0 / 3.0));
        assert!(close(k[1][2], 0.0));
        assert!(close(k[2][2], 1.5));
        for row in &k {
            assert!(close(row.iter().sum::<f64>(), 0.0));
        }
    }

    #[test]
    fn crseg3_trace_basis() {
        let cr = right_triangle_cr();

        // The three facets, indexed by the opposite parent vertex: explicit table.
        // With the facet vertices in ascending parent-index order `lo < hi`, the
        // trace is `shape[facet] = 1`, `shape[lo] = ξ`, `shape[hi] = −ξ`.
        let expected = [
            // facet 0 (edge opposite vertex 1, vertices (1, 2)): [1, ξ, −ξ]
            |ξ: f64| [1.0, ξ, -ξ],
            // facet 1 (vertices (0, 2)): [ξ, 1, −ξ]
            |ξ: f64| [ξ, 1.0, -ξ],
            // facet 2 (vertices (0, 1)): [ξ, −ξ, 1]
            |ξ: f64| [ξ, -ξ, 1.0],
        ];

        for (facet, formula) in expected.iter().enumerate() {
            let trace = CRSeg3::from_crtri3(&cr, facet).unwrap();
            assert_eq!(trace.local_facet(), facet);
            assert_eq!(trace.id(), cr.id(), "the dofs are the parent's facets");
            for ξ in [-1.0, -0.4, 0.0, 0.3, 1.0] {
                let shape = trace.shape([ξ]);
                let expected_shape = formula(ξ);
                for a in 0..3 {
                    assert!(
                        (shape[a] - expected_shape[a]).abs() < 1e-12,
                        "facet {facet}, ξ = {ξ}, basis {a}: {} != {}",
                        shape[a],
                        expected_shape[a]
                    );
                }
                assert!(close(shape.iter().sum::<f64>(), 1.0));
                // The facet's own dof is not interpolatory either: it equals 1 on the
                // whole facet (that is its edge average), while the opposite vertex
                // vanishes.
                assert!(close(shape[facet], 1.0));
            }
        }

        // Out-of-range facet indices are rejected.
        assert!(CRSeg3::from_crtri3(&cr, 3).is_none());
    }

    #[test]
    fn crseg3_matches_the_parent_trace() {
        // The trace element's basis must equal the parent's CR basis evaluated at the
        // corresponding parent barycentric coordinates.
        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        let cr = right_triangle_cr();

        for facet in 0..3 {
            let trace = CRSeg3::from_crtri3(&cr, facet).unwrap();
            for ξ in [-1.0, -0.5, 0.0, 0.25, 1.0] {
                let point = trace.coordinates([ξ]);
                let xi = tri.local_coordinates(point).unwrap();
                let parent = cr.shape(xi);
                let shape = trace.shape([ξ]);
                for a in 0..3 {
                    assert!(
                        (shape[a] - parent[a]).abs() < 1e-12,
                        "facet {facet}, ξ = {ξ}, basis {a}: {} != {}",
                        shape[a],
                        parent[a]
                    );
                }
            }
        }
    }

    #[test]
    fn crseg3_geometry_is_the_facet() {
        // Right triangle (0,0), (3,0), (0,4): facet lengths 5 (opposite vertex 1),
        // 4 (opposite vertex 2), 3 (opposite vertex 3); jacobe = L / 2.
        let cr = right_triangle_cr();
        let lengths = [5.0, 4.0, 3.0];
        let vertices = [
            [[3.0, 0.0], [0.0, 4.0]],
            [[0.0, 0.0], [0.0, 4.0]],
            [[0.0, 0.0], [3.0, 0.0]],
        ];

        for facet in 0..3 {
            let trace = CRSeg3::from_crtri3(&cr, facet).unwrap();
            assert!(close(trace.jacobe([0.0]), lengths[facet] / 2.0));
            let ends = trace.facet_vertices();
            assert_eq!(ends[0].coordinates, vertices[facet][0]);
            assert_eq!(ends[1].coordinates, vertices[facet][1]);
            // Midpoint of the facet at ξ = 0.
            let mid = trace.coordinates([0.0]);
            for (d, value) in mid.iter().enumerate() {
                assert!(close(
                    *value,
                    0.5 * (ends[0].coordinates[d] + ends[1].coordinates[d])
                ));
            }
            // The parent's vertices are kept in local order.
            assert_eq!(trace.parent_vertices()[1].coordinates, [3.0, 0.0]);
        }
    }

    #[test]
    fn crseg3_penalty_assembles_on_the_parent_dofs() {
        // temperature_penalty on a trace element yields α ∫ θₐ θ_b ds over the
        // parent's facet dofs. The row sums are α ∫ θₐ ds = α |facet| δ_{a,facet}.
        use crate::operation::heat::temperature_penalty;

        let cr = right_triangle_cr();
        let facet = CRSeg3::from_crtri3(&cr, 0).unwrap(); // edge (3,0)-(0,4), length 5
        let alpha = 10.0;
        let (triplets, doublets) = temperature_penalty(&vec![facet], alpha, |_| 1.0);

        // The triplets hold the lower triangle; accumulate full row sums.
        let mut row_sum = [0.0_f64; 3];
        for t in &triplets {
            row_sum[t.row] += t.val;
            if t.row != t.col {
                row_sum[t.col] += t.val;
            }
        }
        for (a, sum) in row_sum.iter().enumerate() {
            let expected = if a == 0 { alpha * 5.0 } else { 0.0 };
            assert!(
                (sum - expected).abs() < 1e-10,
                "row {a}: {sum} != {expected}"
            );
        }

        // The load is α ∫ θₐ ū ds = α |facet| δ_{a,facet} for constant ū = 1.
        let mut load = [0.0_f64; 3];
        for d in &doublets {
            load[d.index] += d.value;
        }
        assert!(close(load[0], alpha * 5.0));
        assert!(close(load[1], 0.0));
        assert!(close(load[2], 0.0));
    }

    #[test]
    fn crtri3_from_tri3_assigns_edge_ids_and_midpoints() {
        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );

        // Edge table: ascending vertex pair → global edge dof, numbered after the
        // three vertices.
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        let mut next = 3;
        let mut edge_id = |a: usize, b: usize| {
            *edges.entry((a, b)).or_insert_with(|| {
                let id = next;
                next += 1;
                id
            })
        };

        let cr = CRTri3::from_tri3(&tri, &mut edge_id);
        // Edge i is the edge opposite vertex i: (1,2), (0,2), (0,1).
        assert_eq!(cr.id(), [3, 4, 5]);
        assert_eq!(
            cr.nodes().map(|node| node.coordinates),
            [[1.5, 2.0], [0.0, 2.0], [1.5, 0.0]]
        );
        assert_eq!(cr.triangle_vertices(), tri.vertices_coordinates());
        assert!(close(cr.jacobe([0.0, 0.0]), 12.0));
    }

    #[test]
    fn crtri3_assembles_with_the_heat_conduction_operator() {
        // The CR element is a drop-in for the existing operators: they only ask for
        // `Variable<D, 1, P, N>` (here `P = 2`, `N = 3`). Since ∇θᵢ = −2∇λᵢ, the CR
        // conduction matrix is exactly 4× the `Tri3` one (the integrand is constant,
        // so the 3-point rule reproduces it exactly for both elements).
        use crate::operation::heat::conduction_stiffness;

        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        );
        let tris = vec![tri];
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        let mut edge_id = |a: usize, b: usize| -> usize {
            let next = edges.len();
            *edges.entry((a, b)).or_insert(next + 3)
        };
        let crs: Vec<CRTri3<2>> = tris
            .iter()
            .map(|t| CRTri3::from_tri3(t, &mut edge_id))
            .collect();
        assert_eq!(crs[0].id(), [3, 4, 5]);

        let k_cr = conduction_stiffness(&crs, 1.0);
        let k_p1 = conduction_stiffness(&tris, 1.0);
        assert_eq!(k_cr.len(), k_p1.len());

        let mut cr = [[0.0_f64; 3]; 3];
        for t in &k_cr {
            cr[t.row - 3][t.col - 3] += t.val;
        }
        let mut p1 = [[0.0_f64; 3]; 3];
        for t in &k_p1 {
            p1[t.row][t.col] += t.val;
        }
        for (i, row) in cr.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                assert!(
                    (value - 4.0 * p1[i][j]).abs() < 1e-10,
                    "K_cr[{i}][{j}] = {value} != 4 * {}",
                    p1[i][j]
                );
            }
        }
    }

    #[test]
    fn crtri3_contains_and_local_coordinates() {
        let cr = right_triangle_cr();
        assert!(cr.contains([1.0, 4.0 / 3.0]));
        let lc = cr.local_coordinates([1.0, 4.0 / 3.0]).unwrap();
        assert!(close(lc[0], 1.0 / 3.0) && close(lc[1], 1.0 / 3.0));
        assert!(!cr.contains([3.0, 4.0]));

        // Spatial (D = 3) triangle: area 3 → jacobe 6; off-plane points are rejected.
        let spatial = CRTri3::new(
            Node::new(0, [1.0, 1.5, 0.0]),
            Node::new(1, [0.0, 1.5, 0.0]),
            Node::new(2, [1.0, 0.0, 0.0]),
        );
        assert!(close(spatial.jacobe([0.0, 0.0]), 6.0));
        assert!(spatial.contains([1.0, 1.0, 0.0]));
        assert!(!spatial.contains([1.0, 1.0, 0.1]));
    }
}
