//! Raviart–Thomas element on triangles: [`RTTri3`], the lowest-order
//! `H(div)`-conforming triangle (RT0).
//!
//! The RT0 local space is the three-dimensional space of affine vector fields
//!
//! ```text
//! RT₀(T) = [P₀(T)]² ⊕ x·P₀(T) = { q(x) = a + c x : a ∈ R², c ∈ R },
//! ```
//!
//! and its degrees of freedom are the **normal fluxes** through the edges,
//!
//! ```text
//! Φᵢ(q) = ∫_{Γᵢ} q · nᵢ ds,
//! ```
//!
//! with `Γᵢ` the edge **opposite vertex `i`** and `nᵢ` its **outward** unit
//! normal. The basis function of edge `i` is
//!
//! ```text
//! qᵢ(x) = (x − Aᵢ) / (2 |T|),
//! ```
//!
//! for which `Φⱼ(qᵢ) = δᵢⱼ` and `∇·qᵢ = 1/|T|`. Consequently a global RT0
//! space is **not** contained in `H¹` (the tangential component jumps), but its
//! normal trace is single-valued on every interior edge: the space is
//! `H(div)`-conforming, which is exactly what a mixed (flux/potential) method
//! needs. See `book/src/theory-rt.md` for the theory, and
//! [`crate::element::crtri3`] for the `H¹`-nonconforming counterpart — after
//! static condensation the hybridized RT0 element reproduces the CR stiffness.
//!
//! # Relationship with the element traits
//!
//! The element interface is split into a shared [`ElementBase`] (dof ids, node
//! coordinates, isoparametric geometry) and two specialisations:
//! [`Element`](crate::element::Element) for scalar fields, whose `shape`
//! returns one value per dof, and [`HdivElement`] for vector-valued
//! `H(div)` fields. `RTTri3` implements `ElementBase` and [`HdivElement`]: the
//! vector basis is [`HdivElement::shape`] and its divergence is
//! [`HdivElement::divergence_shape`]. It is **not** a scalar `Element`, so the
//! scalar assembly operators of [`crate::operation`] do not apply: a mixed
//! method assembles the mass matrix ([`RTTri3::mass_matrix`]) and the
//! divergence incidence instead of an element stiffness.
//!
//! # Dof numbering and geometry
//!
//! As in [`crate::element::crtri3::CRTri3`], the dofs belong to **edges**, so
//! `RTTri3` stores the three **edge midpoints** in [`Node`]s whose id is the
//! global edge-dof id; the affine geometry is recovered exactly by reflection,
//! `Aᵢ = Mⱼ + Mₖ − Mᵢ` ([`RTTri3::triangle_vertices`]), and the same edge table
//! (`io::gmsh::EdgeTable`) numbers CR and RT0 unknowns identically.
//!
//! # Basis in parametric coordinates
//!
//! The reference basis is `q̂₁ = (ξ, η)`, `q̂₂ = (ξ − 1, η)`,
//! `q̂₃ = (ξ, η − 1)`; the physical basis follows from the **contravariant
//! Piola** map `q = J q̂ / |det J|`, which preserves normal fluxes and makes
//! `Φᵢ` the natural dof. Because `|det J| = jacobe` (also for a triangle
//! embedded in 3D, where it is the cross-product norm), the map is valid for
//! `D = 2` and `D = 3` alike.

use std::convert::Infallible;

use crate::element::seg2::Seg2;
use crate::element::tri3::Tri3;
use crate::element::{
    Dunavant, ElementBase, FullIntegration, GaussLegendre1D, HdivElement, IntegrationScheme,
    ReducedIntegration,
};
use crate::node::Node;

/// Lowest-order Raviart–Thomas (RT0) triangle.
///
/// `nodes[i]` is the **midpoint of the edge opposite vertex `i`**, carrying the
/// global id of that edge's flux dof. The element has three vector-valued basis
/// functions `qᵢ(x) = (x − Aᵢ) / (2|T|)` with `Φⱼ(qᵢ) = δᵢⱼ` and
/// `∇·qᵢ = 1/|T|`; the global space is `H(div)`-conforming with **one dof per
/// edge**, the same layout as [`crate::element::crtri3::CRTri3`].
///
/// It implements [`ElementBase`] and [`HdivElement`], not the scalar
/// [`Element`](crate::element::Element): its shape functions are vectors (see
/// the module documentation).
#[derive(Clone, Debug, PartialEq)]
pub struct RTTri3<const D: usize> {
    nodes: [Node<D>; 3],
}

impl<const D: usize> RTTri3<D> {
    /// Builds the element from its three edge midpoints.
    ///
    /// `m1` is the midpoint of the edge opposite vertex 1, `m2` of the edge
    /// opposite vertex 2, `m3` of the edge opposite vertex 3. Each node's id is
    /// the global flux dof (edge) id.
    pub fn new(m1: Node<D>, m2: Node<D>, m3: Node<D>) -> Self {
        RTTri3 {
            nodes: [m1, m2, m3],
        }
    }
    /// Builds the RT0 element of the triangle `tri`.
    ///
    /// `edge_id` maps an edge, given by its two end vertex ids in ascending order,
    /// to the global id of that edge's dof. It is called three times per element,
    /// so callers normally pass a closure that looks the edge up in — or inserts it
    /// into — an edge table. The same numbering as for
    /// [`CRTri3`](crate::element::crtri3::CRTri3) can be reused.
    pub fn from_tri3(tri: &Tri3<D>, edge_id: impl FnMut(usize, usize) -> usize) -> Self {
        let mut edge_id = edge_id;
        // `Infallible` makes the fallible path total, so this cannot fail.
        let built = Self::try_from_tri3(tri, |a, b| Ok::<usize, Infallible>(edge_id(a, b)));
        match built {
            Ok(rt) => rt,
            Err(never) => match never {},
        }
    }

    /// Fallible variant of [`RTTri3::from_tri3`], for edge tables that only know
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

        Ok(RTTri3 { nodes })
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
    /// The RT0 basis `qᵢ = (x − Aᵢ)/(2|T|)` is written in terms of the vertex
    /// **opposite** edge `i`, so the vertices (not the midpoints) are the
    /// geometric data the kernel needs.
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
    /// map) and for the parametric coordinates.
    fn as_tri3(&self) -> Tri3<D> {
        let [a, b, c] = self.triangle_vertices();
        Tri3::new(Node::new(0, a), Node::new(1, b), Node::new(2, c))
    }
    /// Area of the triangle (length-scaled: `jacobe / 2`, valid for `D = 3` too).
    pub fn area(&self) -> f64 {
        0.5 * self.jacobe([0.0, 0.0])
    }
    /// Whether the physical point `x` lies inside this element.
    pub fn contains(&self, x: [f64; D]) -> bool {
        self.local_coordinates(x).is_some()
    }
    /// Reference basis `q̂₁ = (ξ, η)`, `q̂₂ = (ξ − 1, η)`, `q̂₃ = (ξ, η − 1)`.
    ///
    /// These are the RT0 basis functions of the reference triangle
    /// `(0,0), (1,0), (0,1)`; [`HdivElement::shape`] maps them to physical
    /// space by the contravariant Piola map.
    pub fn reference_vector_shape(&self, parametric_coordinates: [f64; 2]) -> [[f64; 2]; 3] {
        let [ξ, η] = parametric_coordinates;
        [[ξ, η], [ξ - 1.0, η], [ξ, η - 1.0]]
    }
    /// Divergence of the field with edge-flux coefficients `coefficients`,
    /// `∇·q_h = (1/|T|) Σᵢ cᵢ`.
    pub fn divergence(&self, coefficients: [f64; 3]) -> f64 {
        coefficients.iter().sum::<f64>() / self.area()
    }
    /// Element mass matrix `Mᵢⱼ = ∫_T qᵢ·qⱼ dΩ`, integrated with the 3-point
    /// Dunavant rule (exact: the integrand is quadratic).
    ///
    /// Unlike the CR mass matrix this one is **dense**: `H(div)` conformity is
    /// paid for by a coupled flux mass.
    pub fn mass_matrix(&self) -> [[f64; 3]; 3] {
        let points = <Self as IntegrationScheme<2, 3>>::points();
        let weights = <Self as IntegrationScheme<2, 3>>::weights();
        let mut m = [[0.0_f64; 3]; 3];
        for (point, weight) in points.into_iter().zip(weights) {
            let q = self.shape(point);
            let w = self.jacobe(point) * weight;
            for (i, row) in m.iter_mut().enumerate() {
                for (j, value) in row.iter_mut().enumerate() {
                    let dot: f64 = (0..D).map(|d| q[i][d] * q[j][d]).sum();
                    *value += w * dot;
                }
            }
        }
        m
    }
    /// Outward unit normal of edge `edge` (the edge opposite local vertex
    /// `edge`), or `None` if `edge ≥ 3` or the triangle is degenerate.
    ///
    /// Computed without any dimension-specific formula: the vector from the
    /// opposite vertex to the edge midpoint is orthogonalized against the edge
    /// tangent, which yields the outward normal for `D = 2` and for a triangle
    /// embedded in `D = 3` alike.
    pub fn edge_normal(&self, edge: usize) -> Option<[f64; D]> {
        if edge >= 3 {
            return None;
        }
        let a = self.triangle_vertices();
        let (j, k) = ((edge + 1) % 3, (edge + 2) % 3);
        let tangent: [f64; D] = std::array::from_fn(|d| a[k][d] - a[j][d]);
        let midpoint: [f64; D] = std::array::from_fn(|d| 0.5 * (a[j][d] + a[k][d]));
        // Points from the opposite vertex towards the edge, i.e. outward.
        let out: [f64; D] = std::array::from_fn(|d| midpoint[d] - a[edge][d]);
        let tt: f64 = tangent.iter().map(|v| v * v).sum();
        if tt == 0.0 {
            return None;
        }
        let to: f64 = tangent.iter().zip(&out).map(|(a, b)| a * b).sum();
        let mut normal: [f64; D] = std::array::from_fn(|d| out[d] - to / tt * tangent[d]);
        let norm = normal.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm == 0.0 {
            return None;
        }
        for value in &mut normal {
            *value /= norm;
        }
        Some(normal)
    }
    /// Flux dof `Φᵢ(q) = ∫_{Γᵢ} q·nᵢ ds` of an arbitrary vector field `q`, or
    /// `None` if `edge ≥ 3` or the triangle is degenerate.
    ///
    /// Evaluated with the 2-point Gauss–Legendre rule on the edge, which is exact
    /// for polynomials of degree ≤ 3. For the basis itself this returns
    /// `Φⱼ(qᵢ) = δᵢⱼ`; for the element-wise gradient of a linear field it returns
    /// `ℓᵢ nᵢ·∇u`, the flux dofs used to post-process a CR solution.
    pub fn flux_dof(&self, edge: usize, q: impl Fn([f64; D]) -> [f64; D]) -> Option<f64> {
        let normal = self.edge_normal(edge)?;
        let a = self.triangle_vertices();
        let (j, k) = ((edge + 1) % 3, (edge + 2) % 3);
        let length = (0..D)
            .map(|d| (a[k][d] - a[j][d]).powi(2))
            .sum::<f64>()
            .sqrt();

        let points = <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::points();
        let weights = <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::weights();
        let mut flux = 0.0;
        for (point, weight) in points.into_iter().zip(weights) {
            // Map ξ ∈ [−1, 1] to s ∈ [0, 1] along the edge.
            let s = 0.5 * (point[0] + 1.0);
            let x: [f64; D] = std::array::from_fn(|d| a[j][d] + s * (a[k][d] - a[j][d]));
            let qx = q(x);
            let normal_flux: f64 = (0..D).map(|d| qx[d] * normal[d]).sum();
            flux += 0.5 * length * weight * normal_flux;
        }
        Some(flux)
    }
}

/// Geometry and dof ids of the RT0 triangle; see [`ElementBase`].
impl<const D: usize> ElementBase<D, 2, 3> for RTTri3<D> {
    /// Global ids of the three edge dofs (the flux unknowns).
    fn id(&self) -> [usize; 3] {
        [self.nodes[0].id, self.nodes[1].id, self.nodes[2].id]
    }
    /// Coordinates of the three edge midpoints (the element's nodes).
    ///
    /// These are *not* the triangle vertices; use
    /// [`RTTri3::triangle_vertices`] for those.
    fn vertices_coordinates(&self) -> [[f64; D]; 3] {
        [
            self.nodes[0].coordinates,
            self.nodes[1].coordinates,
            self.nodes[2].coordinates,
        ]
    }
    /// Isoparametric map built from the midpoints, `x(ξ) = Σᵢ θᵢ(ξ) Mᵢ` with
    /// `θᵢ = 1 − 2λᵢ`; it is the affine map of the underlying triangle.
    fn coordinates(&self, parametric_coordinates: [f64; 2]) -> [f64; D] {
        let λ = barycentric(parametric_coordinates);
        let midpoints = self.vertices_coordinates();
        std::array::from_fn(|i| {
            let mut x = 0.0;
            for a in 0..3 {
                x += (1.0 - 2.0 * λ[a]) * midpoints[a][i];
            }
            x
        })
    }
    /// Area scaling `√det(JᵀJ) = 2 × area(T)` (the cross-product norm when the
    /// triangle is embedded in `D = 3`); constant for a straight triangle.
    fn jacobe(&self, parametric_coordinates: [f64; 2]) -> f64 {
        self.as_tri3().jacobe(parametric_coordinates)
    }
    /// `J[0][i] = ∂xᵢ/∂ξ`, `J[1][i] = ∂xᵢ/∂η` (constant for a straight triangle).
    fn jacobe_mat(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 2] {
        self.as_tri3().jacobe_mat(parametric_coordinates)
    }
    /// Inverse isoparametric map: parametric coordinates of `x` inside the
    /// triangle, or `None` if it lies outside.
    fn local_coordinates(&self, x: [f64; D]) -> Option<[f64; 2]> {
        self.as_tri3().local_coordinates(x)
    }
}

/// RT0 as an [`HdivElement`]: vector basis and divergence.
impl<const D: usize> HdivElement<D, 2, 3> for RTTri3<D> {
    /// Physical vector basis `qᵢ(x) = (x − Aᵢ) / (2|T|)`, evaluated through the
    /// contravariant Piola map `qᵢ = J q̂ᵢ / jacobe`.
    ///
    /// Returned as `[edge][spatial]`, one `D`-vector per edge dof.
    fn shape(&self, parametric_coordinates: [f64; 2]) -> [[f64; D]; 3] {
        let j = self.jacobe_mat(parametric_coordinates);
        let jacobe = self.jacobe(parametric_coordinates);
        let q_hat = self.reference_vector_shape(parametric_coordinates);
        std::array::from_fn(|i| {
            std::array::from_fn(|d| (j[0][d] * q_hat[i][0] + j[1][d] * q_hat[i][1]) / jacobe)
        })
    }
    /// Divergence of the basis functions, `∇·qᵢ = 1/|T|` (constant; the
    /// parametric coordinates are part of the [`HdivElement`] interface and are
    /// unused here).
    fn divergence_shape(&self, _parametric_coordinates: [f64; 2]) -> [f64; 3] {
        [1.0 / self.area(); 3]
    }
}

/// Full integration: 3-point Dunavant rule (algebraic degree 2), needed for the
/// quadratic mass integrand. Reduced integration: 1 point (centroid), exact for
/// the constant divergence integrand.
impl<const D: usize> IntegrationScheme<2, 1> for RTTri3<D> {
    fn points() -> [[f64; 2]; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Dunavant<1> as IntegrationScheme<2, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<2, 3> for RTTri3<D> {
    fn points() -> [[f64; 2]; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::points()
    }
    fn weights() -> [f64; 3] {
        <Dunavant<3> as IntegrationScheme<2, 3>>::weights()
    }
}

impl<const D: usize> FullIntegration<2, 3> for RTTri3<D> {}
impl<const D: usize> ReducedIntegration<2, 1> for RTTri3<D> {}

/// Trace of an [`RTTri3`] on one of its three edges: the RT0 boundary-flux
/// element.
///
/// An RT0 flux is `H(div)`-conforming, so only its **normal** trace is
/// single-valued. On the edge opposite parent vertex `facet`, the trace
/// `qᵢ·n` of the parent's three vector basis functions is a polynomial of the
/// edge coordinate (for RT0 it is the constant `δᵢ,facet/ℓ`), and the dofs of
/// this element are the parent's three **edge-flux dofs**. Its purpose is the
/// natural Dirichlet term of the mixed heat problem,
///
/// `ℓᵢ = −∫_{Γ_D} T_D (qᵢ·n) ds`,
///
/// assembled by [`crate::operation::heat::mixed_dirichlet`].
///
/// The element is **not** isoparametric and does not implement [`Element`](crate::element::Element): its
/// geometry is the facet's straight [`Seg2`] (`jacobe = L/2`), while its
/// “shape functions” are the normal-flux traces returned by
/// [`RTSeg3::normal_flux_shape`], not scalar nodal bases. `facet = k` means
/// “the edge opposite parent vertex `k`”, so [`ElementBase::id`] returns the
/// parent's three flux dofs unchanged. The facet's parametric direction is
/// fixed by listing its vertices in ascending parent-index order.
#[derive(Clone, Debug, PartialEq)]
pub struct RTSeg3<const D: usize> {
    facets: [usize; 3],
    vertices: [Node<D>; 3],
    facet: usize,
}

impl<const D: usize> RTSeg3<D> {
    /// Builds a facet element from the parent's three flux dofs, the parent's
    /// three vertices (local order) and the local facet index (`0..3`).
    ///
    /// Returns `None` if `facet` is not in `0..3`.
    pub fn new(facets: [usize; 3], vertices: [Node<D>; 3], facet: usize) -> Option<Self> {
        (facet < 3).then_some(RTSeg3 {
            facets,
            vertices,
            facet,
        })
    }

    /// Traces the RT0 triangle `parent` on its local facet `facet`.
    ///
    /// Returns `None` if `facet` is not in `0..3`.
    pub fn from_rttri3(parent: &RTTri3<D>, facet: usize) -> Option<Self> {
        let [a, b, c] = parent.triangle_vertices();
        Self::new(
            parent.id(),
            [Node::new(0, a), Node::new(1, b), Node::new(2, c)],
            facet,
        )
    }

    /// The local facet index: the edge opposite parent vertex `local_facet()`.
    pub fn local_facet(&self) -> usize {
        self.facet
    }

    /// The parent element's vertices, in local order.
    pub fn parent_vertices(&self) -> &[Node<D>; 3] {
        &self.vertices
    }

    /// The facet's two vertices, in ascending parent index order.
    pub fn facet_vertices(&self) -> [Node<D>; 2] {
        let [lo, hi] = self.facet_vertex_indices();
        [self.vertices[lo], self.vertices[hi]]
    }

    /// Normal-flux trace of the parent's three RT0 basis functions on this
    /// facet, evaluated at the edge coordinate `ξ ∈ [−1, 1]`.
    ///
    /// Index `i` is the parent's edge dof; `n` is the parent's **outward**
    /// normal on this facet. For RT0 the result is constant and equal to
    /// `δᵢ,facet / ℓ_facet`; it is computed from the affine basis at the two
    /// facet endpoints so that it remains correct for any straight triangle.
    pub fn normal_flux_shape(&self, parametric_coordinates: [f64; 1]) -> [f64; 3] {
        let parent = self.parent();
        let Some(normal) = parent.edge_normal(self.facet) else {
            return [0.0; 3];
        };
        let vertices = parent.triangle_vertices();
        let area = parent.area();
        if area == 0.0 {
            return [0.0; 3];
        }
        let [lo, hi] = self.facet_vertex_indices();
        let value = |vertex: usize, dof: usize| -> f64 {
            (0..D)
                .map(|d| (vertices[vertex][d] - vertices[dof][d]) * normal[d])
                .sum::<f64>()
                / (2.0 * area)
        };
        let low = [value(lo, 0), value(lo, 1), value(lo, 2)];
        let high = [value(hi, 0), value(hi, 1), value(hi, 2)];
        let ξ = parametric_coordinates[0];
        std::array::from_fn(|i| 0.5 * (1.0 - ξ) * low[i] + 0.5 * (1.0 + ξ) * high[i])
    }

    /// The two parent vertex indices that span the facet, in ascending order.
    fn facet_vertex_indices(&self) -> [usize; 2] {
        let mut out = [0_usize; 2];
        let mut next = 0;
        for index in 0..3 {
            if index != self.facet {
                out[next] = index;
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

    /// The parent RT0 element reconstructed from the parent vertices and flux
    /// dofs (the midpoints are not needed for the trace).
    fn parent(&self) -> RTTri3<D> {
        let vertices = self.vertices_coordinates();
        let midpoint = |p: usize, q: usize| -> [f64; D] {
            std::array::from_fn(|d| 0.5 * (vertices[p][d] + vertices[q][d]))
        };
        RTTri3::new(
            Node::new(self.facets[0], midpoint(1, 2)),
            Node::new(self.facets[1], midpoint(0, 2)),
            Node::new(self.facets[2], midpoint(0, 1)),
        )
    }
}

impl<const D: usize> ElementBase<D, 1, 3> for RTSeg3<D> {
    fn id(&self) -> [usize; 3] {
        self.facets
    }

    /// The parent's vertices (the RT0 dofs are edge fluxes, not points).
    fn vertices_coordinates(&self) -> [[f64; D]; 3] {
        self.vertices.map(|vertex| vertex.coordinates)
    }

    fn coordinates(&self, parametric_coordinates: [f64; 1]) -> [f64; D] {
        self.as_seg2().coordinates(parametric_coordinates)
    }

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

impl<const D: usize> IntegrationScheme<1, 1> for RTSeg3<D> {
    fn points() -> [[f64; 1]; 1] {
        <GaussLegendre1D<1> as IntegrationScheme<1, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <GaussLegendre1D<1> as IntegrationScheme<1, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<1, 2> for RTSeg3<D> {
    fn points() -> [[f64; 1]; 2] {
        <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::points()
    }
    fn weights() -> [f64; 2] {
        <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::weights()
    }
}

impl<const D: usize> FullIntegration<1, 2> for RTSeg3<D> {}
impl<const D: usize> ReducedIntegration<1, 1> for RTSeg3<D> {}

/// Barycentric coordinates `λ = (1 − ξ − η, ξ, η)` of the reference triangle.
fn barycentric([ξ, η]: [f64; 2]) -> [f64; 3] {
    [1.0 - ξ - η, ξ, η]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ElementBase;
    use crate::element::HdivElement;
    use std::collections::HashMap;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    fn right_triangle() -> Tri3<2> {
        Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 0.0]),
            Node::new(2, [0.0, 4.0]),
        )
    }

    /// RT0 element of the right triangle above (area 6), with edge dofs 3, 4, 5.
    fn right_triangle_rt() -> RTTri3<2> {
        let tri = right_triangle();
        let mut next = 3;
        RTTri3::from_tri3(&tri, |_, _| {
            let id = next;
            next += 1;
            id
        })
    }

    #[test]
    fn rttri3_from_tri3_assigns_edge_ids_and_midpoints() {
        let tri = right_triangle();
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        let mut next = 3;
        let rt = RTTri3::from_tri3(&tri, |a, b| {
            *edges.entry((a, b)).or_insert_with(|| {
                let id = next;
                next += 1;
                id
            })
        });

        // Edge i is the edge opposite vertex i: (1,2), (0,2), (0,1).
        assert_eq!(rt.id(), [3, 4, 5]);
        assert_eq!(
            rt.nodes().map(|node| node.coordinates),
            [[1.5, 2.0], [0.0, 2.0], [1.5, 0.0]]
        );
        assert_eq!(rt.triangle_vertices(), tri.vertices_coordinates());
        assert!(close(rt.area(), 6.0));
        assert!(close(rt.jacobe([0.0, 0.0]), 12.0));
        // Geometry agrees with the underlying linear triangle.
        assert!(rt.contains([1.0, 4.0 / 3.0]));
        assert!(!rt.contains([-1.0, -1.0]));
    }

    #[test]
    fn rttri3_basis_forms_agree() {
        // Physical form qᵢ = (x − Aᵢ)/(2|T|), the barycentric form
        // qᵢ = [λⱼ(xⱼ − xᵢ) + λₖ(xₖ − xᵢ)]/(2|T|), and the reference + Piola form
        // must all give the same vector at every point.
        let rt = right_triangle_rt();
        let a = rt.triangle_vertices();
        let area = rt.area();

        for [ξ, η] in [[0.0, 0.0], [1.0 / 3.0, 1.0 / 3.0], [0.2, 0.5], [1.0, 0.0]] {
            let x = rt.coordinates([ξ, η]);
            let q = rt.shape([ξ, η]);
            let λ = barycentric([ξ, η]);
            for (i, a_i) in a.iter().enumerate() {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                for (d, x_d) in x.iter().enumerate() {
                    let physical = (x_d - a_i[d]) / (2.0 * area);
                    let barycentric =
                        (λ[j] * (a[j][d] - a_i[d]) + λ[k] * (a[k][d] - a_i[d])) / (2.0 * area);
                    assert!(
                        (q[i][d] - physical).abs() < 1e-13,
                        "edge {i}, d = {d}: {} != {physical}",
                        q[i][d]
                    );
                    assert!((q[i][d] - barycentric).abs() < 1e-13);
                }
            }
        }
    }

    fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
        a[0] * b[0] + a[1] * b[1]
    }

    fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    #[test]
    fn rttri3_flux_dofs_are_kronecker_deltas() {
        // Φⱼ(qᵢ) = δᵢⱼ, independently of the point at which the (affine) basis
        // is evaluated: the edge integrand is constant for i = j and affine with
        // zero mean otherwise, so the 2-point rule is exact.
        let rt = right_triangle_rt();
        let a = rt.triangle_vertices();
        let area = rt.area();
        for (i, a_i) in a.iter().enumerate() {
            // The basis function qᵢ = (x − Aᵢ)/(2|T|) as a field.
            let q_i = |x: [f64; 2]| -> [f64; 2] {
                std::array::from_fn(|d| (x[d] - a_i[d]) / (2.0 * area))
            };
            for j in 0..3 {
                let dof = rt.flux_dof(j, q_i).unwrap();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dof - expected).abs() < 1e-12,
                    "Φ_{j}(q_{i}) = {dof} != {expected}"
                );
            }
        }
        assert!(rt.flux_dof(3, |_| [0.0, 0.0]).is_none());
    }

    #[test]
    fn rttri3_divergence_is_unit_and_the_divergence_theorem_holds() {
        // ∇·qᵢ = 1/|T|, and for any coefficients the discrete divergence theorem
        // ∮ q_h·n ds = Σᵢ Φᵢ = ∫_T ∇·q_h dΩ holds exactly on the element.
        let rt = right_triangle_rt();
        let area = rt.area();
        let a = rt.triangle_vertices();
        for value in rt.divergence_shape([0.0, 0.0]) {
            assert!(close(value, 1.0 / area));
        }

        let coefficients = [1.0, -2.0, 0.5];
        let q_h = |x: [f64; 2]| -> [f64; 2] {
            std::array::from_fn(|d| {
                (0..3)
                    .map(|i| coefficients[i] * (x[d] - a[i][d]) / (2.0 * area))
                    .sum()
            })
        };
        let boundary: f64 = (0..3).map(|j| rt.flux_dof(j, q_h).unwrap()).sum::<f64>();
        let integral = area * rt.divergence(coefficients);
        assert!(close(boundary, coefficients.iter().sum::<f64>()));
        assert!(close(integral, boundary));
    }

    #[test]
    fn rttri3_mass_matrix_is_dense_and_matches_the_closed_form() {
        // ∫_T qᵢ·qⱼ dΩ = [Σ_K (x_K−xᵢ)·(x_K−xⱼ) + 9 (c−xᵢ)·(c−xⱼ)] / (48 |T|).
        let rt = right_triangle_rt();
        let a = rt.triangle_vertices();
        let area = rt.area();
        let centroid: [f64; 2] = std::array::from_fn(|d| a.iter().map(|v| v[d]).sum::<f64>() / 3.0);
        let m = rt.mass_matrix();

        for i in 0..3 {
            for j in 0..3 {
                let vertices: f64 = (0..3)
                    .map(|k| {
                        dot2(
                            std::array::from_fn(|d| a[k][d] - a[i][d]),
                            std::array::from_fn(|d| a[k][d] - a[j][d]),
                        )
                    })
                    .sum();
                let center = dot2(
                    std::array::from_fn(|d| centroid[d] - a[i][d]),
                    std::array::from_fn(|d| centroid[d] - a[j][d]),
                );
                let expected = (vertices + 9.0 * center) / (48.0 * area);
                assert!(
                    (m[i][j] - expected).abs() < 1e-12,
                    "M[{i}][{j}] = {} != {expected}",
                    m[i][j]
                );
                assert!(close(m[i][j], m[j][i]));
            }
            assert!(m[i][i] > 0.0);
        }
        // H(div) conformity couples the edges: the mass matrix is not diagonal.
        assert!(m[0][1].abs() > 1e-3);
        assert!(m[0][2].abs() > 1e-3);
        assert!(m[1][2].abs() > 1e-3);
    }

    #[test]
    fn rttri3_cr_gradient_flux_dofs_are_cr_stiffness_rows() {
        // For a linear field u the CR dofs are the midpoint values mⱼ = u(Mⱼ) and
        // the element-wise gradient is the constant vector g. Post-processing it
        // into RT0 gives Φₑ(∇u) = ℓₑ nₑ·g, which must equal Σⱼ Kᶜᴿₑⱼ mⱼ with
        // Kᶜᴿₑⱼ = ℓₑ ℓⱼ nₑ·nⱼ / |T| (the CR element stiffness, theorem 5.1 of
        // book/src/theory-rt.md read in the other direction).
        let rt = right_triangle_rt();
        let a = rt.triangle_vertices();
        let midpoints = rt.vertices_coordinates();
        let area = rt.area();

        let lengths: [f64; 3] = std::array::from_fn(|e| {
            let (j, k) = ((e + 1) % 3, (e + 2) % 3);
            dot2(
                std::array::from_fn(|d| a[k][d] - a[j][d]),
                std::array::from_fn(|d| a[k][d] - a[j][d]),
            )
            .sqrt()
        });
        let normals: [[f64; 2]; 3] = std::array::from_fn(|e| rt.edge_normal(e).unwrap());

        let gradient = [0.7, -1.3];
        // Edge averages of a linear field are its values at the edge midpoints.
        let cr_dofs: [f64; 3] = std::array::from_fn(|i| {
            2.0 + gradient[0] * midpoints[i][0] + gradient[1] * midpoints[i][1]
        });
        let stiffness =
            |e: usize, j: usize| lengths[e] * lengths[j] * dot2(normals[e], normals[j]) / area;

        for e in 0..3 {
            let expected = lengths[e] * dot2(normals[e], gradient);
            let condensed: f64 = (0..3).map(|j| stiffness(e, j) * cr_dofs[j]).sum();
            let post_processed = rt.flux_dof(e, |_| gradient).unwrap();
            assert!(close(condensed, expected));
            assert!(close(post_processed, expected));
        }
    }

    #[test]
    fn rtseg3_traces_the_parent_normal_flux() {
        // Right triangle (0,0), (3,0), (0,4): facet lengths 5 (opposite vertex 1),
        // 4 (opposite vertex 2), 3 (opposite vertex 3); the RT0 normal trace on
        // its own edge is the constant 1/ℓ, and zero on the other two edges.
        let rt = right_triangle_rt();
        let lengths = [5.0, 4.0, 3.0];
        for (facet, length) in lengths.iter().enumerate() {
            let trace = RTSeg3::from_rttri3(&rt, facet).unwrap();
            assert_eq!(trace.local_facet(), facet);
            assert_eq!(trace.id(), rt.id());
            assert!(close(trace.jacobe([0.0]), length / 2.0));

            for ξ in [-1.0, -0.25, 0.0, 0.6, 1.0] {
                let shape = trace.normal_flux_shape([ξ]);
                for (i, value) in shape.iter().enumerate() {
                    let expected = if i == facet { 1.0 / length } else { 0.0 };
                    assert!(
                        (value - expected).abs() < 1e-12,
                        "facet {facet}, ξ = {ξ}, dof {i}: {value} != {expected}"
                    );
                }
            }

            // ∫_{Γ} (q_facet·n) ds = 1 (Kronecker property of the flux dof).
            let points = <RTSeg3<2> as IntegrationScheme<1, 2>>::points();
            let weights = <RTSeg3<2> as IntegrationScheme<1, 2>>::weights();
            let integral: f64 = points
                .into_iter()
                .zip(weights)
                .map(|(xi, w)| trace.normal_flux_shape(xi)[facet] * trace.jacobe(xi) * w)
                .sum();
            assert!(close(integral, 1.0));
        }
        assert!(RTSeg3::from_rttri3(&rt, 3).is_none());
    }

    #[test]
    fn elementbase_is_shared_with_scalar_elements() {
        // The geometry half of the interface is independent of the shape value
        // type: the same generic function accepts a scalar `Element` (Tri3) and
        // the vector-valued `HdivElement` (RTTri3).
        fn geometry<T: ElementBase<2, 2, 3>>(elm: &T) -> ([f64; 2], f64, bool) {
            (
                elm.coordinates([1.0 / 3.0, 1.0 / 3.0]),
                elm.jacobe([0.0, 0.0]),
                elm.contains([0.25, 0.25]),
            )
        }
        fn vector_basis<T: HdivElement<2, 2, 3>>(elm: &T) -> [[f64; 2]; 3] {
            elm.shape([0.0, 0.0])
        }

        let tri = right_triangle();
        let rt = right_triangle_rt();
        let (tri_center, tri_jacobe, tri_inside) = geometry(&tri);
        let (rt_center, rt_jacobe, rt_inside) = geometry(&rt);
        for d in 0..2 {
            assert!(close(tri_center[d], rt_center[d]));
        }
        assert!(close(tri_jacobe, rt_jacobe));
        assert!(tri_inside && rt_inside);

        let through_trait = vector_basis(&rt);
        assert_eq!(through_trait, rt.shape([0.0, 0.0]));
    }

    #[test]
    fn rttri3_works_for_a_triangle_embedded_in_3d() {
        // Spatial triangle (0,0,0), (1,0,0), (0,1,1): area = √2/2, jacobe = √2.
        let tri = Tri3::new(
            Node::new(0, [0.0, 0.0, 0.0]),
            Node::new(1, [1.0, 0.0, 0.0]),
            Node::new(2, [0.0, 1.0, 1.0]),
        );
        let mut next = 3;
        let rt = RTTri3::from_tri3(&tri, |_, _| {
            let id = next;
            next += 1;
            id
        });
        let area = std::f64::consts::FRAC_1_SQRT_2;
        assert!(close(rt.area(), area));
        assert!(close(rt.jacobe([0.0, 0.0]), 2.0 * area));

        // The Piola basis is still qᵢ = (x − Aᵢ)/(2|T|).
        let a = rt.triangle_vertices();
        let x = rt.coordinates([0.3, 0.4]);
        let q = rt.shape([0.3, 0.4]);
        for (i, a_i) in a.iter().enumerate() {
            for (d, x_d) in x.iter().enumerate() {
                assert!((q[i][d] - (x_d - a_i[d]) / (2.0 * area)).abs() < 1e-13);
            }
        }

        // The normals are unit, orthogonal to the plane normal (i.e. they lie in
        // the triangle's plane), point outward, and the flux dofs are Kronecker.
        let t1 = [1.0, 0.0, 0.0];
        let t2 = [0.0, 1.0, 1.0];
        let plane_normal = [
            t1[1] * t2[2] - t1[2] * t2[1],
            t1[2] * t2[0] - t1[0] * t2[2],
            t1[0] * t2[1] - t1[1] * t2[0],
        ];
        for e in 0..3 {
            let n = rt.edge_normal(e).unwrap();
            assert!(close(dot3(n, n), 1.0));
            assert!(close(dot3(n, plane_normal), 0.0));
            // Outward: the opposite vertex lies on the negative side of the edge.
            let opposite = a[e];
            let midpoint: [f64; 3] =
                std::array::from_fn(|d| 0.5 * (a[(e + 1) % 3][d] + a[(e + 2) % 3][d]));
            let to_edge: [f64; 3] = std::array::from_fn(|d| midpoint[d] - opposite[d]);
            assert!(dot3(n, to_edge) > 0.0);
            for (i, a_i) in a.iter().enumerate() {
                let q_i = |x: [f64; 3]| -> [f64; 3] {
                    std::array::from_fn(|d| (x[d] - a_i[d]) / (2.0 * area))
                };
                let dof = rt.flux_dof(e, q_i).unwrap();
                let expected = if i == e { 1.0 } else { 0.0 };
                assert!((dof - expected).abs() < 1e-12);
            }
        }
    }
}
