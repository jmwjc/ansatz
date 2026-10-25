//! Crouzeix–Raviart (CR) nonconforming elements.
//!
//! The CR family keeps the piecewise-linear local space of the conforming simplices
//! ([`crate::element::tri3::Tri3`], [`crate::element::tet4::Tet4`]) but moves the degrees of freedom from the vertices to the
//! **facet averages**: the dof attached to facet `f` is `(1/|f|) ∫_f v ds`, which for
//! a linear function equals its value at the facet's centroid. Basis functions are
//! `θₐ = 1 − D λₐ` (`D` the space dimension, `λₐ` the barycentric coordinates), so
//! `∫_{f_b} θₐ ds = |f_b| δ_ab`.
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
//! Interior elements and their boundary traces live in one module each, because a
//! trace carries the *whole* parent element (the field on a facet is a linear
//! function of all of the parent's facet dofs, not only of the facet's own dof):
//!
//! - [`tri3`]: [`CRTri3`] (2D triangle, 3 edge dofs) with its trace [`CRSeg3`]
//!   (a segment carrying the parent's three edge dofs);
//! - [`tet4`]: [`CRTet4`] (3D tetrahedron, 4 face dofs) with its trace [`CRTri4`]
//!   (a triangle carrying the parent's four face dofs).
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
//! | [`CRTet4`] | tetrahedron | 4 (its faces) | 3D interior element |
//! | [`CRTri4`] | triangle | 4 (inherited from the parent tetrahedron) | 3D boundary trace |
//!
//! In every CR element the number of dofs equals the number of facets of the parent
//! simplex, which is why the traces reuse `3` (2D) and `4` (3D).
//!
//! A dof belongs to a **facet**, not to a mesh vertex, so the global numbering needs a
//! facet table. `io::gmsh` builds one while reading ([`crate::io::gmsh::GmshOptions`])
//! and exposes the elements through `crtri3_elements`, `crtri3_boundary`,
//! `crtet4_elements` and `crtet4_boundary`.

pub mod tet4;
pub mod tri3;

pub use tet4::{CRTet4, CRTri4};
pub use tri3::{CRSeg3, CRTri3};
