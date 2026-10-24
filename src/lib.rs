#![cfg_attr(
    not(test),
    deny(clippy::unwrap_used, clippy::expect_used, clippy::panic)
)]

//! `ansatz` — an approximation layer that turns a continuum mechanics problem into a
//! sparse algebraic system.
//!
//! An *ansatz* is the assumed form of the solution. In a finite element method it is a
//! piecewise polynomial on a mesh; in a meshfree method it is a moving-least-squares or
//! reproducing-kernel approximation on a point cloud. This crate abstracts over that choice:
//! it owns the approximation (shape functions), the numerical integration, and the assembly
//! into sparse triplets/doublets — and nothing else. Solving the resulting system is left to
//! an external sparse backend (`faer` by default), so discretization and linear algebra stay
//! independent.
//!
//! Everything is generic over the spatial dimension `D`, the parametric dimension `P`
//! (`P < D` covers elements embedded in a higher-dimensional space), the number of nodes `N`,
//! and the tensor order `O` of the unknown field (0 = scalar, 1 = vector, …), and the physics
//! is written in index (summation) notation rather than in fixed Voigt matrices. The same code
//! path therefore serves 1D/2D/3D elasticity, heat conduction, and structural elements.
//!
//! Main modules:
//!
//! - [`node`]: node types.
//! - [`element`]: the element trait (shape functions, Jacobian, etc.) and concrete elements (`Seg2`, `Tri3`, `Poi1`).
//! - [`operation`]: assembly abstractions (`Variable`, stiffness/load/penalty operators) and post-processing (`value`, `gradient`, `divergence`, `curl`, nodal recovery).
//! - [`io`]: mesh file reading (GMSH).
//! - [`spatial`]: AABB and broad-phase spatial index shared by FEM point location and meshfree neighbour search.

pub mod element;
pub mod io;
pub mod node;
pub mod operation;
pub mod sparse_vector;
pub mod spatial;
