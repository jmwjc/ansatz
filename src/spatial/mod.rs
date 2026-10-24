//! Spatial data structures shared by FEM point location and meshfree neighbor search.

pub mod aabb;
pub mod index;

pub use aabb::Aabb;
pub use index::SpatialIndex;
