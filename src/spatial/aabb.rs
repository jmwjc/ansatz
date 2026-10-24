//! Axis-aligned bounding boxes.

/// Axis-aligned bounding box in `D` spatial dimensions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb<const D: usize> {
    pub min: [f64; D],
    pub max: [f64; D],
}

impl<const D: usize> Aabb<D> {
    /// An empty box: `min = +∞`, `max = -∞`.
    pub fn empty() -> Self {
        Self {
            min: [f64::INFINITY; D],
            max: [f64::NEG_INFINITY; D],
        }
    }

    /// Whether the box contains no points.
    pub fn is_empty(&self) -> bool {
        (0..D).any(|i| self.min[i] > self.max[i])
    }

    /// Bounding box of a point cloud, or `None` if the cloud is empty.
    pub fn from_points(points: &[[f64; D]]) -> Option<Self> {
        let first = *points.first()?;
        let mut aabb = Self {
            min: first,
            max: first,
        };
        for point in &points[1..] {
            for ((min, max), &coordinate) in aabb.min.iter_mut().zip(aabb.max.iter_mut()).zip(point)
            {
                *min = min.min(coordinate);
                *max = max.max(coordinate);
            }
        }
        Some(aabb)
    }

    /// Box centered at `center` with the given half-width `radius` in every axis.
    pub fn from_center_radius(center: [f64; D], radius: f64) -> Self {
        Self {
            min: std::array::from_fn(|i| center[i] - radius),
            max: std::array::from_fn(|i| center[i] + radius),
        }
    }

    /// Whether `x` lies inside the closed box.
    pub fn contains(&self, x: [f64; D]) -> bool {
        (0..D).all(|i| x[i] >= self.min[i] && x[i] <= self.max[i])
    }

    /// Whether two closed boxes intersect.
    pub fn intersects(&self, other: &Self) -> bool {
        (0..D).all(|i| self.min[i] <= other.max[i] && other.min[i] <= self.max[i])
    }

    /// Smallest box containing both inputs.
    pub fn union(&self, other: &Self) -> Self {
        Self {
            min: std::array::from_fn(|i| self.min[i].min(other.min[i])),
            max: std::array::from_fn(|i| self.max[i].max(other.max[i])),
        }
    }

    /// Box expanded by `margin` in every direction.
    pub fn expanded(&self, margin: f64) -> Self {
        Self {
            min: std::array::from_fn(|i| self.min[i] - margin),
            max: std::array::from_fn(|i| self.max[i] + margin),
        }
    }

    /// Center point of the box.
    pub fn center(&self) -> [f64; D] {
        std::array::from_fn(|i| 0.5 * (self.min[i] + self.max[i]))
    }

    /// Side lengths of the box.
    pub fn size(&self) -> [f64; D] {
        std::array::from_fn(|i| self.max[i] - self.min[i])
    }

    /// Squared Euclidean distance from `x` to the box (zero inside).
    pub fn distance_squared(&self, x: [f64; D]) -> f64 {
        let mut sum = 0.0;
        for (i, &coordinate) in x.iter().enumerate() {
            let d = if coordinate < self.min[i] {
                self.min[i] - coordinate
            } else if coordinate > self.max[i] {
                coordinate - self.max[i]
            } else {
                0.0
            };
            sum += d * d;
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contains_and_intersects() {
        let a = Aabb {
            min: [0.0, 0.0],
            max: [1.0, 2.0],
        };
        let b = Aabb {
            min: [0.5, -1.0],
            max: [2.0, 0.5],
        };
        assert!(a.contains([0.5, 1.0]));
        assert!(!a.contains([1.5, 1.0]));
        assert!(a.intersects(&b));
        assert_eq!(
            a.union(&b),
            Aabb {
                min: [0.0, -1.0],
                max: [2.0, 2.0]
            }
        );
    }

    #[test]
    fn distance_and_expand() {
        let a = Aabb {
            min: [0.0, 0.0],
            max: [1.0, 1.0],
        };
        assert_eq!(a.distance_squared([0.5, 0.5]), 0.0);
        assert_eq!(a.distance_squared([2.0, 2.0]), 2.0);
        assert_eq!(a.expanded(0.5).min, [-0.5, -0.5]);
    }

    #[test]
    fn empty_box() {
        let a = Aabb::<2>::empty();
        assert!(a.is_empty());
        let b = Aabb {
            min: [0.0, 0.0],
            max: [1.0, 1.0],
        };
        assert_eq!(a.union(&b), b);
    }
}
