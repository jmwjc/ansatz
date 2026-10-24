//! Broad-phase spatial index based on AABBs and a uniform grid.

use std::collections::HashMap;

use super::Aabb;

/// Broad-phase index over items with AABBs.
///
/// The current backend is a uniform grid. The public API is intentionally
/// backend-agnostic so it can later be replaced by a BVH.
pub struct SpatialIndex<T, const D: usize> {
    items: Vec<T>,
    boxes: Vec<Aabb<D>>,
    grid: HashMap<[i64; D], Vec<usize>>,
    origin: [f64; D],
    cell_size: [f64; D],
}

impl<T, const D: usize> SpatialIndex<T, D> {
    /// Builds an index over `items`; `aabb_of` returns each item's AABB.
    pub fn new(items: impl IntoIterator<Item = T>, aabb_of: impl Fn(&T) -> Aabb<D>) -> Self {
        let items: Vec<T> = items.into_iter().collect();
        let boxes: Vec<Aabb<D>> = items.iter().map(aabb_of).collect();

        let mut bounds = Aabb::empty();
        for aabb in &boxes {
            bounds = bounds.union(aabb);
        }

        if items.is_empty() || bounds.is_empty() || D == 0 {
            return Self {
                items,
                boxes,
                grid: HashMap::new(),
                origin: [0.0; D],
                cell_size: [1.0; D],
            };
        }

        let per_axis = (items.len().max(1) as f64).powf(1.0 / D as f64).ceil() as usize;
        let per_axis = per_axis.clamp(1, 256);

        let mut cell_size = [0.0; D];
        for (i, size) in cell_size.iter_mut().enumerate() {
            let span = bounds.max[i] - bounds.min[i];
            *size = if span > 0.0 {
                span / per_axis as f64
            } else {
                1.0
            };
        }

        let origin = bounds.min;
        let mut grid: HashMap<[i64; D], Vec<usize>> = HashMap::new();
        for (index, aabb) in boxes.iter().enumerate() {
            let min = cell_of(aabb.min, origin, cell_size);
            let max = cell_of(aabb.max, origin, cell_size);
            for_each_cell(min, max, |cell| {
                grid.entry(cell).or_default().push(index);
            });
        }

        Self {
            items,
            boxes,
            grid,
            origin,
            cell_size,
        }
    }

    /// Number of indexed items.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// Returns items whose AABB contains `x`.
    pub fn query_point(&self, x: [f64; D]) -> impl Iterator<Item = &T> + '_ {
        let cell = cell_of(x, self.origin, self.cell_size);
        let candidates: &[usize] = self.grid.get(&cell).map(Vec::as_slice).unwrap_or(&[]);

        candidates
            .iter()
            .filter_map(move |&index| self.boxes[index].contains(x).then_some(&self.items[index]))
    }

    /// Returns items whose AABB intersects `query`.
    pub fn query_aabb(&self, query: &Aabb<D>) -> impl Iterator<Item = &T> + '_ {
        let query = *query;
        let mut candidates = self.candidates_in(&query);
        candidates.sort_unstable();
        candidates.dedup();

        candidates.into_iter().filter_map(move |index| {
            self.boxes[index]
                .intersects(&query)
                .then_some(&self.items[index])
        })
    }

    /// Returns items whose AABB is within `radius` of `center`.
    pub fn query_radius(&self, center: [f64; D], radius: f64) -> impl Iterator<Item = &T> + '_ {
        let query = Aabb::from_center_radius(center, radius);
        let mut candidates = self.candidates_in(&query);
        candidates.sort_unstable();
        candidates.dedup();

        let radius_squared = radius * radius;
        candidates.into_iter().filter_map(move |index| {
            (self.boxes[index].distance_squared(center) <= radius_squared)
                .then_some(&self.items[index])
        })
    }

    fn candidates_in(&self, query: &Aabb<D>) -> Vec<usize> {
        if self.items.is_empty() || query.is_empty() || D == 0 {
            return Vec::new();
        }

        let min = cell_of(query.min, self.origin, self.cell_size);
        let max = cell_of(query.max, self.origin, self.cell_size);
        let mut candidates = Vec::new();
        for_each_cell(min, max, |cell| {
            if let Some(bucket) = self.grid.get(&cell) {
                candidates.extend_from_slice(bucket);
            }
        });
        candidates
    }
}

fn cell_of<const D: usize>(x: [f64; D], origin: [f64; D], cell_size: [f64; D]) -> [i64; D] {
    std::array::from_fn(|i| ((x[i] - origin[i]) / cell_size[i]).floor() as i64)
}

fn for_each_cell<const D: usize>(min: [i64; D], max: [i64; D], mut f: impl FnMut([i64; D])) {
    let mut cell = min;
    loop {
        f(cell);

        let mut axis = D;
        loop {
            if axis == 0 {
                return;
            }
            axis -= 1;
            if cell[axis] < max[axis] {
                cell[axis] += 1;
                break;
            }
            cell[axis] = min[axis];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point_aabb(x: [f64; 2]) -> Aabb<2> {
        Aabb { min: x, max: x }
    }

    #[test]
    fn point_query_matches_brute_force() {
        let points = vec![[0.0, 0.0], [0.5, 0.5], [1.0, 0.0], [2.0, 2.0]];
        let index = SpatialIndex::new(points.clone(), |x| point_aabb(*x));

        assert_eq!(
            index.query_point([0.5, 0.5]).copied().collect::<Vec<_>>(),
            vec![[0.5, 0.5]]
        );
        assert_eq!(
            index.query_point([3.0, 3.0]).copied().collect::<Vec<_>>(),
            Vec::<[f64; 2]>::new()
        );
    }

    #[test]
    fn aabb_and_radius_queries() {
        let points = vec![[0.0, 0.0], [1.0, 0.0], [5.0, 0.0]];
        let index = SpatialIndex::new(points, |x| point_aabb(*x));

        let query = Aabb {
            min: [-0.5, -0.5],
            max: [1.5, 0.5],
        };
        let mut found: Vec<_> = index.query_aabb(&query).copied().collect();
        found.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
        assert_eq!(found, vec![[0.0, 0.0], [1.0, 0.0]]);

        let mut found: Vec<_> = index.query_radius([0.0, 0.0], 1.1).copied().collect();
        found.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
        assert_eq!(found, vec![[0.0, 0.0], [1.0, 0.0]]);
    }

    #[test]
    fn empty_index() {
        let index: SpatialIndex<[f64; 2], 2> = SpatialIndex::new(Vec::new(), |x| point_aabb(*x));
        assert!(index.is_empty());
        assert_eq!(index.query_point([0.0, 0.0]).count(), 0);
    }
}
