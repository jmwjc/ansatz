use crate::element::{Element, FullIntegration, IntegrationScheme, Point, ReducedIntegration};
use crate::node::Node;

pub struct Poi1<const D: usize> {
    pub nodes: Node<D>,
}

impl<const D: usize> Poi1<D> {
    pub fn new(node: Node<D>) -> Self {
        Poi1 { nodes: node }
    }
}

impl<const D: usize> IntegrationScheme<0, 1> for Poi1<D> {
    fn points() -> [[f64; 0]; 1] {
        <Point as IntegrationScheme<0, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <Point as IntegrationScheme<0, 1>>::weights()
    }
}

// Full and reduced integration coincide for a point element: a single point in zero-dimensional space.
impl<const D: usize> FullIntegration<0, 1> for Poi1<D> {}
impl<const D: usize> ReducedIntegration<0, 1> for Poi1<D> {}

impl<const D: usize> Element<D, 0, 1> for Poi1<D> {
    fn id(&self) -> [usize; 1] {
        [self.nodes.id]
    }

    /// The shape function of a point element is identically 1 (the zero-dimensional parametric space has no coordinates).
    fn shape(&self, _parametric_coordinates: [f64; 0]) -> [f64; 1] {
        [1.0]
    }

    /// The shape function is constant, so its physical gradient is the zero
    /// `D`-dimensional vector.
    fn derivative_shape(&self, _parametric_coordinates: [f64; 0]) -> [[f64; D]; 1] {
        [[0.0; D]; 1]
    }

    fn vertices_coordinates(&self) -> [[f64; D]; 1] {
        [self.nodes.coordinates]
    }

    fn coordinates(&self, _parametric_coordinates: [f64; 0]) -> [f64; D] {
        self.nodes.coordinates
    }

    /// A zero-dimensional point element has no stretching, so the Jacobian is identically 1.
    fn jacobe(&self, _parametric_coordinates: [f64; 0]) -> f64 {
        1.0
    }

    /// The Jacobian of a point element has no parametric directions: an empty
    /// `0 × D` matrix.
    fn jacobe_mat(&self, _parametric_coordinates: [f64; 0]) -> [[f64; D]; 0] {
        []
    }
}

#[cfg(test)]
mod tests {
    use super::Poi1;
    use crate::element::Element;
    use crate::node::Node;

    #[test]
    fn poi1_shape_and_coordinates() {
        let p = Poi1::new(Node::new(3, [1.0, 2.0]));
        assert_eq!(p.id(), [3]);
        assert_eq!(p.shape([]), [1.0]);
        assert_eq!(p.derivative_shape([]), [[0.0, 0.0]]);
        assert_eq!(p.coordinates([]), [1.0, 2.0]);
        assert_eq!(p.jacobe([]), 1.0);
    }
}
