use crate::element::{Element, FullIntegration, GaussLegendre1D, IntegrationScheme, ReducedIntegration};
use crate::node::Node;

pub struct Seg2<const D: usize> {
    nodes: [Node<D>; 2],
}

impl<const D: usize> Seg2<D> {
    pub fn new(n1: Node<D>, n2: Node<D>) -> Self {
        Seg2 { nodes: [n1, n2] }
    }
}

impl<const D: usize> IntegrationScheme<1, 1> for Seg2<D> {
    fn points() -> [[f64; 1]; 1] {
        <GaussLegendre1D<1> as IntegrationScheme<1, 1>>::points()
    }
    fn weights() -> [f64; 1] {
        <GaussLegendre1D<1> as IntegrationScheme<1, 1>>::weights()
    }
}

impl<const D: usize> IntegrationScheme<1, 2> for Seg2<D> {
    fn points() -> [[f64; 1]; 2] {
        <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::points()
    }
    fn weights() -> [f64; 2] {
        <GaussLegendre1D<2> as IntegrationScheme<1, 2>>::weights()
    }
}

impl<const D: usize> ReducedIntegration<1, 1> for Seg2<D> {}
impl<const D: usize> FullIntegration<1, 2> for Seg2<D> {}

impl<const D: usize> Element<D, 1, 2> for Seg2<D> {
    fn id(&self) -> [usize; 2] {
        [self.nodes[0].id, self.nodes[1].id]
    }
    fn shape(&self, parametric_coordinates: [f64; 1]) -> [f64; 2] {
        let ξ = parametric_coordinates[0];
        [0.5 * (1.0 - ξ), 0.5 * (1.0 + ξ)]
    }
    /// Physical gradient of the linear shape functions: `∇ₓN₁ = -(x₂-x₁)/L²`,
    /// `∇ₓN₂ = (x₂-x₁)/L²`, i.e. `±t̂/L` where `t̂` is the unit tangent.
    fn derivative_shape(&self, parametric_coordinates: [f64; 1]) -> [[f64; D]; 2] {
        let l = 2.0 * self.jacobe(parametric_coordinates);
        let l2 = l * l;
        std::array::from_fn(|node| {
            std::array::from_fn(|i| {
                let t = (self.nodes[1].coordinates[i] - self.nodes[0].coordinates[i]) / l2;
                if node == 0 { -t } else { t }
            })
        })
    }
    fn vertices_coordinates(&self) -> [[f64; D]; 2] {
        [self.nodes[0].coordinates, self.nodes[1].coordinates]
    }
    fn coordinates(&self, parametric_coordinates: [f64; 1]) -> [f64; D] {
        let vertics_coordinates = self.vertices_coordinates();
        let shape = self.shape(parametric_coordinates);
        std::array::from_fn(|i| {
            vertics_coordinates[0][i] * shape[0] + vertics_coordinates[1][i] * shape[1]
        })
    }
    fn jacobe(&self, _parametric_coordinates: [f64; 1]) -> f64 {
        (0..D)
            .map(|i| (self.nodes[0].coordinates[i] - self.nodes[1].coordinates[i]).powi(2))
            .sum::<f64>()
            .sqrt()
            * 0.5
    }
    fn jacobe_mat(&self, _parametric_coordinates: [f64; 1]) -> [[f64; D]; 1] {
        [std::array::from_fn(|i| {
            0.5 * (self.nodes[1].coordinates[i] - self.nodes[0].coordinates[i])
        })]
    }
}

#[cfg(test)]
mod tests {
    use super::Seg2;
    use crate::element::Element;
    use crate::node::Node;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn seg2_physical_gradient() {
        // Line element from (0,0) to (3,4), length L = 5.
        // ∇N1 = -(3,4)/25, ∇N2 = (3,4)/25.
        let s = Seg2::new(Node::new(0, [0.0, 0.0]), Node::new(1, [3.0, 4.0]));
        let g = s.derivative_shape([0.0]);
        assert!(close(g[0][0], -3.0 / 25.0));
        assert!(close(g[0][1], -4.0 / 25.0));
        assert!(close(g[1][0], 3.0 / 25.0));
        assert!(close(g[1][1], 4.0 / 25.0));
    }
}
