use crate::element::{Element, FullIntegration, IntegrationScheme};
use crate::operation::{Doublets, Triplets, TripletsAndDoublets, Variable};
use crate::sparse_vector::Doublet;
use faer::sparse::Triplet;

/// Assembles the stiffness matrix and returns sparse triplets for the lower
/// triangle (including the diagonal).
///
/// `ea` is the axial stiffness (Young's modulus × cross-sectional area).
/// `v` is a set of 1D elements (implementing [`Variable<1, 1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`].
pub fn stiffness<T, const P: usize, const N: usize, const G: usize>(v: &T, ea: f64) -> Triplets
where
    T: Variable<1, 1, P, N>,
    T::Item: Element<1, P, N> + FullIntegration<P, G>,
{
    let mut triplets: Triplets = Vec::with_capacity(N * (N + 1) / 2 * G);
    let ξ = T::Item::points();
    let w = T::Item::weights();
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let dshape = elm.derivative_shape(ξg);
            let jacobe = elm.jacobe(ξg);
            for i in 0..N {
                for j in 0..i + 1 {
                    triplets.push(Triplet::new(
                        id[i].max(id[j]),
                        id[i].min(id[j]),
                        ea * dshape[i][0] * dshape[j][0] * jacobe * wg,
                    ));
                }
            }
        }
    }
    triplets
}

/// Applies a concentrated point load, using point
/// elements ([`crate::element::poi1::Poi1`]) to represent the load points.
///
/// `load` is the force vector and `v` is a set of point elements
/// ([`Variable<D, D, 0, 1>`]); each point element corresponds to one load point,
/// where `F_{i,p} += Nᵢ load_p` is applied. Returns a list of [`Doublet`]s that the
/// caller can merge into a load vector via
/// [`crate::sparse_vector::try_new_from_doublet`].
pub fn point_load<T, const D: usize>(v: &T, load: [f64; D]) -> Doublets
where
    T: Variable<D, D, 0, 1>,
    T::Item: Element<D, 0, 1>,
{
    let mut doublets = Vec::with_capacity(D);
    for elm in v.iter() {
        let id = elm.id();
        let shape = elm.shape([]);
        for (component, &value) in load.iter().enumerate() {
            doublets.push(Doublet::new(id[0] * D + component, shape[0] * value));
        }
    }
    doublets
}

/// Applies a body force (distributed load) in operator form (element-wise
/// integration).
///
/// `f` is the body force density (axial force per unit length), integrated over
/// each element:
/// `Fᵢ = ∫ Nᵢ f dx = Σ_g Nᵢ(ξ_g) f(x_g) |J| w_g`.
/// Returns a list of [`Doublet`]s that the caller can merge into a load vector
/// via [`crate::sparse_vector::try_new_from_doublet`].
pub fn body_force<T, const P: usize, const N: usize, const G: usize>(
    v: &T,
    f: impl Fn([f64; 1]) -> f64,
) -> Doublets
where
    T: Variable<1, 1, P, N>,
    T::Item: Element<1, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut doublets = Vec::with_capacity(N * G);
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let fg = f(elm.coordinates(ξg));
            for i in 0..N {
                doublets.push(Doublet::new(id[i], shape[i] * fg * jacobe * wg));
            }
        }
    }
    doublets
}

/// Applies a displacement constraint (essential (Dirichlet) boundary condition)
/// with the penalty method, using point elements
/// ([`crate::element::poi1::Poi1`]) to represent the constrained points.
///
/// `penalty` is the penalty parameter and `ū` the prescribed displacement.
/// `v` is a set of point elements ([`Variable<D, 1, 0, 1>`]); each point element
/// corresponds to one constrained point, where `Kᵢⱼ += α Nᵢ Nⱼ` and
/// `Fᵢ += α Nᵢ ū` are applied. Returns a `(stiffness triplets, [`Doublet`]
/// load terms)` pair.
pub fn point_displacement_penalty<T, const D: usize>(v: &T, penalty: f64, u: f64) -> TripletsAndDoublets
where
    T: Variable<D, 1, 0, 1>,
    T::Item: Element<D, 0, 1>,
{
    let mut triplets = Vec::new();
    let mut doublets = Vec::new();
    for elm in v.iter() {
        let id = elm.id();
        let shape = elm.shape([]);
        for i in 0..1 {
            for j in 0..i + 1 {
                triplets.push(Triplet::new(id[i], id[j], penalty * shape[i] * shape[j]));
            }
            doublets.push(Doublet::new(id[i], penalty * shape[i] * u));
        }
    }
    (triplets, doublets)
}

#[cfg(test)]
mod tests {
    use super::point_load;
    use crate::element::poi1::Poi1;
    use crate::node::Node;
    use crate::sparse_vector::try_new_from_doublet;

    #[test]
    fn point_load_applies_vector_components() {
        let points = vec![Poi1::new(Node::new(2, [0.0, 0.0]))];
        let doublets = point_load(&points, [3.0, 4.0]);
        let force = try_new_from_doublet(6, &doublets).unwrap();

        assert_eq!(force[4], 3.0);
        assert_eq!(force[5], 4.0);
    }
}
