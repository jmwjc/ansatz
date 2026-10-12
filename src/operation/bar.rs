use crate::element::Element;
use crate::element::FullIntegration;
use crate::element::IntegrationScheme;
use crate::operation::Variable;
use crate::sparse_vector::Doublet;
use faer::sparse::Triplet;

/// Assembles the stiffness matrix and returns sparse triplets for the lower
/// triangle (including the diagonal).
///
/// `ea` is the axial stiffness (Young's modulus × cross-sectional area).
/// `v` is a set of 1D elements (implementing [`Variable<1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`].
pub fn stiffness<T, const P: usize, const N: usize, const G: usize>(
    v: &T,
    ea: f64,
) -> Vec<Triplet<usize, usize, f64>>
where
    T: Variable<1, P, N>,
    T::Item: Element<1, P, N> + FullIntegration<P, G>,
{
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::with_capacity(N * (N + 1) / 2 * G);
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
                        id[i],
                        id[j],
                        ea * dshape[i] * dshape[j] * jacobe * wg,
                    ));
                }
            }
        }
    }
    triplets
}

/// Applies a natural boundary condition (concentrated force), using point
/// elements ([`crate::element::poi1::Poi1`]) to represent the load points.
///
/// `p` is the force magnitude and `v` is a set of point elements
/// ([`Variable<D, 0, 1>`]); each point element corresponds to one load point,
/// where `Fᵢ += Nᵢ P` is applied. Returns a list of [`Doublet`]s that the
/// caller can merge into a load vector via
/// [`crate::sparse_vector::try_new_from_doublet`].
pub fn traction<T, const D: usize>(v: &T, p: f64) -> Vec<Doublet<usize, f64>>
where
    T: Variable<D, 0, 1>,
    T::Item: Element<D, 0, 1>,
{
    let mut doublets = Vec::new();
    for elm in v.iter() {
        let id = elm.id();
        let shape = elm.shape([]);
        for i in 0..1 {
            doublets.push(Doublet::new(id[i], shape[i] * p));
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
) -> Vec<Doublet<usize, f64>>
where
    T: Variable<1, P, N>,
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
/// `v` is a set of point elements ([`Variable<D, 0, 1>`]); each point element
/// corresponds to one constrained point, where `Kᵢⱼ += α Nᵢ Nⱼ` and
/// `Fᵢ += α Nᵢ ū` are applied. Returns a `(stiffness triplets, [`Doublet`]
/// load terms)` pair.
pub fn displacement_penalty<T, const D: usize>(
    v: &T,
    penalty: f64,
    u: f64,
) -> (Vec<Triplet<usize, usize, f64>>, Vec<Doublet<usize, f64>>)
where
    T: Variable<D, 0, 1>,
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
