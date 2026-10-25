//! Heat-conduction (heat equation) assembly operators.
//!
//! These are the thermal analogues of the bar operators in [`super::bar`]:
//! they assemble the conduction stiffness, the heat-capacity (mass) matrix,
//! and the load terms arising from volumetric heat sources, prescribed
//! boundary heat fluxes, convection (Robin) boundaries, and prescribed
//! temperatures (Dirichlet) — the latter either approximately, with the penalty
//! method ([`temperature_penalty`]), or exactly, with the Lagrange-multiplier
//! coupling assembled by [`temperature_lagrangian`].
//!
//! All operators reuse the [`crate::element::Element`] / [`Variable`] /
//! [`crate::element::IntegrationScheme`] abstractions; see
//! `book/src/theory-heat.md` for the underlying theory.

use crate::element::{Element, FullIntegration, IntegrationScheme};
use crate::operation::{
    node_coefficients, Doublets, Triplets, TripletsAndDoublets, Variable, VariableError,
};
use crate::sparse_vector::Doublet;
use faer::sparse::Triplet;

/// Assembles the thermal-conduction (diffusion) stiffness matrix and returns
/// sparse triplets for the lower triangle (including the diagonal).
///
/// `k` is the thermal conductivity (W/(m·K)): `Kᵢⱼ = ∫ k ∇Nᵢ·∇Nⱼ dΩ`, where
/// `∇Nᵢ·∇Nⱼ` is the `D`-dimensional dot product of the physical gradients
/// returned by [`Element::derivative_shape`]. For a 1D line element
/// ([`crate::element::seg2::Seg2`], `P = 1`) this differs from
/// [`super::bar::stiffness`] only in the material coefficient (`k` replaces
/// `ea`).
///
/// `v` is a set of elements (implementing [`Variable<D, 1, P, N>`]); each element
/// must implement [`FullIntegration<P, G>`].
pub fn conduction_stiffness<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    k: f64,
) -> Triplets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
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
                    let dot = (0..D).map(|d| dshape[i][d] * dshape[j][d]).sum::<f64>();
                    triplets.push(Triplet::new(
                        id[i].max(id[j]),
                        id[i].min(id[j]),
                        k * dot * jacobe * wg,
                    ));
                }
            }
        }
    }
    triplets
}

/// Assembles the heat-capacity (mass) matrix and returns sparse triplets for
/// the lower triangle (including the diagonal).
///
/// `rho_cp` is the volumetric heat capacity `ρ c_p` (J/(m³·K)). Unlike
/// [`conduction_stiffness`], this uses the shape functions themselves rather
/// than their derivatives: `Mᵢⱼ = ∫ ρ c_p Nᵢ Nⱼ dΩ`. On a 1D `Seg2` element of
/// length `L` this yields the consistent mass matrix
/// `(ρ c_p L / 6) [[2, 1], [1, 2]]`.
///
/// `v` is a set of elements (implementing [`Variable<D, 1, P, N>`]); each element
/// must implement [`FullIntegration<P, G>`].
pub fn capacity<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    rho_cp: f64,
) -> Triplets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let mut triplets: Triplets = Vec::with_capacity(N * (N + 1) / 2 * G);
    let ξ = T::Item::points();
    let w = T::Item::weights();
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            for i in 0..N {
                for j in 0..i + 1 {
                    triplets.push(Triplet::new(
                        id[i].max(id[j]),
                        id[i].min(id[j]),
                        rho_cp * shape[i] * shape[j] * jacobe * wg,
                    ));
                }
            }
        }
    }
    triplets
}

/// Applies a distributed volumetric heat source in operator form (element-wise
/// integration).
///
/// `q` is the volumetric heat-source density (W/m³) as a function of physical
/// coordinates, integrated over each element:
/// `Fᵢ = ∫ Nᵢ q dΩ = Σ_g Nᵢ(ξ_g) q(x_g) |J| w_g`. This is the thermal
/// analogue of [`super::bar::body_force`], generalized to arbitrary space
/// dimension `D` (the source closure receives `[f64; D]` physical coordinates).
///
/// Returns a list of [`Doublet`]s that the caller can merge into a load vector
/// via [`crate::sparse_vector::try_new_from_doublet`].
pub fn heat_source<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    q: impl Fn([f64; D]) -> f64,
) -> Doublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
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
            let qg = q(elm.coordinates(ξg));
            for i in 0..N {
                doublets.push(Doublet::new(id[i], shape[i] * qg * jacobe * wg));
            }
        }
    }
    doublets
}

/// Applies a natural boundary condition: a prescribed outward heat flux,
/// integrated over boundary elements of any parametric dimension (points,
/// edges, or faces).
///
/// `qn` is the prescribed outward heat flux (W/m²); each boundary element
/// contributes `Fᵢ += ∫ Nᵢ (−q_n) dΓ = Σ_g Nᵢ(ξ_g) (−q_n) |J| w_g` to the load
/// vector. The minus sign follows from the convention `k ∂T/∂n = −q_n`, so this
/// is the sign-flipped analogue of [`super::bar::point_load`], generalized to
/// arbitrary boundary dimension. Returns a list of [`Doublet`]s.
///
/// `v` is a set of boundary elements (implementing [`Variable<D, 1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a concentrated flux.
pub fn heat_flux<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    qn: f64,
) -> Doublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
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
            for i in 0..N {
                doublets.push(Doublet::new(id[i], -shape[i] * qn * jacobe * wg));
            }
        }
    }
    doublets
}

/// Applies a temperature constraint (essential (Dirichlet) boundary condition)
/// with the penalty method, integrated over boundary elements of any
/// parametric dimension (points, edges, or faces).
///
/// `penalty` is the penalty parameter and `td` a closure mapping physical
/// coordinates to the prescribed temperature `T_D`. For each boundary element,
/// the penalty term contributes `Kᵢⱼ += α ∫ Nᵢ Nⱼ dΓ` and `Fᵢ += α ∫ Nᵢ T_D dΓ`,
/// evaluated with numerical integration (the temperature is sampled at each
/// quadrature point). Returns a `(stiffness triplets, [`Doublet`] load terms)`
/// pair.
///
/// `v` is a set of boundary elements (implementing [`Variable<D, 1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a single-node constraint.
pub fn temperature_penalty<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    penalty: f64,
    td: impl Fn([f64; D]) -> f64,
) -> TripletsAndDoublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut triplets = Vec::with_capacity(N * (N + 1) / 2 * G);
    let mut doublets = Vec::with_capacity(N * G);
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let td_g = td(elm.coordinates(ξg));
            for i in 0..N {
                for j in 0..i + 1 {
                    triplets.push(Triplet::new(
                        id[i].max(id[j]),
                        id[i].min(id[j]),
                        penalty * shape[i] * shape[j] * jacobe * wg,
                    ));
                }
                doublets.push(Doublet::new(id[i], penalty * shape[i] * td_g * jacobe * wg));
            }
        }
    }
    (triplets, doublets)
}

/// Applies a temperature constraint (essential (Dirichlet) boundary condition)
/// **exactly**, through a Lagrange multiplier, integrated over boundary
/// elements of any parametric dimension (points, edges, or faces).
///
/// This is the multiplier counterpart of [`temperature_penalty`]: the constraint
/// `T = T_D` on `Γ_D` is carried by a multiplier field `λ_h = Σₐ λₐ Λₐ` whose
/// basis `Λₐ` need not coincide with the temperature shape functions `Nⱼ`. For
/// each boundary element and quadrature point the operator assembles the
/// **coupling block**
///
/// `Gₐⱼ = ∫ Λₐ Nⱼ dΓ`  (all `M × N` entries)
///
/// and the **constraint load**
///
/// `gₐ = ∫ Λₐ T_D dΓ`  (all `M` entries),
///
/// where `td` maps physical coordinates to the prescribed temperature `T_D`
/// (sampled at each quadrature point). These enter the saddle-point system
///
/// ```text
/// [ K  Gᵀ ] [ T ]   [ F ]
/// [ G  0  ] [ λ ] = [ g ]
/// ```
///
/// with `K` the conduction (`+` convection) matrix and `F` the source/flux
/// loads; the multiplier is the outward heat flux `λ = q_n` on `Γ_D`, which the
/// penalty method only approximates by `α (T_h − T_D)`.
///
/// # Two discretizations
///
/// `temperature` and `multiplier` are two element sets paired in iteration
/// order: element `k` of `multiplier` provides the multiplier basis
/// (`<L::Item as Element>::shape`) and the global multiplier dof ids
/// (`<L::Item as Element>::id`) for element `k` of `temperature`. They must
/// therefore have the same length and share the parametric domain `P` and the
/// geometry of each boundary element, so that the temperature element can supply
/// the quadrature rule (`T::Item::points` / `weights`, which is why it needs
/// [`FullIntegration<P, G>`]) and the Jacobian. The multiplier element only
/// contributes `id` and `shape`; its own geometry and Jacobian are never
/// queried, so it may be built on stand-in nodes that merely carry the
/// multiplier dof ids.
///
/// Passing the same element set twice is the usual case (the multiplier is
/// discretized on the temperature nodes) and makes `G` the boundary mass
/// matrix. A multiplier with fewer or differently shaped basis functions (for
/// example element-wise constant on each boundary edge) is expressed by a set of
/// elements whose `shape` returns the desired basis on the shared parametric
/// domain.
///
/// # Returned indices and assembly
///
/// Unlike the stiffness-like operators in this module, **every** entry of the
/// coupling block is returned: `G` sits off the diagonal of the saddle-point
/// system and is not described by a symmetric triangle, so returning only a
/// lower triangle would lose entries. Rows are multiplier dofs, columns
/// temperature dofs; the transposed block `Gᵀ` is assembled by swapping the row
/// and column of each triplet (same value).
///
/// The saddle-point matrix is symmetric but indefinite, so it must be factored
/// with sparse LU (`faer::sparse::SparseColMat::sp_lu`) rather than with a
/// Cholesky factorization. Because LU reads the whole matrix, the
/// lower-triangle triplets returned by [`conduction_stiffness`] / [`convection`]
/// must be mirrored as well (their convention is meant for Cholesky with
/// `Side::Lower`). See `book/src/theory-heat.md` §6 for the theory.
///
/// `temperature` holds boundary elements (implementing [`Variable<D, 1, P, N>`]);
/// each must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a single-node constraint,
/// where `G` degenerates to the identity and `λₐ` to the total heat flow (W)
/// leaving through that node.
///
/// # Errors
///
/// Returns [`VariableError::ElementCountMismatch`] if the two element sets hold
/// a different number of elements.
pub fn temperature_lagrangian<
    T,
    L,
    const D: usize,
    const P: usize,
    const N: usize,
    const M: usize,
    const G: usize,
>(
    temperature: &T,
    multiplier: &L,
    td: impl Fn([f64; D]) -> f64,
) -> Result<TripletsAndDoublets, VariableError>
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
    L: Variable<D, 1, P, M>,
    L::Item: Element<D, P, M>,
{
    let (n_temperature, n_multiplier) = (temperature.iter().count(), multiplier.iter().count());
    if n_temperature != n_multiplier {
        return Err(VariableError::ElementCountMismatch {
            first: n_temperature,
            second: n_multiplier,
        });
    }

    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut triplets: Triplets = Vec::with_capacity(M * N * G);
    let mut doublets: Doublets = Vec::with_capacity(M * G);
    for (elm, lelm) in temperature.iter().zip(multiplier.iter()) {
        let id = elm.id();
        let lid = lelm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let lambda_shape = lelm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let td_g = td(elm.coordinates(ξg));
            for a in 0..M {
                for j in 0..N {
                    triplets.push(Triplet::new(
                        lid[a],
                        id[j],
                        lambda_shape[a] * shape[j] * jacobe * wg,
                    ));
                }
                doublets.push(Doublet::new(lid[a], lambda_shape[a] * td_g * jacobe * wg));
            }
        }
    }
    Ok((triplets, doublets))
}

/// Applies a Robin (convective) boundary condition, integrated over boundary
/// elements of any parametric dimension (points, edges, or faces).
///
/// The Newton-cooling condition `−k ∂T/∂n = h(T − T_∞)` contributes to both
/// the stiffness matrix (`Kᵢⱼ += h ∫ Nᵢ Nⱼ dΓ`) and the load vector (`Fᵢ +=
/// h T_∞ ∫ Nᵢ dΓ`). `h` is the convective heat-transfer coefficient (W/(m²·K))
/// and `t_inf` the ambient temperature `T_∞`. Returns a `(stiffness triplets,
/// [`Doublet`] load terms)` pair.
///
/// Note the structural identity with [`temperature_penalty`]: the two operators
/// share the same algebraic form, but the penalty method enforces an essential
/// (Dirichlet) constraint whereas convection models a physical (natural-type)
/// boundary flux.
///
/// `v` is a set of boundary elements (implementing [`Variable<D, 1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a single-node term.
pub fn convection<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    h: f64,
    t_inf: f64,
) -> TripletsAndDoublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut triplets = Vec::with_capacity(N * (N + 1) / 2 * G);
    let mut doublets = Vec::with_capacity(N * G);
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            for i in 0..N {
                for j in 0..i + 1 {
                    triplets.push(Triplet::new(
                        id[i].max(id[j]),
                        id[i].min(id[j]),
                        h * shape[i] * shape[j] * jacobe * wg,
                    ));
                }
                doublets.push(Doublet::new(id[i], h * t_inf * shape[i] * jacobe * wg));
            }
        }
    }
    (triplets, doublets)
}

/// Computes the L² error `‖u_h − u‖_{L²} = (∫ (u_h − u)² dΩ)^{1/2}` between the
/// finite-element field `u_h` and the exact field `u`, element-wise with
/// numerical quadrature.
///
/// `u_h` maps a node id to its nodal value; `u` maps physical coordinates to
/// the exact value.
pub fn l2_error<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    u_h: impl Fn(usize) -> f64,
    u: impl Fn([f64; D]) -> f64,
) -> f64
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut sum = 0.0_f64;
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let uh = (0..N).map(|i| shape[i] * u_h(id[i])).sum::<f64>();
            let ue = u(elm.coordinates(ξg));
            let diff = uh - ue;
            sum += diff * diff * jacobe * wg;
        }
    }
    sum.sqrt()
}

/// Computes the H¹ semi-norm error `|u_h − u|_{H¹} = (∫ |∇u_h − ∇u|² dΩ)^{1/2}`
/// between the finite-element field `u_h` and the exact field `u`, element-wise
/// with numerical quadrature.
///
/// `u_h` maps a node id to its nodal value; `grad_u` maps physical coordinates
/// to the exact gradient.
pub fn h1_semi_error<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    u_h: impl Fn(usize) -> f64,
    grad_u: impl Fn([f64; D]) -> [f64; D],
) -> f64
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut sum = 0.0_f64;
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let dshape = elm.derivative_shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let xg = elm.coordinates(ξg);
            let grad_uh: [f64; D] =
                std::array::from_fn(|d| (0..N).map(|i| dshape[i][d] * u_h(id[i])).sum());
            let grad_ue = grad_u(xg);
            let diff2: f64 = (0..D).map(|d| (grad_uh[d] - grad_ue[d]).powi(2)).sum();
            sum += diff2 * jacobe * wg;
        }
    }
    sum.sqrt()
}

/// Physical temperature gradient at a known element and parametric point.
///
/// # Errors
///
/// Returns [`VariableError::CoefficientOutOfRange`] if the nodal temperature
/// vector is too short for a referenced node.
pub fn gradient<T, const D: usize, const P: usize, const N: usize>(
    elm: &T,
    xi: [f64; P],
    temperature: &[f64],
) -> Result<[f64; D], VariableError>
where
    T: Element<D, P, N>,
{
    let dshape = elm.derivative_shape(xi);
    let node_ids = elm.id();
    let mut grad = [0.0_f64; D];

    for a in 0..N {
        let value = node_coefficients(temperature, node_ids[a], 1)?[0];
        for (i, component) in grad.iter_mut().enumerate() {
            *component += dshape[a][i] * value;
        }
    }

    Ok(grad)
}

/// Heat flux `q = -k ∇T` at a known element and parametric point.
///
/// # Errors
///
/// Returns [`VariableError::CoefficientOutOfRange`] if the nodal temperature
/// vector is too short for a referenced node.
pub fn flux<T, const D: usize, const P: usize, const N: usize>(
    elm: &T,
    xi: [f64; P],
    temperature: &[f64],
    conductivity: f64,
) -> Result<[f64; D], VariableError>
where
    T: Element<D, P, N>,
{
    let grad = gradient(elm, xi, temperature)?;
    Ok(std::array::from_fn(|i| -conductivity * grad[i]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::poi1::Poi1;
    use crate::element::seg2::Seg2;
    use crate::node::Node;
    use crate::sparse_vector::try_new_from_doublet;
    use faer::Side;
    use faer::prelude::*;
    use faer::sparse::SparseColMat;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-5
    }

    #[test]
    fn conduction_stiffness_matches_analytic() {
        // One Seg2 element of length L = 4 with k = 2.
        // K^e = (k/L) [[1, -1], [-1, 1]] = [[1/2, -1/2], [-1/2, 1/2]].
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let triplets = conduction_stiffness(&bars, 2.0);

        let mut k00 = 0.0;
        let mut k10 = 0.0;
        let mut k11 = 0.0;
        for t in &triplets {
            match (t.row, t.col) {
                (0, 0) => k00 += t.val,
                (1, 0) => k10 += t.val,
                (1, 1) => k11 += t.val,
                _ => panic!("unexpected triplet ({}, {})", t.row, t.col),
            }
        }
        assert!(close(k00, 0.5));
        assert!(close(k10, -0.5));
        assert!(close(k11, 0.5));
    }

    #[test]
    fn capacity_matches_consistent_mass_matrix() {
        // One Seg2 element of length L = 4 with rho_cp = 1.
        // Consistent mass matrix: (rho_cp L / 6) [[2, 1], [1, 2]] = [[4/3, 2/3], [2/3, 4/3]].
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let triplets = capacity(&bars, 1.0);

        let mut m00 = 0.0;
        let mut m10 = 0.0;
        let mut m11 = 0.0;
        for t in &triplets {
            match (t.row, t.col) {
                (0, 0) => m00 += t.val,
                (1, 0) => m10 += t.val,
                (1, 1) => m11 += t.val,
                _ => panic!("unexpected triplet ({}, {})", t.row, t.col),
            }
        }
        assert!(close(m00, 4.0 / 3.0));
        assert!(close(m10, 2.0 / 3.0));
        assert!(close(m11, 4.0 / 3.0));
    }

    #[test]
    fn steady_state_dirichlet_linear_profile() {
        // Two Seg2 elements over [0, 4], k = 2, no source, T(0) = 0, T(4) = 1.
        // Steady k T'' = 0 -> linear profile, T(2) = 0.5.
        let k = 2.0;
        let bars: Vec<Seg2<1>> = vec![
            Seg2::new(Node::new(0, [0.0]), Node::new(1, [2.0])),
            Seg2::new(Node::new(1, [2.0]), Node::new(2, [4.0])),
        ];
        let left = Poi1::new(Node::new(0, [0.0]));
        let right = Poi1::new(Node::new(2, [4.0]));

        let mut triplets = conduction_stiffness(&bars, k);
        let alpha = 1e7;
        let (pt, pd) = temperature_penalty(&vec![left], alpha, |_| 0.0);
        triplets.extend(pt);
        let mut doublets = pd;
        let (pt2, pd2) = temperature_penalty(&vec![right], alpha, |_| 1.0);
        triplets.extend(pt2);
        doublets.extend(pd2);

        let kk = SparseColMat::<usize, f64>::try_new_from_triplets(3, 3, &triplets).unwrap();
        let f = try_new_from_doublet(3, &doublets).unwrap();
        let llt = kk.sp_cholesky(Side::Lower).unwrap();
        let d = llt.solve(&f);

        assert!(close(d[0], 0.0));
        assert!(close(d[1], 0.5));
        assert!(close(d[2], 1.0));
    }

    #[test]
    fn steady_state_uniform_source() {
        // Two Seg2 elements over [0, 4], k = 1, uniform source Q = 1, T(0) = T(4) = 0.
        // Analytic: T(x) = -x²/2 + 2x, so T(2) = 2.
        let k = 1.0;
        let bars: Vec<Seg2<1>> = vec![
            Seg2::new(Node::new(0, [0.0]), Node::new(1, [2.0])),
            Seg2::new(Node::new(1, [2.0]), Node::new(2, [4.0])),
        ];
        let left = Poi1::new(Node::new(0, [0.0]));
        let right = Poi1::new(Node::new(2, [4.0]));

        let mut triplets = conduction_stiffness(&bars, k);
        let alpha = 1e7;
        let (pt, pd) = temperature_penalty(&vec![left], alpha, |_| 0.0);
        triplets.extend(pt);
        let mut doublets = pd;
        let (pt2, pd2) = temperature_penalty(&vec![right], alpha, |_| 0.0);
        triplets.extend(pt2);
        doublets.extend(pd2);
        doublets.extend(heat_source(&bars, |_| 1.0));

        let kk = SparseColMat::<usize, f64>::try_new_from_triplets(3, 3, &triplets).unwrap();
        let f = try_new_from_doublet(3, &doublets).unwrap();
        let llt = kk.sp_cholesky(Side::Lower).unwrap();
        let d = llt.solve(&f);

        assert!(close(d[0], 0.0));
        assert!(close(d[1], 2.0));
        assert!(close(d[2], 0.0));
    }

    #[test]
    fn steady_state_heat_flux() {
        // One Seg2 element over [0, 4], k = 1, T(0) = 0, prescribed outward
        // flux qn = -1 (inward heat) at x = 4. Steady k T'' = 0 -> T(x) = x, T(4) = 4.
        let k = 1.0;
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let left = Poi1::new(Node::new(0, [0.0]));
        let right = Poi1::new(Node::new(1, [4.0]));

        let mut triplets = conduction_stiffness(&bars, k);
        let alpha = 1e7;
        let (pt, pd) = temperature_penalty(&vec![left], alpha, |_| 0.0);
        triplets.extend(pt);
        let mut doublets = pd;
        doublets.extend(heat_flux(&vec![right], -1.0));

        let kk = SparseColMat::<usize, f64>::try_new_from_triplets(2, 2, &triplets).unwrap();
        let f = try_new_from_doublet(2, &doublets).unwrap();
        let llt = kk.sp_cholesky(Side::Lower).unwrap();
        let d = llt.solve(&f);

        assert!(close(d[0], 0.0));
        assert!(close(d[1], 4.0));
    }

    #[test]
    fn steady_state_convection() {
        // One Seg2 element over [0, 4], k = 1, T(0) = 0, Robin at x = 4:
        // -k T'(L) = h (T(L) - T_inf) with h = 1, T_inf = 10.
        // Analytic: T(x) = h T_inf/(k + h L) x = 10/5 x = 2x, T(4) = 8.
        let k = 1.0;
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let left = Poi1::new(Node::new(0, [0.0]));
        let right = Poi1::new(Node::new(1, [4.0]));

        let mut triplets = conduction_stiffness(&bars, k);
        let alpha = 1e7;
        let (pt, pd) = temperature_penalty(&vec![left], alpha, |_| 0.0);
        triplets.extend(pt);
        let mut doublets = pd;
        let (ct, cd) = convection(&vec![right], 1.0, 10.0);
        triplets.extend(ct);
        doublets.extend(cd);

        let kk = SparseColMat::<usize, f64>::try_new_from_triplets(2, 2, &triplets).unwrap();
        let f = try_new_from_doublet(2, &doublets).unwrap();
        let llt = kk.sp_cholesky(Side::Lower).unwrap();
        let d = llt.solve(&f);

        assert!(close(d[0], 0.0));
        assert!(close(d[1], 8.0));
    }

    #[test]
    fn heat_flux_integrates_over_segment() {
        // A Seg2<1> boundary element of length L = 4 with outward flux qn.
        // Fᵢ += ∫ Nᵢ (−q_n) dx = (−q_n) (L/2) [1, 1].
        let segs: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let qn = 2.0;
        let doublets = heat_flux(&segs, qn);

        let mut f0 = 0.0;
        let mut f1 = 0.0;
        for d in &doublets {
            match d.index {
                0 => f0 += d.value,
                1 => f1 += d.value,
                _ => panic!("unexpected index {}", d.index),
            }
        }
        assert!(close(f0, -qn * 4.0 / 2.0));
        assert!(close(f1, -qn * 4.0 / 2.0));
    }

    #[test]
    fn temperature_penalty_integrates_over_segment() {
        // A Seg2<1> boundary element of length L = 4, penalty alpha, td = T_D.
        // Penalty stiffness: alpha * ∫ NᵢNⱼ dx = alpha (L/6) [[2, 1], [1, 2]].
        // Penalty load: alpha * T_D * ∫ Nᵢ dx = alpha T_D (L/2) [1, 1].
        let segs: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let alpha = 1e7;
        let td = 3.0;
        let (triplets, doublets) = temperature_penalty(&segs, alpha, |_| td);

        let mut k00 = 0.0;
        let mut k10 = 0.0;
        let mut k11 = 0.0;
        for t in &triplets {
            match (t.row, t.col) {
                (0, 0) => k00 += t.val,
                (1, 0) => k10 += t.val,
                (1, 1) => k11 += t.val,
                _ => panic!("unexpected triplet ({}, {})", t.row, t.col),
            }
        }
        assert!(close(k00, alpha * 4.0 / 6.0 * 2.0));
        assert!(close(k10, alpha * 4.0 / 6.0));
        assert!(close(k11, alpha * 4.0 / 6.0 * 2.0));

        let mut f0 = 0.0;
        let mut f1 = 0.0;
        for d in &doublets {
            match d.index {
                0 => f0 += d.value,
                1 => f1 += d.value,
                _ => panic!("unexpected index {}", d.index),
            }
        }
        assert!(close(f0, alpha * td * 4.0 / 2.0));
        assert!(close(f1, alpha * td * 4.0 / 2.0));
    }

    /// Test-only multiplier element: a constant (`P0`-like) basis on a 1D edge.
    ///
    /// It shares the parametric domain of [`Seg2`] (`P = 1`) but carries a
    /// single degree of freedom whose basis function is identically one. The
    /// dof id is carried by `dof`; the geometry stored in `ends` is only there
    /// to satisfy the trait, because [`temperature_lagrangian`] takes the
    /// quadrature and the Jacobian from the temperature element.
    struct ConstantEdge {
        dof: Node<1>,
        ends: [Node<1>; 2],
    }

    impl ConstantEdge {
        fn new(dof: usize, a: f64, b: f64) -> Self {
            Self {
                dof: Node::new(dof, [a]),
                ends: [Node::new(0, [a]), Node::new(1, [b])],
            }
        }
    }

    impl Element<1, 1, 1> for ConstantEdge {
        fn id(&self) -> [usize; 1] {
            [self.dof.id]
        }
        fn shape(&self, _parametric_coordinates: [f64; 1]) -> [f64; 1] {
            [1.0]
        }
        fn derivative_shape(&self, _parametric_coordinates: [f64; 1]) -> [[f64; 1]; 1] {
            [[0.0]]
        }
        fn vertices_coordinates(&self) -> [[f64; 1]; 1] {
            [self.ends[0].coordinates]
        }
        fn coordinates(&self, _parametric_coordinates: [f64; 1]) -> [f64; 1] {
            self.ends[0].coordinates
        }
        fn jacobe(&self, _parametric_coordinates: [f64; 1]) -> f64 {
            1.0
        }
        fn jacobe_mat(&self, _parametric_coordinates: [f64; 1]) -> [[f64; 1]; 1] {
            [[1.0]]
        }
        fn local_coordinates(&self, x: [f64; 1]) -> Option<[f64; 1]> {
            let (a, b) = (self.ends[0].coordinates[0], self.ends[1].coordinates[0]);
            ((x[0] - a) * (x[0] - b) <= 0.0).then_some([0.0])
        }
    }

    #[test]
    fn temperature_lagrangian_integrates_over_segment() {
        // A Seg2<1> boundary element of length L = 4 with T_D = x/4, and the
        // multiplier discretized on the same nodes.
        // G = ∫ NᵢNⱼ dx = (L/6) [[2, 1], [1, 2]] = [[4/3, 2/3], [2/3, 4/3]]
        // (all four entries are returned, the block is not a symmetric triangle).
        // g = ∫ Nᵢ T_D dx = G [0, 1]ᵀ for this linear T_D = [2/3, 4/3].
        let segs: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let (triplets, doublets) =
            temperature_lagrangian(&segs, &segs, |x: [f64; 1]| x[0] / 4.0).unwrap();

        let mut g = [[0.0_f64; 2]; 2];
        for t in &triplets {
            g[t.row][t.col] += t.val;
        }
        assert!(close(g[0][0], 4.0 / 3.0));
        assert!(close(g[0][1], 2.0 / 3.0));
        assert!(close(g[1][0], 2.0 / 3.0));
        assert!(close(g[1][1], 4.0 / 3.0));

        let mut f = [0.0_f64; 2];
        for d in &doublets {
            f[d.index] += d.value;
        }
        assert!(close(f[0], 2.0 / 3.0));
        assert!(close(f[1], 4.0 / 3.0));
    }

    #[test]
    fn temperature_lagrangian_mixed_basis() {
        // Same boundary edge [0, 4], but the multiplier is constant (1 dof)
        // while the temperature is linear (2 dofs):
        // G = [∫1·N₀ dΓ, ∫1·N₁ dΓ] = [L/2, L/2] = [2, 2],
        // g = ∫1·T_D dΓ = L·T_D = 12 for T_D = 3.
        let segs: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let multipliers: Vec<ConstantEdge> = vec![ConstantEdge::new(7, 0.0, 4.0)];
        let (triplets, doublets) = temperature_lagrangian(&segs, &multipliers, |_| 3.0).unwrap();

        let mut g = [0.0_f64; 2];
        for t in &triplets {
            assert_eq!(t.row, 7, "the multiplier has a single dof");
            g[t.col] += t.val;
        }
        assert!(close(g[0], 2.0));
        assert!(close(g[1], 2.0));

        let mut f = 0.0;
        for d in &doublets {
            assert_eq!(d.index, 7);
            f += d.value;
        }
        assert!(close(f, 12.0));
    }

    #[test]
    fn temperature_lagrangian_rejects_mismatched_element_sets() {
        let segs: Vec<Seg2<1>> = vec![
            Seg2::new(Node::new(0, [0.0]), Node::new(1, [2.0])),
            Seg2::new(Node::new(1, [2.0]), Node::new(2, [4.0])),
        ];
        let one: Vec<ConstantEdge> = vec![ConstantEdge::new(5, 0.0, 4.0)];
        let err = temperature_lagrangian(&segs, &one, |_| 0.0).unwrap_err();
        assert_eq!(
            err,
            VariableError::ElementCountMismatch {
                first: 2,
                second: 1
            }
        );
        let message = err.to_string();
        assert!(message.contains("cannot pair 2 elements with 1 elements"));
    }

    #[test]
    fn lagrangian_saddle_point_is_exact() {
        // Two Seg2 elements over [0, 4], k = 2, no source, T(0) = 0, T(4) = 1.
        // Steady k T'' = 0 -> T = x/4, so T(2) = 0.5 exactly (no penalty error).
        // The multipliers are the outward fluxes: at the cold end heat leaves
        // (λ = +k/4 = 0.5), at the hot end it enters (λ = -k/4 = -0.5), so
        // Σλ = 0 as required by the global heat balance.
        let k = 2.0;
        let bars: Vec<Seg2<1>> = vec![
            Seg2::new(Node::new(0, [0.0]), Node::new(1, [2.0])),
            Seg2::new(Node::new(1, [2.0]), Node::new(2, [4.0])),
        ];
        // Temperature boundary points at x = 0 and x = 4, and multiplier
        // elements carrying the constrained dof ids 3 and 4 (after the 3
        // temperature unknowns).
        let ends: [Vec<Poi1<1>>; 2] = [
            vec![Poi1::new(Node::new(0, [0.0]))],
            vec![Poi1::new(Node::new(2, [4.0]))],
        ];
        let multipliers: [Vec<Poi1<1>>; 2] = [
            vec![Poi1::new(Node::new(3, [0.0]))],
            vec![Poi1::new(Node::new(4, [4.0]))],
        ];

        // Sparse LU reads the whole matrix, so mirror the lower-triangle stiffness
        // triplets (the `Side::Lower` Cholesky convention does not apply here).
        let mut triplets: Triplets = Vec::new();
        for t in conduction_stiffness(&bars, k) {
            if t.row != t.col {
                triplets.push(Triplet::new(t.col, t.row, t.val));
            }
            triplets.push(t);
        }
        let mut doublets: Doublets = Vec::new();
        for ((temperature, multiplier), td) in ends.iter().zip(&multipliers).zip([0.0, 1.0]) {
            let (gt, gd) = temperature_lagrangian(temperature, multiplier, |_| td).unwrap();
            for d in &gd {
                doublets.push(*d);
            }
            for t in &gt {
                triplets.push(Triplet::new(t.row, t.col, t.val));
                triplets.push(Triplet::new(t.col, t.row, t.val));
            }
        }

        let a = SparseColMat::<usize, f64>::try_new_from_triplets(5, 5, &triplets).unwrap();
        let f = try_new_from_doublet(5, &doublets).unwrap();
        let d = a.sp_lu().unwrap().solve(&f);

        assert!(close(d[0], 0.0));
        assert!(close(d[1], 0.5));
        assert!(close(d[2], 1.0));
        assert!(close(d[3], 0.5));
        assert!(close(d[4], -0.5));
        assert!(close(d[3] + d[4], 0.0));
    }

    #[test]
    fn convection_integrates_over_segment() {
        // A Seg2<1> boundary element of length L = 4, h, t_inf = T_∞.
        // Stiffness: h * ∫ NᵢNⱼ dx = h (L/6) [[2, 1], [1, 2]].
        // Load: h * T_∞ * ∫ Nᵢ dx = h T_∞ (L/2) [1, 1].
        let segs: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let h = 5.0;
        let t_inf = 20.0;
        let (triplets, doublets) = convection(&segs, h, t_inf);

        let mut k00 = 0.0;
        let mut k10 = 0.0;
        let mut k11 = 0.0;
        for t in &triplets {
            match (t.row, t.col) {
                (0, 0) => k00 += t.val,
                (1, 0) => k10 += t.val,
                (1, 1) => k11 += t.val,
                _ => panic!("unexpected triplet ({}, {})", t.row, t.col),
            }
        }
        assert!(close(k00, h * 4.0 / 6.0 * 2.0));
        assert!(close(k10, h * 4.0 / 6.0));
        assert!(close(k11, h * 4.0 / 6.0 * 2.0));

        let mut f0 = 0.0;
        let mut f1 = 0.0;
        for d in &doublets {
            match d.index {
                0 => f0 += d.value,
                1 => f1 += d.value,
                _ => panic!("unexpected index {}", d.index),
            }
        }
        assert!(close(f0, h * t_inf * 4.0 / 2.0));
        assert!(close(f1, h * t_inf * 4.0 / 2.0));
    }

    #[test]
    fn conduction_stiffness_segment_in_2d() {
        // A Seg2<2> line element in 2D space, nodes (0,0) and (3,4), length L = 5.
        // K^e = (k/L) [[1, -1], [-1, 1]] = [[0.2, -0.2], [-0.2, 0.2]] for k = 1.
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 4.0]),
        )];
        let triplets = conduction_stiffness(&bars, 1.0);

        let mut k00 = 0.0;
        let mut k10 = 0.0;
        let mut k11 = 0.0;
        for t in &triplets {
            match (t.row, t.col) {
                (0, 0) => k00 += t.val,
                (1, 0) => k10 += t.val,
                (1, 1) => k11 += t.val,
                _ => panic!("unexpected triplet ({}, {})", t.row, t.col),
            }
        }
        assert!(close(k00, 0.2));
        assert!(close(k10, -0.2));
        assert!(close(k11, 0.2));
    }

    #[test]
    fn heat_source_segment_in_2d() {
        // A Seg2<2> line element of length L = 5 with uniform source q = 1.
        // Fᵢ += ∫ Nᵢ q dx = q (L/2) [1, 1] = [2.5, 2.5].
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [3.0, 4.0]),
        )];
        let doublets = heat_source(&bars, |_| 1.0);

        let mut f0 = 0.0;
        let mut f1 = 0.0;
        for d in &doublets {
            match d.index {
                0 => f0 += d.value,
                1 => f1 += d.value,
                _ => panic!("unexpected index {}", d.index),
            }
        }
        assert!(close(f0, 2.5));
        assert!(close(f1, 2.5));
    }

    #[test]
    fn l2_and_h1_errors_analytic() {
        // One Seg2<1> element over [0, 4], u_h = 0, exact u = x.
        // L2² = ∫ x² dx = 64/3, H1² = ∫ 1² dx = 4.
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let uh = |_id: usize| 0.0;
        let l2 = l2_error(&bars, uh, |x: [f64; 1]| x[0]);
        let h1 = h1_semi_error(&bars, uh, |_: [f64; 1]| [1.0]);
        assert!(close(l2, (64.0_f64 / 3.0).sqrt()));
        assert!(close(h1, 2.0));
    }

    #[test]
    fn gradient_and_flux_of_linear_temperature() {
        let bars = [Seg2::new(Node::new(0, [0.0]), Node::new(1, [2.0]))];
        let temperature = [0.0, 2.0]; // T = x
        let grad = gradient(&bars[0], [0.0], &temperature).unwrap();
        assert_eq!(grad, [1.0]);
        let q = flux(&bars[0], [0.0], &temperature, 3.0).unwrap();
        assert_eq!(q, [-3.0]);
    }
}
