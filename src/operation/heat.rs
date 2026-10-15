//! Heat-conduction (heat equation) assembly operators.
//!
//! These are the thermal analogues of the bar operators in [`super::bar`]:
//! they assemble the conduction stiffness, the heat-capacity (mass) matrix,
//! and the load terms arising from volumetric heat sources, prescribed
//! boundary heat fluxes, convection (Robin) boundaries, and prescribed
//! temperatures (Dirichlet, via the penalty method).
//!
//! All operators reuse the [`crate::element::Element`] / [`Variable`] /
//! [`crate::element::IntegrationScheme`] abstractions; see
//! `book/src/theory-heat.md` for the underlying theory.

use crate::element::{Element, FullIntegration, IntegrationScheme};
use crate::operation::{Doublets, Triplets, TripletsAndDoublets, Variable};
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
/// `v` is a set of elements (implementing [`Variable<D, P, N>`]); each element
/// must implement [`FullIntegration<P, G>`].
pub fn conduction_stiffness<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    k: f64,
) -> Triplets
where
    T: Variable<D, P, N>,
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
/// `v` is a set of elements (implementing [`Variable<D, P, N>`]); each element
/// must implement [`FullIntegration<P, G>`].
pub fn capacity<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    rho_cp: f64,
) -> Triplets
where
    T: Variable<D, P, N>,
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
    T: Variable<D, P, N>,
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
/// is the sign-flipped analogue of [`super::bar::traction`], generalized to
/// arbitrary boundary dimension. Returns a list of [`Doublet`]s.
///
/// `v` is a set of boundary elements (implementing [`Variable<D, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a concentrated flux.
pub fn heat_flux<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    qn: f64,
) -> Doublets
where
    T: Variable<D, P, N>,
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
/// `v` is a set of boundary elements (implementing [`Variable<D, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a single-node constraint.
pub fn temperature_penalty<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    penalty: f64,
    td: impl Fn([f64; D]) -> f64,
) -> TripletsAndDoublets
where
    T: Variable<D, P, N>,
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
/// `v` is a set of boundary elements (implementing [`Variable<D, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a single-node term.
pub fn convection<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    h: f64,
    t_inf: f64,
) -> TripletsAndDoublets
where
    T: Variable<D, P, N>,
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
    T: Variable<D, P, N>,
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
pub fn h1_semierror<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    u_h: impl Fn(usize) -> f64,
    grad_u: impl Fn([f64; D]) -> [f64; D],
) -> f64
where
    T: Variable<D, P, N>,
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
        let f = try_new_from_doublet(3, &doublets);
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
        let f = try_new_from_doublet(3, &doublets);
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
        let f = try_new_from_doublet(2, &doublets);
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
        let f = try_new_from_doublet(2, &doublets);
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
        let h1 = h1_semierror(&bars, uh, |_: [f64; 1]| [1.0]);
        assert!(close(l2, (64.0_f64 / 3.0).sqrt()));
        assert!(close(h1, 2.0));
    }
}
