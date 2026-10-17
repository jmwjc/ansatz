//! Linear-elasticity assembly operators.
//!
//! These are the vector-field (displacement) analogues of the bar operators in
//! [`super::bar`]: they assemble the elastic stiffness, the volumetric body
//! force, the surface traction, and the prescribed-displacement penalty terms.
//!
//! The displacement field is vector-valued, so every operator is expressed over
//! [`Variable<D, 1, P, N>`] (tensor order `O = 1`): node `a` owns `D` global
//! degrees of freedom `a·D + p` (`p = 0..D`).
//!
//! Following the tensor-notation convention of this crate, the physics inside
//! the operators is written directly in index (summation) form rather than by
//! assembling explicit Voigt `B`/`D` matrices. See
//! `book/src/theory-elasticity.md` for the underlying theory.

use crate::element::{Element, FullIntegration, IntegrationScheme};
use crate::operation::{Doublets, Triplets, TripletsAndDoublets, Variable};
use crate::sparse_vector::Doublet;
use faer::sparse::Triplet;

/// Converts the engineering constants — Young's modulus `e` and Poisson's ratio
/// `nu` — into the Lamé constants `(λ, μ)` of 3D / plane-strain isotropic
/// linear elasticity:
///
/// `λ = Eν / ((1 + ν)(1 − 2ν))`, `μ = E / (2(1 + ν))`.
pub fn lame_constants(e: f64, nu: f64) -> (f64, f64) {
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    (lambda, mu)
}

/// Effective first Lamé constant `λ* = Eν / (1 − ν²)` for plane stress.
///
/// In plane stress the out-of-plane stress is condensed out, which changes the
/// in-plane volumetric coupling from `λ` to `λ*` while leaving the shear
/// modulus `μ = E / (2(1 + ν))` unchanged. Pass this value as `lambda` to
/// [`stiffness`] for plane-stress problems.
pub fn plane_stress_lambda(e: f64, nu: f64) -> f64 {
    e * nu / (1.0 - nu * nu)
}

/// Isotropic linear-elastic stiffness kernel (integrand density): the
/// bilinear-form value `a(Nᵢ eₐ, Nⱼ e_b)` per unit measure, before the `|J| w_g`
/// quadrature weight is applied.
///
/// In index (summation) form,
///
/// `λ (∂ₐNᵢ)(∂_b Nⱼ) + μ ( δ_ab (∂ₖNᵢ)(∂ₖNⱼ) + (∂_b Nᵢ)(∂ₐNⱼ) )`
///
/// where `∂ₖNᵢ = dshape[i][k]` is the physical gradient returned by
/// [`Element::derivative_shape`], `i, j` are nodes, `a, b` are displacement
/// components, and `k` is summed over the `D` spatial directions. This is the
/// isotropic form `ε(u):C:ε(v) = λ(∇·u)(∇·v) + μ(∇u:∇v + ∇u:∇vᵀ)` evaluated
/// on the single-component test functions `u = Nᵢ eₐ`, `v = Nⱼ e_b`.
fn isotropic_kernel<const D: usize, const N: usize>(
    dshape: [[f64; D]; N],
    lambda: f64,
    mu: f64,
    i: usize,
    j: usize,
    a: usize,
    b: usize,
) -> f64 {
    let grad_dot = (0..D)
        .map(|k| dshape[i][k] * dshape[j][k])
        .sum::<f64>();
    let delta = if a == b { 1.0 } else { 0.0 };
    lambda * dshape[i][a] * dshape[j][b]
        + mu * (delta * grad_dot + dshape[i][b] * dshape[j][a])
}

/// Assembles the linear-elastic stiffness matrix and returns sparse triplets for
/// the lower triangle (including the diagonal).
///
/// `lambda` and `mu` are the Lamé constants. The element contribution is
/// written in tensor (summation) form over repeated indices:
///
/// `K_{(a,p),(b,q)} = ∫ [ λ (∂ₚNₐ)(∂qN_b) + μ ( δ_pq (∂ₖNₐ)(∂ₖN_b) + (∂qNₐ)(∂ₚN_b) ) ] dΩ`
///
/// where `∂ₖNₐ = dshape[a][k]` is the physical gradient returned by
/// [`Element::derivative_shape`], `a, b` are nodes, `p, q` are displacement
/// components, and `k` is summed over the `D` spatial directions. This is the
/// isotropic form `ε(u) : C : ε(v) = λ (∇·u)(∇·v) + μ (∇u:∇v + ∇u:∇vᵀ)` in
/// components, and holds for any spatial dimension `D` (2D plane strain, full
/// 3D, and shell-like elements embedded in 3D). For plane stress, substitute
/// [`plane_stress_lambda`] for `lambda`; `mu` is unchanged.
///
/// `v` is a set of elements (implementing [`Variable<D, 1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`].
pub fn stiffness<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    lambda: f64,
    mu: f64,
) -> Triplets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let nd = D * N; // DOFs per element (node-major: dof = a·D + p)
    let mut triplets: Triplets = Vec::with_capacity(nd * (nd + 1) / 2 * G);
    let ξ = T::Item::points();
    let w = T::Item::weights();
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let dshape = elm.derivative_shape(ξg);
            let jacobe = elm.jacobe(ξg);
            // Node-block assembly, with the diagonal node blocks (i, i)
            // handled in a separate loop: each diagonal block emits the lower
            // triangle of its D×D block (b ≤ a, so row ≥ col holds), while each
            // off-diagonal block (i, j), j < i, is emitted in full.
            for i in 0..N {
                // Diagonal node block (i, i): lower triangle (b ≤ a).
                for a in 0..D {
                    for b in 0..=a {
                        let g_a = id[i] * D + a;
                        let g_b = id[i] * D + b;
                        triplets.push(Triplet::new(
                            g_a,
                            g_b,
                            isotropic_kernel(dshape, lambda, mu, i, i, a, b) * jacobe * wg,
                        ));
                    }
                }
                // Off-diagonal node blocks (i, j), j < i: full D×D block.
                for j in 0..i {
                    for a in 0..D {
                        for b in 0..D {
                            let g_a = id[i] * D + a;
                            let g_b = id[j] * D + b;
                            triplets.push(Triplet::new(
                                g_a.max(g_b),
                                g_a.min(g_b),
                                isotropic_kernel(dshape, lambda, mu, i, j, a, b) * jacobe * wg,
                            ));
                        }
                    }
                }
            }
        }
    }
    triplets
}

/// Applies a distributed volumetric body force, in operator form (element-wise
/// integration).
///
/// `f` maps physical coordinates to the body-force density vector `f(x) ∈ R^D`
/// (force per unit volume). Each element contributes
/// `F_{a,p} = ∫ Nₐ f_p dΩ = Σ_g Nₐ(ξ_g) f_p(x_g) |J| w_g` for each node `a` and
/// component `p`.
///
/// Returns a list of [`Doublet`]s that the caller can merge into a load vector
/// via [`crate::sparse_vector::try_new_from_doublet`].
///
/// `v` is a set of elements (implementing [`Variable<D, 1, P, N>`]); each
/// element must implement [`FullIntegration<P, G>`].
pub fn body_force<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    f: impl Fn([f64; D]) -> [f64; D],
) -> Doublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut doublets = Vec::with_capacity(D * N * G);
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let fg = f(elm.coordinates(ξg));
            for a in 0..N {
                for (p, &fp) in fg.iter().enumerate() {
                    // Vector load component: F_{a,p} += Nₐ f_p |J| w.
                    doublets.push(Doublet::new(id[a] * D + p, shape[a] * fp * jacobe * wg));
                }
            }
        }
    }
    doublets
}

/// Applies a natural boundary condition: a prescribed surface traction,
/// integrated over boundary elements of any parametric dimension (points,
/// edges, or faces).
///
/// `t` maps physical coordinates to the traction vector `t̄(x) ∈ R^D` (force
/// per unit surface). Each boundary element contributes
/// `F_{a,p} = ∫ Nₐ t̄_p dΓ = Σ_g Nₐ(ξ_g) t̄_p(x_g) |J| w_g`. This is the
/// vector-valued analogue of [`super::bar::traction`], generalized to arbitrary
/// boundary dimension; for a point element
/// ([`crate::element::poi1::Poi1`]) it reduces to a concentrated force.
///
/// Returns a list of [`Doublet`]s that the caller can merge into a load vector
/// via [`crate::sparse_vector::try_new_from_doublet`].
///
/// `v` is a set of boundary elements (implementing [`Variable<D, 1, P, N>`]);
/// each element must implement [`FullIntegration<P, G>`].
pub fn traction<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    t: impl Fn([f64; D]) -> [f64; D],
) -> Doublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut doublets = Vec::with_capacity(D * N * G);
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let tg = t(elm.coordinates(ξg));
            for a in 0..N {
                for (p, &tp) in tg.iter().enumerate() {
                    // Surface-load component: F_{a,p} += Nₐ t̄_p |J| w.
                    doublets.push(Doublet::new(id[a] * D + p, shape[a] * tp * jacobe * wg));
                }
            }
        }
    }
    doublets
}

/// Applies a displacement constraint (essential (Dirichlet) boundary condition)
/// with the penalty method, integrated over boundary elements of any
/// parametric dimension (points, edges, or faces).
///
/// `penalty` is the penalty parameter and `u` a closure mapping physical
/// coordinates to the prescribed displacement vector `ū(x) ∈ R^D`. For each
/// boundary element, the penalty term contributes, in tensor (summation) form,
///
/// `K_{(a,p),(b,q)} += α ∫ Nₐ N_b δ_pq dΓ`, `F_{a,p} += α ∫ Nₐ ū_p dΓ`,
///
/// evaluated with numerical integration (the prescribed displacement is sampled
/// at each quadrature point). Returns a `(stiffness triplets, [`Doublet`] load
/// terms)` pair.
///
/// `v` is a set of boundary elements (implementing [`Variable<D, 1, P, N>`]);
/// each element must implement [`FullIntegration<P, G>`]. For a point element
/// ([`crate::element::poi1::Poi1`]) this reduces to a single-node vector
/// constraint.
pub fn displacement_penalty<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    penalty: f64,
    u: impl Fn([f64; D]) -> [f64; D],
) -> TripletsAndDoublets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + FullIntegration<P, G>,
{
    let ξ = T::Item::points();
    let w = T::Item::weights();
    let mut triplets = Vec::with_capacity(N * (N + 1) / 2 * D * G);
    let mut doublets = Vec::with_capacity(D * N * G);
    for elm in v.iter() {
        let id = elm.id();
        for g in 0..G {
            let ξg = ξ[g];
            let wg = w[g];
            let shape = elm.shape(ξg);
            let jacobe = elm.jacobe(ξg);
            let ug = u(elm.coordinates(ξg));
            for a in 0..N {
                for b in 0..=a {
                    for p in 0..D {
                        // K_{(a,p),(b,p)} += α Nₐ N_b δ_pp |J| w.
                        let g_a = id[a] * D + p;
                        let g_b = id[b] * D + p;
                        triplets.push(Triplet::new(
                            g_a.max(g_b),
                            g_a.min(g_b),
                            penalty * shape[a] * shape[b] * jacobe * wg,
                        ));
                    }
                }
                for (p, &up) in ug.iter().enumerate() {
                    // F_{a,p} += α Nₐ ū_p |J| w.
                    doublets.push(Doublet::new(
                        id[a] * D + p,
                        penalty * shape[a] * up * jacobe * wg,
                    ));
                }
            }
        }
    }
    (triplets, doublets)
}

/// Computes the L² error `‖u_h − u‖_{L²} = (∫ |u_h − u|² dΩ)^{1/2}` between the
/// finite-element displacement field `u_h` and the exact displacement field
/// `u`, element-wise with numerical quadrature.
///
/// `u_h` maps a `(node id, component)` pair to the nodal displacement value;
/// `u` maps physical coordinates to the exact displacement vector. The pointwise
/// difference uses the Euclidean norm `|·|` over the `D` components.
///
/// `v` is a set of elements (implementing [`Variable<D, 1, P, N>`]); each element
/// must implement [`FullIntegration<P, G>`].
pub fn l2_error<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    u_h: impl Fn(usize, usize) -> f64,
    u: impl Fn([f64; D]) -> [f64; D],
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
            let exact = u(elm.coordinates(ξg));
            let mut diff2 = 0.0_f64;
            for (i, &exact_i) in exact.iter().enumerate() {
                let uh = (0..N).map(|a| shape[a] * u_h(id[a], i)).sum::<f64>();
                diff2 += (uh - exact_i).powi(2);
            }
            sum += diff2 * jacobe * wg;
        }
    }
    sum.sqrt()
}

/// Computes the energy-norm error of the displacement field,
/// `‖u_h − u‖_E = (∫ (ε_h − ε) : C : (ε_h − ε) dΩ)^{1/2}`, between the
/// finite-element field `u_h` and the exact field `u`, element-wise with
/// numerical quadrature.
///
/// The finite-element strain `ε_h = sym ∇u_h` is evaluated at each quadrature
/// point from the physical gradients returned by [`Element::derivative_shape`];
/// the exact strain tensor is supplied by `strain`. Writing `Δ = ε_h − ε` and
/// using the isotropic constitutive law, the integrand is evaluated in tensor
/// (summation) form as
///
/// `Δ : C : Δ = λ (Δ_kk)² + 2μ (Δ_ij Δ_ij)`
///
/// with `λ`, `μ` the Lamé constants, so the result is a genuine energy norm. For
/// plane stress pass `plane_stress_lambda(E, ν)` as `lambda` (and leave `mu`
/// unchanged), together with the in-plane exact strain components.
///
/// `v` is a set of elements (implementing [`Variable<D, 1, P, N>`]); each element
/// must implement [`FullIntegration<P, G>`].
pub fn energy_error<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    u_h: impl Fn(usize, usize) -> f64,
    strain: impl Fn([f64; D]) -> [[f64; D]; D],
    lambda: f64,
    mu: f64,
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
            // Nodal displacement values of this element.
            let un: [[f64; D]; N] =
                std::array::from_fn(|a| std::array::from_fn(|i| u_h(id[a], i)));
            // Finite-element strain ε_h = sym ∇u_h (components over D×D).
            let eh: [[f64; D]; D] = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    0.5 * (0..N)
                        .map(|a| dshape[a][j] * un[a][i] + dshape[a][i] * un[a][j])
                        .sum::<f64>()
                })
            });
            let ee = strain(elm.coordinates(ξg));
            // Δ = ε_h − ε, then Δ : C : Δ = λ (Δ_kk)² + 2μ (Δ_ij Δ_ij).
            let mut trace = 0.0_f64;
            let mut contraction = 0.0_f64;
            for i in 0..D {
                trace += eh[i][i] - ee[i][i];
                for j in 0..D {
                    let delta = eh[i][j] - ee[i][j];
                    contraction += delta * delta;
                }
            }
            sum += (lambda * trace * trace + 2.0 * mu * contraction) * jacobe * wg;
        }
    }
    sum.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::poi1::Poi1;
    use crate::element::seg2::Seg2;
    use crate::element::tri3::Tri3;
    use crate::node::Node;
    use crate::sparse_vector::try_new_from_doublet;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-8
    }

    /// Voigt `B` matrix (3 × 2) of a single node for a 2D displacement field.
    fn voigt_b(g: [f64; 2]) -> [[f64; 2]; 3] {
        [[g[0], 0.0], [0.0, g[1]], [g[1], g[0]]]
    }

    /// Reference 2D plane-strain stiffness block `B_aᵀ D B_b` (2 × 2), computed
    /// independently in Voigt form, for cross-checking the tensor-form operator.
    fn voigt_stiffness(ga: [f64; 2], gb: [f64; 2], d: [[f64; 3]; 3]) -> [[f64; 2]; 2] {
        let ba = voigt_b(ga);
        let bb = voigt_b(gb);
        let mut k = [[0.0; 2]; 2];
        for (p, kp) in k.iter_mut().enumerate() {
            for (q, kpq) in kp.iter_mut().enumerate() {
                let mut s = 0.0;
                for (r, ba_r) in ba.iter().enumerate() {
                    for (c, bb_c) in bb.iter().enumerate() {
                        s += ba_r[p] * d[r][c] * bb_c[q];
                    }
                }
                *kpq = s;
            }
        }
        k
    }

    #[test]
    fn isotropic_kernel_matches_reference() {
        // For λ = μ = 1 and the Tri3 gradient ∇N₁ = (-1, -1):
        //   ∂ₖN₁∂ₖN₁ = 2, so the (0,0) diagonal entry is
        //   λ·(-1)(-1) + μ(1·2 + (-1)(-1)) = 1 + 3 = 4,
        // and the (0,1) entry is λ·(-1)(-1) + μ(0·2 + (-1)(-1)) = 1 + 1 = 2.
        let dshape = [[-1.0_f64, -1.0]; 1];
        assert!(close(isotropic_kernel(dshape, 1.0, 1.0, 0, 0, 0, 0), 4.0));
        assert!(close(isotropic_kernel(dshape, 1.0, 1.0, 0, 0, 0, 1), 2.0));
    }

    #[test]
    fn stiffness_reduces_to_bar_in_1d() {
        // One Seg2<1> element of length L = 4, λ = 1, μ = 2.
        // In 1D the isotropic law gives E = λ + 2μ = 5, so
        // K^e = (E/L) [[1, -1], [-1, 1]] = [[1.25, -1.25], [-1.25, 1.25]].
        let bars: Vec<Seg2<1>> = vec![Seg2::new(Node::new(0, [0.0]), Node::new(1, [4.0]))];
        let triplets = stiffness(&bars, 1.0, 2.0);

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
        assert!(close(k00, 1.25));
        assert!(close(k10, -1.25));
        assert!(close(k11, 1.25));
    }

    #[test]
    fn stiffness_matches_voigt_plane_strain_tri3() {
        // Unit right triangle (0,0), (1,0), (0,1); area = 1/2, jacobe = 1.
        // Constant-strain element, so the tensor-form operator must reproduce
        // the Voigt reference K_ab = B_aᵀ D B_b · area exactly.
        let tris: Vec<Tri3<2>> = vec![Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
            Node::new(2, [0.0, 1.0]),
        )];
        let lambda = 1.0;
        let mu = 1.0;
        let triplets = stiffness(&tris, lambda, mu);

        let d = [
            [lambda + 2.0 * mu, lambda, 0.0],
            [lambda, lambda + 2.0 * mu, 0.0],
            [0.0, 0.0, mu],
        ];
        let grads: [[f64; 2]; 3] = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
        let area = 0.5;

        let mut reference = [[0.0; 6]; 6];
        for a in 0..3 {
            for b in 0..3 {
                let k = voigt_stiffness(grads[a], grads[b], d);
                for p in 0..2 {
                    for q in 0..2 {
                        reference[a * 2 + p][b * 2 + q] = k[p][q] * area;
                    }
                }
            }
        }

        let mut assembled = [[0.0; 6]; 6];
        for t in &triplets {
            assembled[t.row][t.col] += t.val;
        }
        for r in 0..6 {
            for c in 0..=r {
                assert!(
                    close(assembled[r][c], reference[r][c]),
                    "mismatch at ({r}, {c}): {} != {}",
                    assembled[r][c],
                    reference[r][c]
                );
            }
        }
    }

    #[test]
    fn body_force_maps_vector_components() {
        // One Seg2<2> element from (0,0) to (4,0), length L = 4, uniform
        // body force f = (1, 2). F_{a,p} = f_p · L/2 = (2, 4) per node.
        let bars: Vec<Seg2<2>> = vec![Seg2::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [4.0, 0.0]),
        )];
        let doublets = body_force(&bars, |_| [1.0, 2.0]);

        let mut f = [0.0; 4];
        for d in &doublets {
            f[d.index] += d.value;
        }
        assert!(close(f[0], 2.0));
        assert!(close(f[1], 4.0));
        assert!(close(f[2], 2.0));
        assert!(close(f[3], 4.0));
    }

    #[test]
    fn traction_reduces_to_concentrated_force() {
        // A single point element at node 5, prescribed traction t = (3, 4).
        // F_{5D+p} += N · t_p = t_p.
        let pts: Vec<Poi1<2>> = vec![Poi1::new(Node::new(5, [0.0, 0.0]))];
        let doublets = traction(&pts, |_| [3.0, 4.0]);

        let f = try_new_from_doublet(12, &doublets);
        assert!(close(f[10], 3.0)); // node 5, component 0
        assert!(close(f[11], 4.0)); // node 5, component 1
    }

    #[test]
    fn displacement_penalty_maps_vector_components() {
        // A single point element at node 3, prescribed displacement u = (1, 2),
        // penalty alpha. K_{3D+p,3D+p} += alpha, F_{3D+p} += alpha · u_p.
        let pts: Vec<Poi1<2>> = vec![Poi1::new(Node::new(3, [0.0, 0.0]))];
        let alpha = 1e7;
        let (triplets, doublets) = displacement_penalty(&pts, alpha, |_| [1.0, 2.0]);

        let mut rows = Vec::new();
        for t in &triplets {
            assert_eq!(t.row, t.col);
            assert!(close(t.val, alpha));
            rows.push(t.row);
        }
        rows.sort();
        // Node 3 owns DOFs 6 and 7 (id·D + p).
        assert_eq!(rows, vec![6, 7]);

        let f = try_new_from_doublet(8, &doublets);
        assert!(close(f[6], alpha));
        assert!(close(f[7], 2.0 * alpha));
    }

    #[test]
    fn plane_stress_lamé_round_trip() {
        // For E = 200e9, ν = 0.3 the conversions must match the known formulas.
        let e = 200e9;
        let nu = 0.3;
        let (lambda, mu) = lame_constants(e, nu);
        assert!(close(lambda, e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))));
        assert!(close(mu, e / (2.0 * (1.0 + nu))));
        assert!(close(plane_stress_lambda(e, nu), e * nu / (1.0 - nu * nu)));
    }

    /// Unit right triangle (0,0), (1,0), (0,1), area 1/2, forming a single-
    /// element mesh used by the error-operator tests.
    fn unit_triangle() -> Vec<Tri3<2>> {
        vec![Tri3::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
            Node::new(2, [0.0, 1.0]),
        )]
    }

    #[test]
    fn l2_error_is_zero_for_exactly_represented_field() {
        // A linear displacement field u = (1 + 2x, 3 + x + y) is reproduced
        // exactly by the constant-strain triangle, so the error must vanish.
        let tris = unit_triangle();
        let exact = |x: [f64; 2]| [1.0 + 2.0 * x[0], 3.0 + x[0] + x[1]];
        let nodal = [[1.0, 3.0], [3.0, 4.0], [1.0, 4.0]]; // exact at the 3 nodes
        let l2 = l2_error(&tris, |id, p| nodal[id][p], exact);
        assert!(l2 < 1e-12, "L2 error should vanish, got {l2:e}");
    }

    #[test]
    fn l2_error_matches_analytic_constant_field() {
        // u_h = 0, exact u = (1, 0): ∫|u|² dΩ = area = 1/2, so the L² error is
        // sqrt(1/2).
        let tris = unit_triangle();
        let l2 = l2_error(&tris, |_, _| 0.0, |_| [1.0, 0.0]);
        assert!(close(l2, 0.5_f64.sqrt()));
    }

    #[test]
    fn energy_error_is_zero_for_exactly_represented_field() {
        // Same linear field: its strain is constant, so the energy-norm error
        // vanishes for any λ, μ.
        let tris = unit_triangle();
        let nodal = [[1.0, 3.0], [3.0, 4.0], [1.0, 4.0]];
        let exact_strain = |_: [f64; 2]| [[2.0, 0.5], [0.5, 1.0]];
        let e = energy_error(&tris, |id, p| nodal[id][p], exact_strain, 1.0, 2.0);
        assert!(e < 1e-12, "energy error should vanish, got {e:e}");
    }

    #[test]
    fn energy_error_matches_analytic_uniaxial_field() {
        // u_h = 0, exact u = (x, 0), so ε = [[1, 0], [0, 0]] and
        // Δ : C : Δ = λ (1)² + 2μ (1) = λ + 2μ. Over area 1/2 the energy-norm
        // error is sqrt((λ + 2μ)/2).
        let tris = unit_triangle();
        let exact_strain = |_: [f64; 2]| [[1.0, 0.0], [0.0, 0.0]];
        let e = energy_error(&tris, |_, _| 0.0, exact_strain, 3.0, 2.0);
        assert!(close(e, ((3.0 + 2.0 * 2.0) / 2.0_f64).sqrt()));
    }

    #[test]
    fn cantilever_beam_matches_analytic_solution() {
        // End-to-end verification on the GMSH plane-stress cantilever mesh:
        // the L² and energy-norm errors of the finite-element solution against
        // the exact Saint-Venant solution must be small. See
        // `examples/cantilever.rs` and `book/src/example-cantilever-model.md`.
        use crate::element::quad4::Quad4;
        use crate::io::gmsh::read_gmsh;
        use faer::Side;
        use faer::prelude::*;
        use faer::sparse::SparseColMat;

        const L: f64 = 8.0;
        const H: f64 = 1.0;
        const E: f64 = 1000.0;
        const NU: f64 = 0.3;
        const P: f64 = 1.0;
        let i = H.powi(3) / 12.0;

        let exact = |x: [f64; 2]| -> [f64; 2] {
            let (xx, y) = (x[0], x[1]);
            let u = -(P * y / (6.0 * E * i))
                * ((6.0 * L - 3.0 * xx) * xx + (2.0 + NU) * (y * y - H * H / 4.0));
            let v = (P / (6.0 * E * i))
                * (3.0 * NU * y * y * (L - xx)
                    + (4.0 + 5.0 * NU) * H * H * xx / 4.0
                    + (3.0 * L - xx) * xx * xx);
            [u, v]
        };
        let exact_strain = |x: [f64; 2]| -> [[f64; 2]; 2] {
            let (xx, y) = (x[0], x[1]);
            let exx = -(P / (E * i)) * (L - xx) * y;
            let eyy = NU * (P / (E * i)) * (L - xx) * y;
            let exy = 0.5 * (P / (E * i)) * (1.0 + NU) * (H * H / 4.0 - y * y);
            [[exx, exy], [exy, eyy]]
        };
        let exact_traction = |x: [f64; 2]| -> [f64; 2] {
            let y = x[1];
            [0.0, P / (2.0 * i) * (H * H / 4.0 - y * y)]
        };

        let mesh = read_gmsh::<2>("examples/msh/cantilever.msh");
        let units: Vec<Quad4<2>> = mesh.elements("Ω");
        let left: Vec<Seg2<2>> = mesh.elements("Γu");
        let right: Vec<Seg2<2>> = mesh.elements("Γt");

        let max_id = mesh.nodes.keys().copied().max().unwrap();
        let n_dofs = (max_id + 1) * 2;
        let lambda = plane_stress_lambda(E, NU);
        let (_, mu) = lame_constants(E, NU);

        let mut triplets = stiffness(&units, lambda, mu);
        let (pt, mut doublets) = displacement_penalty(&left, 1e7, exact);
        triplets.extend(pt);
        doublets.extend(traction(&right, exact_traction));

        let k =
            SparseColMat::<usize, f64>::try_new_from_triplets(n_dofs, n_dofs, &triplets).unwrap();
        let f = try_new_from_doublet(n_dofs, &doublets);
        let d = k.sp_cholesky(Side::Lower).unwrap().solve(&f);

        let l2 = l2_error(&units, |id, p| d[id * 2 + p], exact);
        let energy = energy_error(&units, |id, p| d[id * 2 + p], exact_strain, lambda, mu);

        assert!(l2 < 5e-2, "cantilever L2 error too large: {l2:e}");
        assert!(energy < 2e-1, "cantilever energy error too large: {energy:e}");
    }
}
