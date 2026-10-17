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

use crate::element::{Element, FullIntegration, IntegrationScheme, ReducedIntegration};
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

/// Bulk modulus in `D` spatial dimensions: `K = λ + 2μ/D`.
///
/// For `D = 3` this is the standard `K = λ + 2μ/3`; for 2D plane strain it is
/// `K = λ + μ`. It is the coefficient of the **dilatational** part of the
/// isotropic bilinear form (see [`dilatational_stiffness`]); together with the
/// shear modulus `μ` it gives the volumetric/deviatoric split
/// `C = K I⊗I + 2μ (I − I⊗I/D)`.
///
pub fn bulk_modulus<const D: usize>(lambda: f64, mu: f64) -> f64 {
    lambda + 2.0 * mu / D as f64
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
/// This operator equals the sum of its volumetric and deviatoric parts,
/// `stiffness(v, λ, μ) == dilatational_stiffness(v, K) + deviatoric_stiffness(v, μ)`
/// with `K =` [`bulk_modulus::<D>`](bulk_modulus).
///
/// Diagonal node blocks `(i, i)` emit the lower triangle of their `D×D` block
/// (`b ≤ a`), while off-diagonal blocks `(i, j)`, `j < i`, are emitted in full.
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
            for i in 0..N {
                // Diagonal node block (i, i): lower triangle (b ≤ a).
                let grad_dot = (0..D).map(|k| dshape[i][k] * dshape[i][k]).sum::<f64>();
                for a in 0..D {
                    for b in 0..=a {
                        let delta = if a == b { 1.0 } else { 0.0 };
                        let value = lambda * dshape[i][a] * dshape[i][b]
                            + mu * (delta * grad_dot + dshape[i][b] * dshape[i][a]);
                        triplets.push(Triplet::new(
                            id[i] * D + a,
                            id[i] * D + b,
                            value * jacobe * wg,
                        ));
                    }
                }
                // Off-diagonal node blocks (i, j), j < i: full D×D block.
                for j in 0..i {
                    let grad_dot = (0..D).map(|k| dshape[i][k] * dshape[j][k]).sum::<f64>();
                    for a in 0..D {
                        for b in 0..D {
                            let delta = if a == b { 1.0 } else { 0.0 };
                            let value = lambda * dshape[i][a] * dshape[j][b]
                                + mu * (delta * grad_dot + dshape[i][b] * dshape[j][a]);
                            let g_a = id[i] * D + a;
                            let g_b = id[j] * D + b;
                            triplets.push(Triplet::new(
                                g_a.max(g_b),
                                g_a.min(g_b),
                                value * jacobe * wg,
                            ));
                        }
                    }
                }
            }
        }
    }
    triplets
}

/// Assembles the **dilatational (volumetric) stiffness** matrix and returns
/// sparse triplets for the lower triangle (including the diagonal).
///
/// The isotropic bilinear form splits additively into a volumetric and a
/// deviatoric part, `a(u, v) = a_dil(u, v) + a_dev(u, v)`, where
///
/// `a_dil(u, v) = ∫ K (∇·u)(∇·v) dΩ`
///
/// is the dilatational part, with `K = λ + 2μ/D` the bulk modulus
/// (see [`bulk_modulus`]). In components, with `u = Nᵢ eₐ` and `v = Nⱼ e_b`,
///
/// `K_{(i,a),(j,b)} = ∫ K (∂ₐNᵢ)(∂_b Nⱼ) dΩ`,
///
/// so this operator only involves the divergence of the shape functions.
///
/// Unlike [`stiffness`] and [`deviatoric_stiffness`], which use full
/// integration, this operator deliberately uses **reduced integration**
/// ([`ReducedIntegration<P, G>`]): selective reduced integration of the
/// volumetric term is the standard remedy for **volumetric locking** in the
/// nearly incompressible limit (`ν → 1/2`, `K → ∞`).
///
/// For elements with a constant gradient (e.g. the three-node triangle) one
/// reduced point already integrates the dilatational term exactly, so
/// `stiffness(v, λ, μ) == dilatational_stiffness(v, K) + deviatoric_stiffness(v, μ)`
/// still holds; for bilinear quadrilaterals the reduced rule differs from full
/// integration by design.
///
/// `bulk` is the bulk modulus `K`. `v` is a set of elements (implementing
/// [`Variable<D, 1, P, N>`]); each element must implement
/// [`ReducedIntegration<P, G>`].
///
/// Diagonal node blocks `(i, i)` emit the lower triangle of their `D×D` block
/// (`b ≤ a`), while off-diagonal blocks `(i, j)`, `j < i`, are emitted in full.
pub fn dilatational_stiffness<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
    bulk: f64,
) -> Triplets
where
    T: Variable<D, 1, P, N>,
    T::Item: Element<D, P, N> + ReducedIntegration<P, G>,
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
            for i in 0..N {
                // Diagonal node block (i, i): lower triangle (b ≤ a).
                for a in 0..D {
                    for b in 0..=a {
                        // K_{(i,a),(i,b)} += K (∂ₐNᵢ)(∂_b Nᵢ) |J| w.
                        let value = bulk * dshape[i][a] * dshape[i][b];
                        triplets.push(Triplet::new(
                            id[i] * D + a,
                            id[i] * D + b,
                            value * jacobe * wg,
                        ));
                    }
                }
                // Off-diagonal node blocks (i, j), j < i: full D×D block.
                for j in 0..i {
                    for a in 0..D {
                        for b in 0..D {
                            // K_{(i,a),(j,b)} += K (∂ₐNᵢ)(∂_b Nⱼ) |J| w.
                            let value = bulk * dshape[i][a] * dshape[j][b];
                            let g_a = id[i] * D + a;
                            let g_b = id[j] * D + b;
                            triplets.push(Triplet::new(
                                g_a.max(g_b),
                                g_a.min(g_b),
                                value * jacobe * wg,
                            ));
                        }
                    }
                }
            }
        }
    }
    triplets
}

/// Assembles the **deviatoric stiffness** matrix and returns sparse triplets for
/// the lower triangle (including the diagonal).
///
/// This is the second part of the volumetric/deviatoric split,
/// `a_dev(u, v) = ∫ 2μ ε'(u) : ε'(v) dΩ`, where
/// `ε' = ε − (1/D)(tr ε) I` is the deviatoric strain and `μ` the shear
/// modulus. In components, with `u = Nᵢ eₐ` and `v = Nⱼ e_b` (summation over
/// `k`),
///
/// `K_{(i,a),(j,b)} = ∫ { μ [ δ_ab (∂ₖNᵢ)(∂ₖNⱼ) + (∂_b Nᵢ)(∂ₐNⱼ) ]
///                        − (2μ/D) (∂ₐNᵢ)(∂_b Nⱼ) } dΩ`.
///
/// It is invariant to the volumetric (dilatational) part of the displacement
/// and, together with [`dilatational_stiffness`], reproduces [`stiffness`]:
/// `stiffness(v, λ, μ) == dilatational_stiffness(v, K) + deviatoric_stiffness(v, μ)`.
///
/// `mu` is the shear modulus. `v` is a set of elements (implementing
/// [`Variable<D, 1, P, N>`]); each element must implement
/// [`FullIntegration<P, G>`].
///
/// Diagonal node blocks `(i, i)` emit the lower triangle of their `D×D` block
/// (`b ≤ a`), while off-diagonal blocks `(i, j)`, `j < i`, are emitted in full.
pub fn deviatoric_stiffness<T, const D: usize, const P: usize, const N: usize, const G: usize>(
    v: &T,
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
            for i in 0..N {
                // Diagonal node block (i, i): lower triangle (b ≤ a).
                let grad_dot = (0..D).map(|k| dshape[i][k] * dshape[i][k]).sum::<f64>();
                for a in 0..D {
                    for b in 0..=a {
                        let delta = if a == b { 1.0 } else { 0.0 };
                        let value = mu * (delta * grad_dot + dshape[i][b] * dshape[i][a])
                            - (2.0 * mu / D as f64) * dshape[i][a] * dshape[i][b];
                        triplets.push(Triplet::new(
                            id[i] * D + a,
                            id[i] * D + b,
                            value * jacobe * wg,
                        ));
                    }
                }
                // Off-diagonal node blocks (i, j), j < i: full D×D block.
                for j in 0..i {
                    let grad_dot = (0..D).map(|k| dshape[i][k] * dshape[j][k]).sum::<f64>();
                    for a in 0..D {
                        for b in 0..D {
                            let delta = if a == b { 1.0 } else { 0.0 };
                            let value = mu * (delta * grad_dot + dshape[i][b] * dshape[j][a])
                                - (2.0 * mu / D as f64) * dshape[i][a] * dshape[j][b];
                            let g_a = id[i] * D + a;
                            let g_b = id[j] * D + b;
                            triplets.push(Triplet::new(
                                g_a.max(g_b),
                                g_a.min(g_b),
                                value * jacobe * wg,
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
            let un: [[f64; D]; N] = std::array::from_fn(|a| std::array::from_fn(|i| u_h(id[a], i)));
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
    use crate::element::quad4::Quad4;
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
    fn stiffness_entries_match_reference() {
        // Unit triangle, λ = μ = 1. Node 0 has ∇N₀ = (-1, -1), so the element
        // (0,0) entry is [λ(-1)(-1) + μ(1·2 + (-1)(-1))]·area = 4·(1/2) = 2,
        // and the (0,1) entry is [λ(-1)(-1) + μ(0 + (-1)(-1))]·area = 2·(1/2) = 1.
        let tris = unit_triangle();
        let triplets = stiffness(&tris, 1.0, 1.0);
        let (mut k00, mut k10) = (0.0, 0.0);
        for t in &triplets {
            match (t.row, t.col) {
                (0, 0) => k00 += t.val,
                (1, 0) => k10 += t.val,
                _ => {}
            }
        }
        assert!(close(k00, 2.0));
        assert!(close(k10, 1.0));
    }

    /// Evaluates the quadratic form `dᵀ K d` from the lower-triangle triplets.
    fn quadratic(triplets: &Triplets, d: &[f64]) -> f64 {
        triplets
            .iter()
            .map(|t| {
                let weight = if t.row == t.col { 1.0 } else { 2.0 };
                weight * t.val * d[t.row] * d[t.col]
            })
            .sum()
    }

    #[test]
    fn dilatational_plus_deviatoric_equals_isotropic() {
        // The volumetric/deviatoric split must reproduce the full stiffness for
        // K = λ + 2μ/D. The dilatational part uses reduced integration, but on a
        // constant-gradient triangle one point is exact, so the identity holds.
        // Compare the quadratic forms for a few random-ish nodal displacement
        // vectors, which exercises every matrix entry.
        let tris = unit_triangle();
        let (lambda, mu) = (1.0, 2.0);
        let bulk = bulk_modulus::<2>(lambda, mu);

        let iso = stiffness(&tris, lambda, mu);
        let mut split = dilatational_stiffness(&tris, bulk);
        split.extend(deviatoric_stiffness(&tris, mu));

        for d in [
            [1.0, 2.0, 3.0, -1.0, 0.5, 2.0],
            [-2.0, 0.3, 1.5, 4.0, -0.7, 1.1],
            [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
        ] {
            let expected = quadratic(&iso, &d);
            let actual = quadratic(&split, &d);
            assert!(
                close(actual, expected),
                "split quadratic form {actual} != isotropic {expected}"
            );
        }
    }

    #[test]
    fn dilatational_stiffness_uses_reduced_integration() {
        // Unit-square Quad4 with nodes (0,0), (1,0), (1,1), (0,1). The reduced
        // 1-point rule sits at the reference centroid with weight 4 and
        // |J| = 1/4, so |J| w = 1; the centroid physical gradients are ±0.5 and
        // the dilatational stiffness is K gᵢ⊗gⱼ. This differs from full 2×2
        // integration (e.g. the (N₀,x)-(N₁,x) entry is -K/4 here, while the
        // exact integral is -K/3), which proves the reduced rule is in use.
        let quads = vec![Quad4::new(
            Node::new(0, [0.0, 0.0]),
            Node::new(1, [1.0, 0.0]),
            Node::new(2, [1.0, 1.0]),
            Node::new(3, [0.0, 1.0]),
        )];
        let bulk = 3.0;
        let triplets = dilatational_stiffness(&quads, bulk);

        let g = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]];
        let mut dense = [[0.0; 8]; 8];
        for t in &triplets {
            dense[t.row][t.col] += t.val;
        }
        for i in 0..4 {
            for j in 0..4 {
                for a in 0..2 {
                    for b in 0..2 {
                        let row = i * 2 + a;
                        let col = j * 2 + b;
                        if row >= col {
                            let expected = bulk * g[i][a] * g[j][b];
                            assert!(
                                close(dense[row][col], expected),
                                "entry ({row},{col}): {} != {expected}",
                                dense[row][col]
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn deviatoric_stiffness_ignores_pure_dilation() {
        // u = α(x, y) is a pure dilation: ε = α I, tr ε = 2α, ε' = 0.
        // Therefore a_dev(u, u) = 0 and a_dil(u, u) = ∫ K (2α)² dΩ = 2 K α²
        // (the reference triangle has area 1/2).
        let tris = unit_triangle();
        let (lambda, mu) = (1.0, 2.0);
        let bulk = bulk_modulus::<2>(lambda, mu);
        let alpha = 0.3;
        // Node-major DOFs: (u,v) at (0,0), (1,0), (0,1).
        let d = [0.0, 0.0, alpha, 0.0, 0.0, alpha];

        let dev = quadratic(&deviatoric_stiffness(&tris, mu), &d);
        assert!(
            dev.abs() < 1e-12,
            "deviatoric energy should vanish: {dev:e}"
        );

        let dil = quadratic(&dilatational_stiffness(&tris, bulk), &d);
        assert!(close(dil, 2.0 * bulk * alpha * alpha));
    }

    #[test]
    fn dilatational_stiffness_ignores_pure_shear() {
        // u = γ(y, 0) is a pure shear: ε_xx = ε_yy = 0, ε_xy = γ/2, tr ε = 0.
        // Therefore a_dil(u, u) = 0 and a_dev(u, u) = ∫ 2μ (γ²/2) dΩ = μ γ²/2.
        let tris = unit_triangle();
        let (lambda, mu) = (1.0, 2.0);
        let bulk = bulk_modulus::<2>(lambda, mu);
        let gamma = 0.4;
        // u = (γ y, 0): only the top-right node (0,1) has u_x = γ.
        let d = [0.0, 0.0, 0.0, 0.0, gamma, 0.0];

        let dil = quadratic(&dilatational_stiffness(&tris, bulk), &d);
        assert!(
            dil.abs() < 1e-12,
            "dilatational energy should vanish: {dil:e}"
        );

        let dev = quadratic(&deviatoric_stiffness(&tris, mu), &d);
        assert!(close(dev, mu * gamma * gamma / 2.0));
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
        assert!(
            energy < 2e-1,
            "cantilever energy error too large: {energy:e}"
        );
    }
}
