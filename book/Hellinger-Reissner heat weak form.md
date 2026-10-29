# Hellinger–Reissner weak form of the heat conduction problem

> Fourth note of the series. It uses the notation of
> [[Crouzeix-Raviart triangular element]], [[Variational consistency of the CR element]],
> [[Legendre moment discretization and the consistency remainder]] and
> [[Raviart-Thomas element]], and it is the variational *parent* of the last two:
> the two-field (Hellinger–Reissner) principle produces the CR method when the flux is
> element-local, the RT0 mixed method when the flux is $H(\mathrm{div})$-conforming, and —
> after hybridization — both.
>
> Verified numerically on the GMSH meshes of `msh/`:
> `cargo run --bin hellinger_reissner_heat` (RT0/P0 saddle point, $O(h)$ in both fields) and
> `cargo run --bin complementary_multiplier_heat` (the multiplier form of §4.3, where the
> multiplier is the temperature), and `cargo run --bin rt0_cr_mixed` (the consistency *and*
> stability check of (B1)–(B2) for the pairing CR × RT0, §7.1) and `cargo run --bin rt1_cr_mixed`
> (the same for CR × RT1, §7.2).
>
> **Citations.** Keys such as `[Re50]`, `[BF91]`, `[Bre74]`, `[For77]` are defined in
> [[Raviart-Thomas element]] §10; the per-result mapping for this note is §11 at the end.

---

## 1. Strong form, signs, and the two first-order equations

Steady heat conduction in $\Omega\subset\mathbb R^2$:

$$
\boxed{\;\boldsymbol q+k\nabla T=\boldsymbol 0\;}
\quad\text{(Fourier law)},\qquad
\boxed{\;\nabla\cdot\boldsymbol q+f=0\;}
\quad\text{(energy balance)},
$$

$$
T=T_D\ \text{on }\Gamma_D,\qquad
\boldsymbol q\cdot\boldsymbol n=q_n\ \text{on }\Gamma_N,\qquad
\boldsymbol q\cdot\boldsymbol n=h\,(T-T_\infty)\ \text{on }\Gamma_R .
$$

The sign convention matches `ansatz::operation::heat`: with $\boldsymbol n$ the **outward** normal
and $k>0$, `q_n = q·n = −k ∂T/∂n` (so a positive `q_n` is heat leaving the body); `ansatz`
integrates `−q_n` in the load vector for exactly this reason. Eliminating $\boldsymbol q$ gives
$k\Delta T=f$, i.e. with $f>0$ the body is *heated at the boundary* and has a maximum inside — note
the sign, which is easy to flip.

Everything below is the two-field version of the first note's model problem (there $k=1$ and
$\boldsymbol q=-\nabla u$).

---

## 2. Four variational principles for one problem

| principle | independent fields | functional | type |
|---|---|---|---|
| Dirichlet (energy / primal) | $T\in H^1_{\Gamma_D}$ | $\Pi(T)=\frac12\int_\Omega k\lvert\nabla T\rvert^2-\int_\Omega fT-\int_{\Gamma_N}q_nT$ | minimum |
| complementary (dual) | $\boldsymbol q$, with $\nabla\cdot\boldsymbol q+f=0$ imposed | $\Pi^{*}(\boldsymbol q)=-\frac{1}{2k}\int_\Omega\lvert\boldsymbol q\rvert^2+\int_{\Gamma_D}T_D\,\boldsymbol q\cdot\boldsymbol n$ | maximum |
| **Hellinger–Reissner (H–R)** | $T,\ \boldsymbol q$ **independent** | $\Pi_{HR}(T,\boldsymbol q)=\int_\Omega\bigl(\boldsymbol q\cdot\nabla T-\frac{1}{2k}\lvert\boldsymbol q\rvert^2-fT\bigr)-\int_{\Gamma_N}q_nT$ | **saddle** |
| complementary with **multiplier** constraints (§4.3) | $\boldsymbol q,\ T,\ \mu$ | $\mathcal L=-\frac{1}{2k}\int\lvert\boldsymbol q\rvert^2+\int_{\Gamma_D}T_D\boldsymbol q\cdot\boldsymbol n-\int T(\nabla\cdot\boldsymbol q+f)-\int_{\Gamma_N}\mu(\boldsymbol q\cdot\boldsymbol n-q_n)$, with $\mu=T|_{\Gamma_N}$ | **saddle** |
| Hu–Washizu | $T,\ \boldsymbol\gamma,\ \boldsymbol q$ | relaxes both Fourier law and the gradient relation | saddle |

The name is the heat analogue of the elasticity two-field principle in stress and displacement
(Reissner 1950; Hellinger 1914). It is **not** a minimum principle, and that single fact dictates
almost everything about its discretization (§5, §7): the discrete operator is indefinite and needs
an inf-sup condition, and the linear solver must factor an indefinite matrix (`sp_lu`, exactly as
noted in `ansatz`'s heat module for its Lagrange-multiplier variant).

---

## 3. Derivation by convex duality (Legendre transform)

Write the primal energy density as a function of the gradient,

$$
g(\boldsymbol\gamma)=\frac{k}{2}\lvert\boldsymbol\gamma\rvert^{2},
\qquad
g^{*}(\boldsymbol q)=\sup_{\boldsymbol\gamma}\bigl[\boldsymbol q\cdot\boldsymbol\gamma-g(\boldsymbol\gamma)\bigr]
=\frac{1}{2k}\lvert\boldsymbol q\rvert^{2}
\quad\Longleftrightarrow\quad
\boldsymbol q=k\boldsymbol\gamma=\arg\max ,
$$

i.e. $g$ and $g^{*}$ are Legendre transforms of each other, and for every $\boldsymbol\gamma$

$$
g(\boldsymbol\gamma)=\sup_{\boldsymbol q}\bigl[\boldsymbol q\cdot\boldsymbol\gamma-g^{*}(\boldsymbol q)\bigr]
\qquad(\text{the sup is attained at }\boldsymbol q=k\boldsymbol\gamma).
$$

Insert $\boldsymbol\gamma=\nabla T$ into $\Pi$:

$$
\Pi(T)=\inf_{\boldsymbol q}\ \underbrace{\int_\Omega\Bigl[\boldsymbol q\cdot\nabla T-\frac{1}{2k}\lvert\boldsymbol q\rvert^{2}-fT\Bigr]d\Omega
-\int_{\Gamma_N}q_nT\,ds}_{=:\ \Pi_{HR}(T,\boldsymbol q)} ,
$$

so the primal minimisation is the $\inf_T\sup_{\boldsymbol q}$ of the H–R functional:

$$
\boxed{\;
\Pi_{HR}(T,\boldsymbol q)=\int_\Omega\Bigl(\boldsymbol q\cdot\nabla T-\frac{1}{2k}\lvert\boldsymbol q\rvert^{2}-fT\Bigr)d\Omega
-\int_{\Gamma_N}q_n\,T\,ds
\;}
\qquad
T\big|_{\Gamma_D}=T_D .
$$

(The Neumann datum enters as a boundary term; the Dirichlet datum is an *essential* constraint of the
functional — which is exactly opposite to what happens after integrating by parts, see §4.2.)

**Stationarity.**

- $\delta\boldsymbol q$: $\displaystyle\int_\Omega\Bigl(\nabla T-\frac{1}{k}\boldsymbol q\Bigr)\cdot\delta\boldsymbol q\,d\Omega=0$
  $\ \Longrightarrow\ $ Fourier law $\boldsymbol q=-k\nabla T$ (since $\delta\boldsymbol q$ is free in $L^2$);
- $\delta T$ (with $\delta T=0$ on $\Gamma_D$): $\displaystyle-\int_\Omega(\nabla\cdot\boldsymbol q+f)\,\delta T\,d\Omega
+\int_{\Gamma_N}\bigl(\boldsymbol q\cdot\boldsymbol n-q_n\bigr)\delta T\,ds=0$
  $\ \Longrightarrow\ $ balance plus the Neumann condition.

So the Euler–Lagrange equations of $\Pi_{HR}$ are precisely the strong form; the functional is
*consistent by construction* — there is no "consistency crime" at the continuous level, only at the
discrete level (§6).

**Equivalence with the primal.** Because $g$ is convex and quadratic, the saddle point exists and is
unique; eliminating $\boldsymbol q$ from the stationarity equations returns $\Pi$ (and hence the
primal weak form of the first note). Conversely eliminating $T$ returns the dual principle of §2.

---

## 4. Two equivalent weak forms, and why the second one is the one to discretize

### 4.1 Form A ("Reissner form", keeps $\nabla T$)

Multiply the two first-order equations by test functions $\boldsymbol r$, $w$ and integrate with **no**
integration by parts:

$$
\text{(A1)}\quad\int_\Omega\Bigl(\frac{1}{k}\boldsymbol q\cdot\boldsymbol r+\nabla T\cdot\boldsymbol r\Bigr)d\Omega=0 ,
\qquad
\text{(A2)}\quad\int_\Omega\bigl(\nabla\cdot\boldsymbol q+f\bigr)w\,d\Omega=0 .
$$

This needs $T\in H^1$ (it contains $\nabla T$), so the temperature must be **globally continuous** —
and that is fatal for low-order pairs: if $T_h$ is piecewise constant then $\nabla T_h\equiv0$ element
by element and (A1) forces $\boldsymbol q_h\equiv0$. Equal-order continuous pairs avoid the trivial
kernel but are still inf-sup unstable (checkerboard). Form A is therefore only useful with
stabilization.

### 4.2 Form B ("mixed form", $T$ as the multiplier of the balance)

Integrate (A1) by parts **element-wise** (legitimate because $\boldsymbol r\in H(\mathrm{div})$ has a
continuous normal trace, so interior contributions cancel):

$$
\text{(B1)}\quad\int_\Omega\frac{1}{k}\boldsymbol q_h\cdot\boldsymbol r\,d\Omega
-\int_\Omega T_h\,\nabla\cdot\boldsymbol r\,d\Omega
=-\int_{\Gamma_D}T_D\,\boldsymbol r\cdot\boldsymbol n\,ds ,
$$

$$
\text{(B2)}\quad-\int_\Omega\nabla\cdot\boldsymbol q_h\,w\,d\Omega=\int_\Omega f\,w\,d\Omega ,
$$

for all $(\boldsymbol r,w)\in H(\mathrm{div})\times L^2$ (with $\boldsymbol r\cdot\boldsymbol n=q_n$ on
$\Gamma_N$). Now $T$ is an **$L^2$ multiplier** — it may jump — and the flux carries the conformity.
Sentences worth remembering:

- **the Dirichlet datum became *natural*** (it sits in the right-hand side of (B1), no constraint on $T$);
- **the Neumann datum became *essential*** (it constrains the flux space);
- the two equations are the Euler equations of the *dual* member of the H–R family: eliminating
  $\boldsymbol q$ from (B1)–(B2) returns the primal form, and eliminating $T$ returns the
  complementary principle.

This is the form used by the RT mixed method and hence by the last note of the series.

---

### 4.3 The same principle read from the complementary side: constraints relaxed by multipliers

The complementary principle of §2 maximizes the stored (dual) energy over **equilibrated** fluxes,
i.e. over fields that satisfy the balance equation and the Neumann data *exactly*:

$$
\Pi^{*}(\boldsymbol q)=-\frac{1}{2k}\int_\Omega\lvert\boldsymbol q\rvert^{2}d\Omega
+\int_{\Gamma_D}T_D\,\boldsymbol q\cdot\boldsymbol n\,ds
\qquad\text{s.t.}\quad \nabla\cdot\boldsymbol q+f=0,\quad \boldsymbol q\cdot\boldsymbol n=q_n\ \text{on }\Gamma_N .
$$

Its admissible set is impossible to discretize directly for low-order fluxes (an element-constant
$\boldsymbol q_h$ cannot satisfy both conditions). The classical remedy is the **Lagrangian (multiplier)
form**: drop the constraints and let them be enforced by multipliers, one for the balance equation and
one for the Neumann datum:

$$
\boxed{\;
\mathcal L(\boldsymbol q,T,\mu)=-\frac{1}{2k}\int_\Omega\lvert\boldsymbol q\rvert^{2}d\Omega
+\int_{\Gamma_D}T_D\,\boldsymbol q\cdot\boldsymbol n\,ds
-\int_\Omega T\bigl(\nabla\cdot\boldsymbol q+f\bigr)d\Omega
-\int_{\Gamma_N}\mu\bigl(\boldsymbol q\cdot\boldsymbol n-q_n\bigr)ds
\;}
$$

with $\boldsymbol q$ now in $L^{2}(\Omega)^2$ — **nothing constrains the flux**.

**Stationarity.**
$\delta T$ returns the balance equation, $\delta\mu$ returns the Neumann datum, and $\delta\boldsymbol q$
gives, after one integration by parts on $\Omega$,

$$
\int_\Omega\Bigl(\nabla T-\frac{1}{k}\boldsymbol q\Bigr)\cdot\delta\boldsymbol q\,d\Omega
+\int_{\Gamma_D}\bigl(T_D-T\bigr)\,\delta\boldsymbol q\cdot\boldsymbol n\,ds
+\int_{\Gamma_N}\bigl(T-\mu\bigr)\,\delta\boldsymbol q\cdot\boldsymbol n\,ds=0 ,
$$

so, since $\delta\boldsymbol q$ is free in $L^2$ and its normal trace is free on each boundary part,

$$
\boxed{\;\boldsymbol q=-k\nabla T\ \text{in }\Omega,\qquad
T=T_D\ \text{on }\Gamma_D,\qquad \mu=T\big|_{\Gamma_N}\;}
$$

**The multiplier of the Neumann datum is the temperature.** Two comments:

1. This is not an extra hypothesis — it is *derived*: the boundary multiplier equals the trace of the
   bulk multiplier, so a single temperature field serves both. The phrasing "the multiplier is the
   temperature" is the statement of this identification.
2. The Dirichlet datum, which was essential in (A) and natural in (B), appears here as the *condition
   that makes the multiplier field admissible* ($T=T_D$ on $\Gamma_D$): all three readings describe the
   same saddle point.

Eliminating $\mu$ (substituting $\mu=T$ and moving the boundary term to the left) returns exactly the
H–R weak form (B1)–(B2) of §4.2, so this is a re-reading of the same principle, not a new method — but
a very useful one, because it says *what the temperature means*: it is the reaction to the balance
equation, and it is the natural unknown with which to relax Neumann data.

### 4.4 Discrete form of the multiplier principle

With $\boldsymbol q_h\in\mathcal{RT}_0$, $T_h\in\mathcal{P}_0$ and a multiplier space $M_h$ on $\Gamma_N$:

$$
\text{(E1)}\ \int_\Omega\frac{1}{k}\boldsymbol q_h\cdot\boldsymbol r-\sum_T\int_T T_h\nabla\cdot\boldsymbol r
=-\int_{\Gamma_D}T_D\,\boldsymbol r\cdot\boldsymbol n-\int_{\Gamma_N}\mu_h\,\boldsymbol r\cdot\boldsymbol n ,
\qquad\forall\boldsymbol r\in\mathcal{RT}_0,
$$

$$
\text{(E2)}\ -\sum_T\int_T\nabla\cdot\boldsymbol q_h\,w=\int_\Omega f w,\ \ \forall w\in\mathcal{P}_0,
\qquad
\text{(E3)}\ \int_{\Gamma_N}\bigl(\boldsymbol q_h\cdot\boldsymbol n-q_n\bigr)\zeta\,ds=0,\ \ \forall\zeta\in M_h .
$$

Three structural facts (all verified in §9.2):

**(a) The multiplier of RT0 is algebraic, not an extra constraint.** On $\Gamma_e$ the normal trace of
$\mathcal{RT}_0$ is a **single constant**, so
$\int_{\Gamma_e}(\boldsymbol q_h\cdot\boldsymbol n-q_n)\zeta\,ds$ for $\zeta=1$ reads
$\Phi_e=\int_{\Gamma_e}q_n$ — *the strong imposition itself*. Hence with $M_h=\{$constants on each
edge$\}$ the fields $(\boldsymbol q_h,T_h)$ obtained with (E3) coincide **identically** with those of the
strong-imposition method, and $\mu_e$ is merely the reaction read off from (E1):
$\mu_e=-\frac{1}{\sigma_e}\bigl[\ldots\bigr]$ on the edge's own row.

**(b) $\mu_h$ converges to the temperature, at second order.** Since $\mu=T|_{\Gamma_N}$ continuously and
the multiplier is an edge constant, $\mu_e$ approximates the *edge average* of $T$; for smooth data the
observed rate is $O(h^{2})$ (measured: 3.2\,e-2 → 7.8\,e-3 → 1.9\,e-3 → 3.6\,e-4, rate ≈ 2), while
$\boldsymbol q_h$ and $T_h$ keep their $O(h)$ rates.

**(c) An over-rich multiplier is *not* solvable — the moment obstruction.** If $M_h$ contains $\xi$
(i.e. a linear multiplier), then (E3) with $\zeta=\xi$ reads

$$
\int_{\Gamma_e}\boldsymbol q_h\cdot\boldsymbol n\,\xi\,ds=\int_{\Gamma_e}q_n\xi\,ds
\iff
\underbrace{\Phi_e}_{\text{constant}}\int_{\Gamma_e}\xi\,ds=0=\int_{\Gamma_e}q_n\xi\,ds ,
$$

so the constraint contains **no unknown** and the discrete problem is solvable only if the data satisfy
the compatibility condition $\int_{\Gamma_e}q_n\xi\,ds=0$ on every Neumann edge; otherwise the system is
singular. Verified: with $q_n=\pi\sin(\pi x)$ (whose first moment is $\ne0$, $\max_e\lvert\int q_n\xi\rvert
=2.5\times10^{-2}$) the linear-multiplier system is indeed unsolvable. This is the same phenomenon as
the design rule of [[Legendre moment discretization and the consistency remainder]] §10: *the multiplier
may not be richer than the trace space of the flux*; to use a linear multiplier one needs an enriched
flux (RT1), whose normal trace carries the $\mu_{e,1}$ moment.

**(d) Pure Neumann problems.** If $\Gamma_D=\emptyset$, (E2) is solvable only for compatible data
($\int_\Omega f+\int_{\Gamma_N}q_n=0$) and $T_h$ (hence $\mu_h$) is determined only up to a constant —
the multiplier inherits the null space of the continuous problem, and a normalization has to be added
(e.g. $\int_\Omega T_h=0$). This is the natural place to see why the reaction $\mu$ and the additive
constant are the same object.

**Why use it then?** Not to replace the strong imposition (for RT0 the two coincide), but because

- the temperature on $\Gamma_N$ becomes an **output** (the reaction) instead of being unavailable;
- data that do not lie in the flux trace space (measured fluxes, data given on a finer space, Robin
  conditions linearised about a measured $T_\infty$) can be admitted variationally, with the compatibility
  condition of (c) as the price;
- it is the clean statement of the duality "temperature $=$ multiplier of the balance", which is what
  makes hypercircle-type a-posteriori bounds work (§10).

## 5. Well-posedness of the two-field problem

With $V=H_{\Gamma_N}(\mathrm{div}):=\{\boldsymbol r\in H(\mathrm{div}):\boldsymbol r\cdot\boldsymbol n=0
\text{ on }\Gamma_N\}$, $Q=L^2(\Omega)$, $a(\boldsymbol q,\boldsymbol r)=\int\frac1k\boldsymbol q\cdot\boldsymbol r$,
$b(\boldsymbol r,w)=-\int w\nabla\cdot\boldsymbol r$, $\ell(\boldsymbol r)=-\int_{\Gamma_D}T_D\boldsymbol r\cdot\boldsymbol n$,
$m(w)=\int fw$, the problem reads

$$
a(\boldsymbol q,\boldsymbol r)+b(\boldsymbol r,T)=\ell(\boldsymbol r)\ \ \forall\boldsymbol r,
\qquad
b(\boldsymbol q,w)=m(w)\ \ \forall w ,
$$

and Brezzi's theory asks for exactly two things (note $a$ is symmetric positive definite, so this is a
saddle point of a *constrained* minimisation of $\int\lvert\boldsymbol q\rvert^2/2k$):

1. **kernel coercivity**: $a(\boldsymbol r,\boldsymbol r)\ge\alpha\lVert\boldsymbol r\rVert_{H(\mathrm{div})}^{2}$
   on $\ker B=\{\boldsymbol r:\nabla\cdot\boldsymbol r=0,\ \boldsymbol r\cdot\boldsymbol n=0\text{ on }\Gamma_N\}$
   — true with $\alpha=1/k$;
2. **inf-sup (LBB)**: $\sup_{\boldsymbol r}\frac{b(\boldsymbol r,w)}{\lVert\boldsymbol r\rVert_{H(\mathrm{div})}}\ge\beta\lVert w\rVert_{L^2}$
   — the classical $\mathrm{div}$-surjectivity of $H(\mathrm{div})$ onto $L^2$.

The inf-sup condition is *not* inherited by subspaces: it must be checked for the chosen discrete pair.
$(\mathcal{RT}_0,\mathcal{P}_0)$ passes, equal-order $(\mathcal{P}_1,\mathcal{P}_1)$ fails. The
stability constant $\beta$ is the reason the mixed method's error bounds for $\boldsymbol q$ and $T$ are
coupled (§7), and it is why the discrete system is indefinite and must be factored by LU
(`SparseColMat::sp_lu`), never Cholesky.

---

## 6. Consistency of the H–R weak form

Because the H–R functional *is* the strong form (§3), the exact solution $(\boldsymbol q,T)$ satisfies
(B1)–(B2) with every test function in $H(\mathrm{div})\times L^2$: **consistency holds exactly, with no
residual**. Compare with the non-conforming CR method of [[Variational consistency of the CR element]],
whose residual is the jump functional $\sum_e\int_e[\![v_h]\!]\partial_{\boldsymbol n}u$ of size
$O(h)$. The price moves elsewhere:

| method | what is violated | size of the defect |
|---|---|---|
| CR (primal, non-conforming) | continuity of $T_h$ across edges | $R(v_h)=O(h)$, needs Strang's lemma |
| H–R / RT0 mixed | nothing at the continuous level; the discrete space must satisfy inf-sup | $O(h)$ **approximation** error of $\boldsymbol q$ and $T$ |
| H–R in moment language | the flux moments $\mu_{e,n}$, $n\ge1$, are not represented by RT0 | $O(h)$ |

In the moment language of [[Legendre moment discretization and the consistency remainder]] this last
line is the right way to read the H–R discretization: the relaxed Fourier law is tested only against
fluxes whose edge data are the moments the discrete space can carry, so the missing moment
$\mu_{e,1}=\frac{\ell_e}{6}D^{2}T[\boldsymbol t_e,\boldsymbol n_e]$ is precisely the $O(h)$ term.

---

## 7. Discretization: one family, three methods

Take finite-dimensional $V_h\subset H(\mathrm{div})$, $Q_h\subset L^2$ and use Form B:

$$
\begin{bmatrix}\boldsymbol A & \boldsymbol C^{T}\\ \boldsymbol C & \boldsymbol 0\end{bmatrix}
\begin{bmatrix}\boldsymbol q\\ \boldsymbol T\end{bmatrix}
=\begin{bmatrix}\boldsymbol \ell\\ \boldsymbol m\end{bmatrix},
\qquad
\boldsymbol A_{ij}=\int_\Omega\frac{1}{k}\boldsymbol r_i\cdot\boldsymbol r_j ,
\quad
\boldsymbol C_{ce}=-\int_{T_c}\nabla\cdot\boldsymbol r_e .
$$

The **same principle** produces the three methods of the series, depending on where the enrichment is
put:

| discrete fields | elimination of $\boldsymbol q$ gives | conforming in |
|---|---|---|
| $\boldsymbol q_h$ **element-local** constant, $T_h\in$ CR | $\boldsymbol q_h=\overline{\nabla_hT_h}$; the CR bilinear form $a_h(T_h,\delta T_h)=\sum_T\int_T\nabla T_h\cdot\nabla\delta T_h$ | $H^1$ **non**-conforming |
| $\boldsymbol q_h\in\mathcal{RT}_0$, $T_h\in\mathcal{P}_0$ | Schur complement $-\boldsymbol C\boldsymbol A^{-1}\boldsymbol C^{T}$ on $T_h$ (cell-centred scheme) | $H(\mathrm{div})$ conforming |
| RT0 $+$ multiplier $\lambda\in\mathcal{P}_0(\partial T)$ (= edge constants = **CR dofs**) | condensation gives exactly $\boldsymbol K^{CR}$ ([[Raviart-Thomas element]] Thm 5.1) | both, through $\lambda$ |

The first row is the reason CR appears in this note: with $T_h$ linear per element, $\nabla T_h$ is
element-constant, so the "flux" of the CR method is literally the *element-local* gradient — the H–R
flux equation is satisfied pointwise inside each element and violated across edges, which is the
non-conforming crime. The third row shows that the CR stiffness matrix is also the Schur complement of
a *conforming* H–R discretization once the trace is promoted to an unknown (hybridization). The circle
closes: **CR, RT0 and hybrid RT0 are three faces of the same two-field principle**, and this is why they
share the "one dof per edge" data layout.

### 7.1 Consistency check of (B1)–(B2) with $T$ in CR and $\boldsymbol q$ in RT0

The pairing "CR temperature + RT0 flux" is tempting because both discretizations have **one dof per
edge** (the same data structure, the same edge table). It is worth checking the two properties the
theory needs — consistency and stability — because the outcome is asymmetric:
**consistency holds exactly, stability fails completely.**

The statements are numbered (i)–(iv) below and proved in (v); the tables are evidence, not proofs.

**(i) Consistency: exact, no Strang term.** The continuous test spaces of (B1)–(B2) are
$H(\mathrm{div})\times L^2$, and

$$
\mathcal{RT}_0\subset H(\mathrm{div}),\qquad
CR\subset L^2
$$

— the CR element is *non-conforming only with respect to $H^1$* (its jumps are in the trace of the
gradient), which is irrelevant here: the mixed form never takes a gradient of $T$, only $\int T\nabla\cdot\boldsymbol r$
with $\nabla\cdot\boldsymbol r\in L^2$. Hence the exact solution satisfies (B1)–(B2) with *discrete*
test functions identically:

$$
R_1(\boldsymbol r_e)=\int\frac1k\boldsymbol q\cdot\boldsymbol r_e-\int T\nabla\cdot\boldsymbol r_e
+\int_{\Gamma_D}T_D\boldsymbol r_e\cdot\boldsymbol n=0\ \ \forall e,
\qquad
R_2(w_h)=-\int\nabla\cdot\boldsymbol q\,w_h-\int f w_h=0\ \ \forall w_h\in CR .
$$

There is no consistency functional at all (contrast with the *jump* functional of
[[Variational consistency of the CR element]]: the CR element is the offending party there only because
it is used as an $H^1$ *trial* space). Measured: $\max_e\lvert R_1\rvert=6.7\times10^{-16}$,
$\max\lvert R_2\rvert=0$.

**(ii) Stability: fails, and the failure does not disappear under refinement.** Since
$\nabla\cdot\mathcal{RT}_0=\mathcal{P}_0$ (constant per element), both equations see only the element
**average** of the temperature:

$$
\int_T T_h\,\nabla\cdot\boldsymbol r_i=\frac{\sigma_i}{A_T}\int_T T_h
=\frac{\sigma_i}{3}\bigl(m_1+m_2+m_3\bigr),
\qquad
\int_T\nabla\cdot\boldsymbol q_h\,\theta_I=\frac13\sum_{i}\sigma_i\Phi_i ,
$$

so the three CR dofs of an element enter (B1) only through their sum, and the three (B2) rows of an
element are *identical*. Consequently every $m$ with

$$
\sum_{I\in T}m_I=0\quad\forall T,\qquad m_I=0\ \text{on }\partial\Omega
$$

(the **shape modes** of the CR element, i.e. its $L^2$-orthogonal complement of $\mathcal{P}_0$)
satisfies $M(0,m)=0$: the assembled matrix is singular. Equivalently, the discrete inf-sup condition

$$
\sup_{\boldsymbol r_h\in\mathcal{RT}_0}\frac{\int q_h\,\nabla\cdot\boldsymbol r_h}{\lVert\boldsymbol r_h\rVert_{H(\mathrm{div})}}
\ \ge\ \beta\,\lVert q_h\rVert_{L^2}\qquad\forall q_h\in Q_h
$$

cannot hold with any $\beta>0$: the supremum only sees $\Pi_0q_h$, so it vanishes on the shape modes
($\ker(\mathrm{div})^{0}\supset\{q:\int_Tq=0\}$). A richer temperature space than $\mathcal{P}_0$ adds
exactly the modes that the $P_0$ divergence cannot see.

Measured on the four meshes of `msh/` (`cargo run --bin rt0_cr_mixed`):

| mesh | triangles | interior CR dofs | $\dim\ker$ | $\max_e\lvert R_1\rvert$ | LU | $\lVert T_h-T\rVert_0/\lVert T\rVert_0$ |
|---|---|---|---|---|---|---|
| `patchtest_n4` | 32 | 40 | $9$ | $6.7\times10^{-16}$ | "succeeds" | $3.2\times10^{4}$ |
| `patchtest_irregular_n4` | 54 | 73 | $19$ | $6.7\times10^{-16}$ | fails | — |
| `patchtest_n8` | 128 | 176 | $49$ | $6.7\times10^{-16}$ | "succeeds" | $1.3\times10^{4}$ |
| `patchtest_irregular_n8` | 238 | 341 | $103$ | $8.9\times10^{-16}$ | fails | — |

The kernel dimension is of the order of the number of interior vertices and **grows proportionally to
the number of elements**, i.e. like $h^{-2}$: refinement does not cure it. When the LU happens to return a vector it satisfies the equations to
$10^{-12}$ but carries an arbitrary kernel component, so the computed temperature is garbage
(relative $L^2$ error $\sim10^4$).

**(v) Proofs of (i) and (ii).** The numerical evidence above is not needed: both statements are
elementary, and they are proved here. Notation (all of it already introduced): edge dofs of $\mathcal{RT}_0$
are $\Phi_e=\int_{\Gamma_e}\boldsymbol q_h\cdot\boldsymbol n_e\,ds$ with the canonical normal, the assembled
local field is $\boldsymbol q_h|_T=\sum_{e\ni T}\sigma_e\Phi_e\boldsymbol\psi_e$ with orientation signs
$\sigma_e^{(T)}=\pm1$, and $\nabla\cdot\boldsymbol q_h|_T=\frac{1}{A_T}\sum_{e\ni T}\sigma_e^{(T)}\Phi_e$
([[Raviart-Thomas element]] §3.2); the CR dofs are $m_I$ and $\theta_I=1-2\lambda_I$
([[Crouzeix-Raviart triangular element]] §4).

> **Lemma 7.1.1 (both equations see only the element averages).** For every element $T$, every edge
> $e\subset\partial T$ and every local CR dof $I$,
> $$
> \int_T\theta_I\,d\Omega=\frac{A_T}{3},\qquad
> \int_T\theta_I\,\nabla\cdot\boldsymbol r_e\,d\Omega=\frac{\sigma_e^{(T)}}{3},\qquad
> \int_T\nabla\cdot\boldsymbol q_h\,\theta_I\,d\Omega=\frac13\sum_{e\subset\partial T}\sigma_e^{(T)}\Phi_e ,
> $$
> and in particular the last two are **independent of $I$**. Moreover
> $\frac{1}{A_T}\int_TT_h\,d\Omega=\frac13\sum_{I\in T}m_I$.
>
> *Proof.* $\int_T\lambda_I=A_T/3$ (CR notes §2), hence $\int_T\theta_I=A_T-\frac23A_T=\frac{A_T}{3}$
> and $\sum_I\theta_I=3-2\sum_I\lambda_I=1$; consequently
> $\frac{1}{A_T}\int_TT_h=\frac{1}{A_T}\sum_Im_I\int_T\theta_I=\frac13\sum_Im_I$. Since
> $\nabla\cdot\boldsymbol r_e=\sigma_e^{(T)}/A_T$ and $\nabla\cdot\boldsymbol q_h$ are **constants** on $T$,
> both remaining identities follow by multiplying by $\frac{A_T}{3}$. $\blacksquare$

> **Theorem 7.1.2 (the coupling block factorises through the element-average operator).** Let
> $$
> \mathcal E:Q_h\longrightarrow\mathbb R^{F},\qquad
> (\mathcal E q)_T:=\frac{1}{A_T}\int_Tq\,d\Omega=\frac13\sum_{I\in T}m_I
> $$
> ($F$ = number of elements) and let $B\in\mathbb R^{n_{\text{flux}}\times F}$ be the divergence incidence
> matrix, $B_{\boldsymbol r_e,T}=\sigma_e^{(T)}$ ($=0$ unless $e\subset\partial T$). Then the assembled
> matrix of (B1)–(B2) is
> $$
> \boldsymbol M=\begin{bmatrix}\boldsymbol K & -B\mathcal E\\[2pt]
> -\mathcal E^{\!\top}B^{\!\top} & \boldsymbol 0\end{bmatrix},
> \qquad
> \boldsymbol K_{ef}=\int_\Omega\frac1k\boldsymbol r_e\cdot\boldsymbol r_f\,d\Omega\ \ (\text{symmetric, positive definite}),
> $$
> hence
> $$
> \boxed{\ \ker\boldsymbol M\ \supset\ \ker\mathcal E\times\{0\}\ }
> $$
> and $\boldsymbol M$ is singular whenever $\mathcal E$ is not injective.
>
> *Proof.* By Lemma 7.1.1 the (B1) row of $\boldsymbol r_e$ receives from element $T$ the CR contribution
> $-\sum_{I\in T}m_I\frac{\sigma^{(T)}_e}{3}=-\sigma^{(T)}_e\cdot\frac13\sum_{I\in T}m_I=-\sigma^{(T)}_e(\mathcal E m)_T$;
> summing over the (one or two) elements $T\ni e$ gives $-(B\mathcal E m)_{\boldsymbol r_e}$. The (B2) row of
> the CR dof $I\in T$ receives from $T$ the flux term
> $-\frac13\sum_{e\subset\partial T}\sigma_e^{(T)}\Phi_e=-(B^{\!\top}\boldsymbol\Phi)_T$, the **same value for
> each of the three** dofs — i.e. the block is $-\mathcal E^{\!\top}B^{\!\top}$, the transpose of the (B1)
> coupling, as it must be for a Galerkin saddle point. The
> $\boldsymbol K$-block is the RT0 mass matrix. Taking $\Phi=0$ and $m\in\ker\mathcal E$ gives
> $\boldsymbol M(0;m)=0$. $\blacksquare$

> **Corollary 7.1.3 (singular for every mesh; the kernel grows like $h^{-2}$).** Let $n_{CR}^{\mathrm{int}}$
> be the number of *unknown* (interior) CR dofs. Since $\operatorname{rank}\mathcal E\le F$,
> $$
> \dim\ker\mathcal E\ \ge\ n_{CR}^{\mathrm{int}}-F .
> $$
> For a simply connected triangulation, Euler's formula gives $E=V+F-1$, and $E_b=V_b$ (boundary edges
> $=$ boundary vertices), so $n_{CR}^{\mathrm{int}}=E-E_b=V-V_b-1+F$ and
> $$
> \dim\ker\boldsymbol M\ \ge\ V_{\mathrm{int}}-1\ \ge\ 1
> $$
> whenever the mesh has an interior vertex, with $V_{\mathrm{int}}-1\sim h^{-2}$: **refinement never
> removes the kernel.** $\blacksquare$
> (Numerical check of the bound: `patchtest_n4` has $F=32$, $n_{CR}^{\mathrm{int}}=40$, so
> $\dim\ker\ge8$; measured $9$. `patchtest_n8`: $\ge48$, measured $49$.)

> **Corollary 7.1.4 (the discrete inf-sup constant is zero — with a criterion).** For any $q_h\in Q_h$ with
> $(\mathcal E q_h)_T=0$ for all $T$ (such $q_h\ne0$ exist by Corollary 7.1.3),
> $$
> b(\boldsymbol r_h,q_h)=-\int_\Omega q_h\,\nabla\cdot\boldsymbol r_h\,d\Omega
> =-\sum_T\nabla\cdot\boldsymbol r_h\big|_{T}\int_Tq_h\,d\Omega=0\qquad\forall\boldsymbol r_h\in W_h ,
> $$
> because $\nabla\cdot\boldsymbol r_h$ is constant on each element. Hence
> $$
> \beta_h:=\inf_{q_h\in Q_h}\sup_{\boldsymbol r_h\in W_h}
> \frac{b(\boldsymbol r_h,q_h)}{\lVert\boldsymbol r_h\rVert_{H(\mathrm{div})}\lVert q_h\rVert_{L^2}}=0 .
> $$
> The obstruction is exactly the one the Fortin criterion detects (Remark 7.1.8): the projection
> underlying (RT0, ·) commutes with $P_{\mathcal P_0}$ only, so $Q_h$ may not contain anything beyond
> element constants. Conversely, when $\mathcal E$ *is* injective — i.e. when the element averages form a
> complete set of functionals on $Q_h$, which for polynomial spaces means $\dim Q_h=F$, one constant per
> element — the standard theory applies, and that is precisely the $(\mathcal{RT}_0,\mathcal P_0)$ pair of
> §4. $\blacksquare$

> **Theorem 7.1.5 (consistency is exact: no Strang term).** Let $T\in H^2(\Omega)$,
> $\boldsymbol q=-k\nabla T$, $f=k\Delta T\in L^2$, $T|_{\Gamma_D}=T_D$ and $\Gamma_N=\emptyset$ (or
> $\boldsymbol q\cdot\boldsymbol n=q_n$ on $\Gamma_N$). Then
> $$
> R_1(\boldsymbol r_h):=\int_\Omega\frac1k\boldsymbol q\cdot\boldsymbol r_h-\int_\Omega T\nabla\cdot\boldsymbol r_h
> +\int_{\Gamma_D}T_D\,\boldsymbol r_h\cdot\boldsymbol n=0\qquad\forall\boldsymbol r_h\in W_h,
> $$
> $$
> R_2(w_h):=-\int_\Omega\nabla\cdot\boldsymbol q\,w_h-\int_\Omega f\,w_h=0\qquad\forall w_h\in Q_h .
> $$
>
> *Proof.* For $\boldsymbol r\in H(\mathrm{div})$ and $T\in H^1$ the product rule
> $\nabla\cdot(T\boldsymbol r)=T\nabla\cdot\boldsymbol r+\boldsymbol r\cdot\nabla T$ holds in
> $\mathcal D'(\Omega)$ and $T\,\boldsymbol r\cdot\boldsymbol n\in H^{-1/2}(\partial\Omega)$, so
> $\int_\Omega\nabla T\cdot\boldsymbol r=-\int_\Omega T\nabla\cdot\boldsymbol r+\int_{\partial\Omega}T\,\boldsymbol r\cdot\boldsymbol n$.
> Testing the Fourier law $\frac1k\boldsymbol q+\nabla T=\boldsymbol 0$ with $\boldsymbol r$ therefore gives
> $$
> 0=\int_\Omega\Bigl(\frac1k\boldsymbol q+\nabla T\Bigr)\cdot\boldsymbol r
> =\int_\Omega\frac1k\boldsymbol q\cdot\boldsymbol r-\int_\Omega T\nabla\cdot\boldsymbol r
> +\int_{\partial\Omega}T\,\boldsymbol r\cdot\boldsymbol n ,
> $$
> and $\int_{\partial\Omega}T\boldsymbol r\cdot\boldsymbol n=\int_{\Gamma_D}T_D\boldsymbol r\cdot\boldsymbol n$
> since $\Gamma_N=\emptyset$ (otherwise the $\Gamma_N$ contribution cancels against the essential datum
> $\boldsymbol q\cdot\boldsymbol n=q_n$ imposed on $W_h$). This is exactly $R_1(\boldsymbol r)=0$. For $R_2$,
> $\nabla\cdot\boldsymbol q=-k\Delta T=-f$ holds a.e., so $R_2(w)=0$ for every $w\in L^2$.
> Finally $\mathcal{RT}_0\subset H(\mathrm{div})$ and $CR\subset L^2$, so both identities hold for all
> $(\boldsymbol r_h,w_h)\in W_h\times Q_h$, i.e. the consistency functional
> $\sup_{(\boldsymbol r_h,w_h)}\lvert R_1+R_2\rvert/\lVert(\boldsymbol r_h,w_h)\rVert$ vanishes identically.
> $\blacksquare$

> **Corollary 7.1.6 (why consistency is not enough).** Strang's/Brezzi's error bound reads
> $$
> \lVert\boldsymbol q-\boldsymbol q_h\rVert_{H(\mathrm{div})}+\lVert T-T_h\rVert_{L^2}
> \le C_{\beta_h}\Bigl(\inf_{\boldsymbol r_h\in W_h}\lVert\boldsymbol q-\boldsymbol r_h\rVert_{H(\mathrm{div})}
> +\inf_{w_h\in Q_h}\lVert T-w_h\rVert_{L^2}+\underbrace{\sup\frac{\lvert R_1+R_2\rvert}{\lVert\cdot\rVert}}_{=\,0\ \text{by Thm 7.1.5}}\Bigr),
> $$
> where $C_{\beta_h}$ is the constant of Brezzi's theorem, which blows up as
> $\beta_h\to0$ (the exact power of $\beta_h^{-1}$ is irrelevant here). So consistency is exact and the
> approximation errors are $O(h)$, yet $\beta_h=0$ by Corollary 7.1.4 makes the bound void; indeed
> 1. **the discrete solution is not unique**: if $(\boldsymbol q_h,T_h)$ solves the system, so does
>    $(\boldsymbol q_h,T_h+m)$ for every $m\in\ker\mathcal E$;
> 2. **no stability bound exists**: $\lVert T_h+\lambda m\rVert_{L^2}$ is unbounded in $\lambda$ while the
>    residual stays zero, so $\lVert T_h\rVert\le C\lVert\text{data}\rVert$ is false;
> 3. **the right-hand side need not be attainable**: $b$ need not lie in $\operatorname{range}\boldsymbol M$,
>    which is why the sparse LU either fails (irregular meshes) or returns an arbitrary kernel component
>    together with a solution of residual $10^{-12}$ but relative error $\sim10^{4}$ (uniform meshes). $\blacksquare$

> **Remark 7.1.7 (why the same element is innocent in the primal method).** In
> [[Variational consistency of the CR element]] the CR element is the *trial and test* space of a
> second-order problem, so its test functions must lie in $H^1_0$; since $CR\not\subset H^1_0$, integration
> by parts leaves the jump functional $\sum_e\int_e[\![v_h]\!]\partial_{\boldsymbol n}u\,ds=O(h)$. In the
> mixed form of §4.2 the CR space appears only in $\int_\Omega q_h\nabla\cdot\boldsymbol r$ with
> $q_h\in L^2$: *no derivative of $q_h$ is taken*, so no jump term can arise and consistency is exact.
> The primal consistency tests $\nabla$, the mixed one tests $\nabla\cdot$ — different operators,
> different consistency statements about the same element.

> **Remark 7.1.8 (the Fortin criterion, and why RT1 cures it).** A classical sufficient condition for
> $\beta_h\ge\beta>0$ is (Fortin): there is $\Pi_h:H(\mathrm{div})\cap\prod_TH^1(T)\to W_h$ with
> $\lVert\Pi_h\boldsymbol v\rVert_{H(\mathrm{div})}\le C\lVert\boldsymbol v\rVert_{H^1}$ and
> $\int_T\nabla\cdot(\boldsymbol v-\Pi_h\boldsymbol v)\,q\,d\Omega=0$ for every $q\in Q_h$.
> *Proof of the criterion.* Given $q_h\in Q_h$ pick $\boldsymbol v$ with $\nabla\cdot\boldsymbol v=q_h$,
> $\lVert\boldsymbol v\rVert_{H^1}\le C\lVert q_h\rVert_{L^2}$ (surjectivity of the divergence); then
> $b(\Pi_h\boldsymbol v,q_h)=\int_\Omega q_h\nabla\cdot\boldsymbol v=\lVert q_h\rVert^2_{L^2}$ and the
> stability of $\Pi_h$ gives $\beta\ge1/C$. $\blacksquare$
> The classical RT interpolant satisfies the commuting property with $P_{\mathcal P_k}$ only, so the
> criterion requires $Q_h|_T\subset\mathcal P_k(T)$. With $W_h=\mathcal{RT}_0$ ($k=0$) this excludes CR —
> Theorem 7.1.2/Corollary 7.1.4 again — while with $W_h=\mathcal{RT}_1$ ($\nabla\cdot\mathcal{RT}_1=\mathcal P_1$)
> it holds for CR as well, which is why the pair of §7.2 is stable.

**Literature status (added).** I am not aware of a paper that uses $(\mathcal{RT}_0,CR)$ as the mixed pair
of (B1)–(B2) — consistent with Corollary 7.1.4, which makes the pair unusable without stabilization. The
"CR + RT0" literature is instead about the *equivalence* of the two methods (Arnold–Brezzi 1985;
Brezzi–Fortin 1991; Marini 1985; Cockburn–Gopalakrishnan 2004), about postprocessing/error estimation
from a CR solution to an RT0 flux, and about discrete complexes containing both; a grouped list with the
caveat that its bibliographic details still need checking is in
[[Raviart-Thomas element]] §10.

**(iii) What the correct pairing is.** Since $\nabla\cdot\mathcal{RT}_k=\mathcal{P}_k$, the natural
partner of $\mathcal{RT}_k$ is $\mathcal{P}_k$-type: (RT0, $P_0$) is the stable pair of §9 (verified
$O(h)$ in both fields), and the CR enrichment of $T$ is a *strict loss* — its element averages reproduce
exactly the $P_0$ solution while the extra modes are invisible. To pair a richer $T$ one must enrich the
flux divergence as well (e.g. $\mathcal{RT}_1$, since $\nabla\cdot\mathcal{RT}_1=\mathcal{P}_1=CR|_T$ per
element, so the shape modes are no longer invisible — whether the pair $(\mathcal{RT}_1,CR)$ is
*uniformly* stable is a separate question that is **not** checked here); this is the mixed counterpart of
the moment rule of
[[Legendre moment discretization and the consistency remainder]] §10 — *a field may only be enriched in
the modes that the pairing can see*.

**(iv) The lesson.** Consistency and stability are independent requirements, and for mixed methods the
binding one is stability. The CR element is *consistent and convergent* as an $H^1$-nonconforming method
(the other notes) and *perfectly consistent but singular* in the mixed pairing above; neither statement
can be inferred from the element alone — only from the pair of spaces and the operator that couples them.

### 7.2 The cured pairing: $T$ in CR with $\boldsymbol q$ in RT1

Raising the flux to $\mathcal{RT}_1$ restores stability, because $\nabla\cdot\mathcal{RT}_1=\mathcal{P}_1$
and $CR|_T=\mathcal{P}_1(T)$: the divergence now *sees* the shape modes of §7.1(ii). Formally the Fortin
argument behind the classical $(\mathcal{RT}_k,\mathcal{P}_k)$ theory goes through — it needs only

$$
\text{(i)}\quad Q_h|_T\subset\mathcal{P}_k(T)
\qquad\text{and}\qquad
\text{(ii)}\quad \nabla\cdot\Pi_h^{\mathcal{RT}_k}\boldsymbol v=P_{\mathcal{P}_k}(\nabla\cdot\boldsymbol v)
$$

with $P_{\mathcal{P}_k}$ the *element-wise* $L^2$ projection; for $k=1$ both hold with $Q_h=CR$ (for $k=0$
condition (i) fails, which is exactly the singular case of §7.1). Hence (RT1, CR) is stable.

`cargo run --bin rt1_cr_mixed` verifies, on the meshes of `msh/`:

**(i) The element.** $\max_{i,k}\lvert\text{dof}_k(\psi_i)-\delta_{ki}\rvert=7\times10^{-16}$ (the RT1
moment dofs are unisolvent), the divergence theorem per basis function holds to $10^{-14}$, and the
$H(\mathrm{div})$ normal trace matches across interior edges to $3\times10^{-14}$ (both Legendre modes;
mode 0 flips with the edge orientation, mode 1 is invariant).

**(ii) Consistency is still exact.** $R_1=O(10^{-12})$, $R_2=0$ for polynomial data — as in §7.1, since
$CR\subset L^2$; nothing changes on that side.

**(iii) Stability is restored.** The very shape mode that annihilated the (RT0, CR) operator now gives

$$
\lVert M(0,\dots,0,m)\rVert=3.95\times10^{1}\quad\text{instead of }0 ,
$$

the system is non-singular and the LU solve has residual $5\times10^{-15}$.

**(iv) Convergence** ($T=\sin\pi x\sin\pi y+x^2$, $k=1$, all-Dirichlet natural data):

| $Q_h$ | quantity | $n=4$ | $n=8$ | $n=16$ | $n=32$ | order |
|---|---|---|---|---|---|---|
| CR | $\lVert\boldsymbol q-\boldsymbol q_h\rVert_0$ | $1.857\times10^{-1}$ | $5.034\times10^{-2}$ | $1.290\times10^{-2}$ | $3.259\times10^{-3}$ | **1.98** |
| CR | $\lVert T-T_h\rVert_0$ | $1.983\times10^{-2}$ | $4.686\times10^{-3}$ | $1.156\times10^{-3}$ | $2.883\times10^{-4}$ | **2.00** |
| CR | $\lVert\nabla\cdot\boldsymbol q_h+f\rVert_0$ | $5.487$ | $2.992$ | $1.524$ | $0.768$ | **0.99** |
| $\mathcal{P}_1$ disc. | $\lVert\boldsymbol q-\boldsymbol q_h\rVert_0$ | $5.527\times10^{-2}$ | $1.397\times10^{-2}$ | $3.511\times10^{-3}$ | $8.799\times10^{-4}$ | **2.00** |
| $\mathcal{P}_1$ disc. | $\lVert T-T_h\rVert_0$ | $1.802\times10^{-2}$ | $4.583\times10^{-3}$ | $1.151\times10^{-3}$ | $2.880\times10^{-4}$ | **2.00** |
| $\mathcal{P}_1$ disc. | $\lVert\nabla\cdot\boldsymbol q_h+f\rVert_0$ | $3.858\times10^{-1}$ | $9.779\times10^{-2}$ | $2.453\times10^{-2}$ | $6.139\times10^{-3}$ | **2.00** |

So both fields converge at $O(h^{2})$ with the CR multiplier — the **flux** error even improves by a
factor $\approx3.4$ over the $(RT_0,P_0)$ pair of §9 — but the **divergence** gains only one order.

**Why the divergence lags behind.** $g:=\nabla\cdot\boldsymbol q_h+f$ is piecewise $\mathcal{P}_1$ and
(B2) gives only $g\perp Q_h$:

$$
\int_\Omega g\,q_h\,d\Omega=0\qquad\forall q_h\in Q_h .
$$

With $Q_h=\mathcal{P}_1^{\text{disc}}$ this is $g\perp\mathcal{P}_1(T)$ for every $T$; since
$\nabla\cdot\boldsymbol q_h\in\mathcal{P}_1$, it forces $g=0$, i.e.
$\nabla\cdot\boldsymbol q_h=-\Pi_{\mathcal{P}_1}f$ **exactly** — which is why the control column
coincides with $\lVert f-\Pi_{\mathcal{P}_1}f\rVert_0=(0.386,\,0.098,\,0.025,\,0.0061)$ to all digits.
With $Q_h=CR$ the orthogonality is imposed only against the $n_e$-dimensional CR space instead of the
$3n_{el}$-dimensional piecewise-$\mathcal{P}_1$ space, so $g$ retains the component of $f$ orthogonal to
CR *within* the piecewise-$\mathcal{P}_1$ space; that component is the mesh-scale "checkerboard" part,
whose pairing with smooth $f$ is $O(h)$ — hence $\lVert g\rVert_0=O(h)$.

**Recommendation.** (RT1, CR) is a legitimate, stable, $O(h^2)$ scheme for $\boldsymbol q$ and $T$; use
it when the temperature is wanted on the edges (it is the multiplier, and it needs no continuity), but
if the *discrete divergence* must be second order (local conservation of order $h^2$), take
$\mathcal{P}_1^{\text{disc}}$ — or post-process $\nabla\cdot\boldsymbol q_h$ by the local
$\mathcal{P}_1$ projection, which restores the second order at no cost.

**Robustness caveat.** On the perturbed meshes of `msh/` the temperature keeps its $O(h^2)$ rate but the
flux rate degrades (1.24, 1.41, 1.09) and the divergence error stagnates ($5.19\to3.59$), whereas the
$\mathcal{P}_1^{\text{disc}}$ control is unaffected. The CR multiplier is therefore noticeably more
sensitive to mesh irregularity than the discontinuous one; this warrants a stability-constant study
before using (RT1, CR) on strongly distorted meshes.

> **Implementation notes (three traps, all found while building this check).**
> 1. **Source term**: the load of a *shared* temperature dof must be accumulated (`rhs[I] += …`) over the
>    two elements. Assignment keeps only the last element and silently halves the forcing — invisible in a
>    patch test, where $f=0$.
> 2. **Dirichlet data are natural in the mixed form**, so *no* temperature dof may be prescribed:
>    prescribing the boundary $T$ dofs as well as putting $T_D$ in the right-hand side of (B1) counts the
>    datum twice and drops the local divergence conditions of the boundary elements (observed: the
>    divergence error decayed only like $h^{0.6}$).
> 3. **Quadrature normalisation**: the commonly tabulated 7-point Dunavant weights sum to $1$, i.e. twice
>    the area of the reference triangle. Using them without the $1/2$ over-scales *every* volume integral
>    by two — harmless for the discrete solution (the factor cancels between matrix and load), fatal for
>    any residual that mixes volume and boundary integrals, such as $R_1$ of §7.1.

### Boundary conditions in practice

| datum | Form A | Form B (= RT mixed) | `ansatz` operator |
|---|---|---|---|
| $T=T_D$ on $\Gamma_D$ | essential on $T$ | **natural**: $\ell_e$ in the RHS, $-\frac{1}{\ell_e}\int_{\Gamma_e}T_Dds$ | `temperature_lagrangian` (multiplier $=$ flux), `temperature_penalty` |
| $\boldsymbol q\cdot\boldsymbol n=q_n$ on $\Gamma_N$ | natural | **essential** on the flux dofs: fix $\Phi_e=\ell_eq_n$ | `heat_flux` |
| Robin on $\Gamma_R$ | natural | natural: substitute $\boldsymbol r\cdot\boldsymbol n=h(T-T_\infty)$ in the boundary term of (B1) | `convection` ($h\int N_iN_j$, $hT_\infty\int N_i$) |

---

## 8. Error estimates

For a stable pair, Brezzi's theory gives

$$
\lVert\boldsymbol q-\boldsymbol q_h\rVert_{H(\mathrm{div})}+\lVert T-T_h\rVert_{L^2}
\le C\Bigl(\inf_{\boldsymbol r_h\in V_h}\lVert\boldsymbol q-\boldsymbol r_h\rVert_{H(\mathrm{div})}
+\inf_{w_h\in Q_h}\lVert T-w_h\rVert_{L^2}\Bigr),
$$

and for $(\mathcal{RT}_0,\mathcal{P}_0)$ on shape-regular meshes, with $\boldsymbol q\in H^1$,
$f\in H^1$:

$$
\lVert\boldsymbol q-\boldsymbol q_h\rVert_{L^2}\le Ch\lvert\boldsymbol q\rvert_{1},
\qquad
\lVert T-T_h\rVert_{L^2}\le Ch\lvert T\rvert_{1},
\qquad
\lVert f-\nabla\cdot\boldsymbol q_h\rVert_{L^2}=\lVert f-\Pi_{0}f\rVert_{L^2}\le Ch\lvert f\rvert_{1}.
$$

The last identity is exact — **local conservation holds cell by cell** — and it is the reason the RT0
divergence is called superconvergent: the discrete divergence is the $P_0$ projection of $f$,
independently of the flux error.

---

## 9. Numerical verification

`cargo run --bin hellinger_reissner_heat` solves Form B with $(\mathcal{RT}_0,\mathcal{P}_0)$,
$k=1$, on the uniform and irregular meshes of `msh/`, for the manufactured solution

$$
T=\sin(\pi x)\sin(\pi y)+x
\qquad\Longrightarrow\qquad
T_D\ne0\ \text{on }\partial\Omega,
\quad k\Delta T=f=-2\pi^{2}\sin(\pi x)\sin(\pi y),
$$

with the Dirichlet datum entering exactly as $-\frac{1}{\ell_e}\int_{\Gamma_e}T_D\,ds$ (verified
below), and the discrete flux oriented by the canonical edge normal (the two triangles of an edge have
opposite outward normals, so a sign $\sigma_{e,T}=\pm1$ accompanies every assembly entry — this is
the implementational trap of RT0).

| mesh | $n$ | $\lVert\boldsymbol q-\boldsymbol q_h\rVert_{0}$ | rate | $\lVert T-T_h\rVert_{0}$ | rate | $\lVert\nabla\cdot\boldsymbol q_h+f\rVert_{0}$ | rate |
|---|---|---|---|---|---|---|---|
| uniform | 4 | $4.998\times10^{-1}$ | — | $1.419\times10^{-1}$ | — | $2.543$ | — |
| uniform | 8 | $2.514\times10^{-1}$ | 0.99 | $7.157\times10^{-2}$ | 0.99 | $1.287$ | 0.98 |
| uniform | 16 | $1.259\times10^{-1}$ | 1.00 | $3.586\times10^{-2}$ | 1.00 | $6.453\times10^{-1}$ | 1.00 |
| uniform | 32 | $6.295\times10^{-2}$ | 1.00 | $1.794\times10^{-2}$ | 1.00 | $3.229\times10^{-1}$ | 1.00 |
| irregular | 8 | $2.420\times10^{-1}$ | — | $5.104\times10^{-2}$ | — | $9.162\times10^{-1}$ | — |
| irregular | 16 | $1.189\times10^{-1}$ | 1.02 | $2.541\times10^{-2}$ | 1.01 | $4.649\times10^{-1}$ | 0.98 |
| irregular | 32 | $5.899\times10^{-2}$ | 1.01 | $1.296\times10^{-2}$ | 0.97 | $2.340\times10^{-1}$ | 0.99 |

All three quantities converge at $O(h)$, as predicted by §8; the divergence column is the $P_0$
interpolation error of $f$ and converges with the same rate but a much larger constant (it does not
benefit from the flux superconvergence).

### 9.2 The multiplier principle (§4.3–4.4)

`cargo run --bin complementary_multiplier_heat` solves (E1)–(E3) with $\Gamma_N=$ top edge,
$\Gamma_D=$ the other three edges, $k=1$, and the manufactured solution
$T=\sin(\pi x)\sin(\pi y)+x^{2}$ (so $T_D$ is *quadratic* on the bottom edge and
$q_n=-k\partial_yT=+\pi\sin(\pi x)$ on $\Gamma_N$, with a non-zero first moment).

| mesh | $n$ | $\lVert\boldsymbol q-\boldsymbol q_h\rVert_0$ | rate | $\lVert T-T_h\rVert_0$ | rate | $\max_{\Gamma_N}\lvert\mu-T\rvert$ | rate |
|---|---|---|---|---|---|---|---|
| uniform | 4 | $5.194\times10^{-1}$ | — | $1.261\times10^{-1}$ | — | $3.243\times10^{-2}$ | — |
| uniform | 8 | $2.613\times10^{-1}$ | 0.99 | $6.340\times10^{-2}$ | 0.99 | $7.809\times10^{-3}$ | 2.05 |
| uniform | 16 | $1.309\times10^{-1}$ | 1.00 | $3.174\times10^{-2}$ | 1.00 | $1.875\times10^{-3}$ | 2.06 |
| uniform | 32 | $6.548\times10^{-2}$ | 1.00 | $1.588\times10^{-2}$ | 1.00 | $3.605\times10^{-4}$ | 2.38 |
| irregular | 8 | $2.456\times10^{-1}$ | — | $4.468\times10^{-2}$ | — | $4.666\times10^{-3}$ | — |
| irregular | 16 | $1.206\times10^{-1}$ | 1.00 | $2.254\times10^{-2}$ | 1.11 | $9.317\times10^{-4}$ | 1.56 |
| irregular | 32 | $5.982\times10^{-2}$ | 1.03 | $1.146\times10^{-2}$ | 0.99 | $2.717\times10^{-4}$ | 2.32 |

and, at $n=8$ on both meshes:

| check | uniform | irregular |
|---|---|---|
| $\lvert(\boldsymbol q,T)_{M_h}\text{-constant}-(\boldsymbol q,T)_{\text{strong}}\rvert$ | $(1.5,\,1.6)\times10^{-15}$ | $(2.5,\,1.1)\times10^{-15}$ |
| linear $\mu_h$: solvable? | **no** | **no** |
| $\max_e\lvert\int_{\Gamma_e}q_n\xi\,ds\rvert$ (must vanish for solvability) | $2.511\times10^{-2}$ | $2.511\times10^{-2}$ |

So: the multiplier form is *consistent* with the strong imposition (the fields agree to round-off, as
§4.4(a) predicts), the multiplier itself is the temperature with a second-order error, and the over-rich
multiplier is unsolvable exactly when the predicted moment of the data fails to vanish.

> **Implementation note (a trap worth recording).** The boundary terms of (E1) must be assembled for
> **all three** $\mathcal{RT}_0$ bases of the adjacent element, not only for the edge's own dof: the trace
> of $\boldsymbol r_i$ on a *neighbouring* edge of the same element has zero **mean** (it is zero through
> the Kronecker property) but a non-zero **first moment**, so a weighted integral
> $\int_{\Gamma_e}T_D\,\boldsymbol r_i\cdot\boldsymbol n\,ds$ does not vanish. Assembling only the own dof
> leaves an $O(h)$ inconsistency: in this example it inflated the flux error from
> $2.6\times10^{-1}$ to $9.6\times10^{-1}$ at $n=8$, and made the multiplier come out as $-T$ instead of
> $+T$. The helper `cr_verify::rt0_boundary_moments` assembles the correct three moments.

### A one-dimensional sanity check (exact, by hand)

In 1D the same discretization is classical. With $q_i$ constant per cell and $T_i$ the multiplier,
(B1) for a cell reduces to $q_{i+1/2}=-(T_{i+1}-T_i)/h$ and (B2) to
$q_{i+1/2}-q_{i-1/2}=-\int_{T_i}f$, i.e. the three-point finite-volume scheme. For $f\equiv1$,
$T(0)=T(1)=0$, $h=1/n$ its solution is $T_i=-\frac{h^{2}}{2}i(n-i)$, and

$$
-\frac{h^{2}}{2}i(n-i)=\frac{x_i(x_i-1)}{2}=T(x_i)\quad\text{with } x_i=ih ,
$$

so the mixed 1D values are **exact at the nodes** (and the fluxes are exact at the midpoints) for
constant $f$ — the 1D counterpart of the cell-by-cell exact conservation above.

---

## 10. Remarks and open items

1. **Two-sided bounds.** Because $\Pi_{HR}$ is a saddle, its value brackets the energy
   (Prager–Synge hypercircle): admissible fluxes give an upper bound, admissible temperatures a lower
   bound. This is the classical use of H–R and it is available here for error estimation, but it
   requires *equilibrated* fluxes ($\nabla\cdot\boldsymbol q_h+f=0$ exactly, which RT0 gives).
   The multiplier $\mu_h$ of §4.3 is the complementary object on $\Gamma_N$: a lower-bound
   (temperature-side) datum, which is what one needs to close the bounds on a Neumann boundary.
2. **Equal-order and stabilized pairs.** Form A with $(\mathcal{P}_1,\mathcal{P}_1)$ is unstable; the
   standard cure is to add a residual-based term, which is the same "least-squares" idea as the
   `Least square` section of [[Crouzeix-Raviart triangular element]] — the H–R functional plus
   $\frac{\gamma}{2}\int\lvert\frac1k\boldsymbol q+\nabla T\rvert^{2}$ is a stabilized H–R principle.
3. **Why this note is the parent of the others.** The moment dictionary of
   [[Legendre moment discretization and the consistency remainder]] says: the H–R flux equation is
   tested against fluxes whose edge moments are the dofs; the CR remainder prices the missing
   $\mu_{e,1}$, the RT0 interpolant misses it too, and RT1 supplies it. Raising the order of either
   method means supplying that moment.
4. **Open items.** (a) the Hu–Washizu three-field form and its relation to the mixed/enhanced
   variational formulations used for nearly-incompressible media; (b) a hypercircle-based a-posteriori
   estimator on the `msh/` meshes; (c) the extension to transient conduction, where H–R must be paired
   with the capacity operator of `ansatz::operation::heat`.

---

## 11. Literature annotations (which reference supports which result)

Citation keys as in [[Raviart-Thomas element]] §10.

| result of this note | literature | comment |
|---|---|---|
| §2 the four principles; §3 $\Pi_{HR}$ and its stationarity | **[He14]**, **[Re50]**, **[BF91]** Ch. VI, **[ABD84]** | the two-field principle in its mechanical form is Hellinger's/Reissner's; Brezzi–Fortin and PEERS are the FE treatment |
| §3 Legendre transform $g^{*}$, min–sup equivalence | **[He14]**, **[Re50]** | convex duality is how the two-field form is derived |
| §4.2 Form B (mixed weak form), Dirichlet natural / Neumann essential | **[Bre74]**, **[RT77]**, **[BF91]**, **[AB85]** | the classical mixed formulation |
| §5 inf-sup conditions, kernel coercivity, saddle structure, LU not Cholesky | **[Bre74]**, **[BBF13]**, **[Ste84]** | Brezzi's theory |
| §6 "mixed = consistent by construction" | **[BF91]** Ch. III, **[AB85]** | the mixed method has no consistency error |
| §7.1 Corollary 7.1.4 ($\beta_h=0$ for any element-wise constant divergence) | — (own) | the criterion behind it is Fortin's, **[For77]** |
| §7.1 Theorem 7.1.2 / Corollary 7.1.3 ($\ker\ge V_{\mathrm{int}}-1$) | — (own) | elementary, but I have not seen it stated for this pair |
| §7.1 Remark 7.1.7 (same element innocent in the primal method) | **[AB85]** | published comparison of the two consistency errors |
| §7.1 Remark 7.1.8 (Fortin criterion; why RT1 cures it) | **[For77]**, **[BF91]**, **[BBF13]** | the criterion and the $(\mathcal{RT}_k,\mathcal{P}_k)$ stability |
| §7.2 (RT1, CR) stable, $O(h^2)$ for $\boldsymbol q,T$, $O(h)$ for $\nabla\cdot\boldsymbol q_h$ | — (own), with **[BBF13]** | the rates follow from the $(\mathcal{RT}_1,\mathcal{P}_1)$ theory; the divergence deficiency is the CR-space effect proved here |
| §9.2 multiplier = temperature; over-rich multiplier obstruction | — (own), with **[CG04]**, **[BF91]** | the identification $\mu=T|_{\Gamma_N}$ is classical for H–R; the "multiplier may not be richer than the flux trace" rule is the notes' moment-visibility statement |
| §10 item 1 (hypercircle bounds) | **[DM99]**, **[EV15]** | equilibrated fluxes ⇒ two-sided bounds |
