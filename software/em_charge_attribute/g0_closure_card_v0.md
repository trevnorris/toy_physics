# G0 shared minimal closure card v0 — DRAFT FOR REVIEW

**Status:** `DRAFT v0`, authored 2026-07-20 for the Codex → Grok → Codex gauntlet and human review. This is a **postulated closure**, not a derivation and not a final model. A future card may replace any `[POSTULATE]`; no entry below may be changed after a force sign is seen.

**Preregistered label:** `G0-v0 / fixed-Σ / fixed-source / two_sided_R_w_postulate / E3-lead / E2-contrast / stationary-as-V→0` `[CONVENTION]`, dimension `1`; provenance: the branch, ensemble, ambient, and anti-frozen discipline required by `g0_minimal_closure_card_scope.md` §§2–5.

**What varies:** only the sleeve-contact operator `𝔅` in §8: E3 is lead and E2 is contrast. All equations, coefficients, geometry, mouth data, source laws, return control, and force bookkeeping in §§1–7 are byte-identical parent data.

**Warning:** this card does not compute, presume, or name an attract/repel result. In particular, it contains no Maxwell field, gauge field, Gauss law, point source, `j∝sV`, Coulomb sign, or like-repel prior. `pathA_38` appears only in the regression paragraph (§10.3), after the parent functional and source sign have been fixed independently.

## 1. Dimensions, fields, domains, and fixed geometry

Write dimensions as `(L,T,M)` and let `[E]=(2,-2,1)`. The spatial coordinate is `y=(x,w)∈R³×R`, `r_i=|x-X_i|`, body orientation is `s_i=±1`, and `z_i=s_iw`. The medium-rest frame and outward-normal convention are used throughout. A jump is `[q]_-^+=q(0^+)-q(0^-)`.

The retained fields are bulk number density `ρ` `(-4,0,0)`, bulk carrier phase `θ` `(0,0,0)`, wall amplitude `r_B∈[0,1]` `(0,0,0)`, parent localized scalar `H` `(-1,0,0)`, reduced scalar `h` `(0,0,0)`, and longitudinal brane displacement `u_L` `(1,0,0)`. The stationary bulk phase is

`θ(y,t)=−μ_∞t/ℏ+φ(y)`, `v=(ℏ/m)∇₄φ`.

This is stationary even though `∂tθ≠0`.

Under the declared body-conjugating `R_w`, `s→−s`; `ρ,θ,r_B` transform by ordinary pullback, while the retained **orientation-channel** fields transform as `H(x,w)→−H(x,−w)`, `h(x)→−h(x)`, and `u_L(x)→−u_L(x)`. This makes both `Q_χH` and `C_hu∇u_L·∇h` invariant. It is `[POSTULATE]`, field dimensions as above; provenance: the minimal common parity representation compatible with a nonzero `C_hu`. A charge-even longitudinal `u_L` or scalar `H` would instead forbid this shared constant mixing on an exactly `R_w`-symmetric ambient and is flagged for review.

### 1.1 Frozen, body-attached sleeve

Before collar rounding, the sleeve surface for body `i` is the cylinder plus hemispherical cap

```text
Σ_i^cyl : r_i=a,                     0≤z_i≤L−a,
Σ_i^cap : r_i²+[z_i−(L−a)]²=a²,      L−a≤z_i≤L.
```

The open mouth is `M_i={w=0,r_i<a}` and its collar is `C_i={w=0,r_i=a}`. At the collar, hard `min/max` intersections are replaced by

`smax_δ(p,q)=δ log(exp(p/δ)+exp(q/δ))`, `smin_δ(p,q)=−smax_δ(−p,−q)`

with `δ=ℓ`; this fixes a smooth signed-distance level set without adding a fit parameter. `Σ_i` means the resulting rounded lateral-plus-cap surface, never the mouth. The lumen and ambient remain connected through `M_i`.

| Entry | Value | Tag, dimension, one-line provenance |
|---|---:|---|
| mouth radius `a` | base unit; `a*=1` | `[CONVENTION]`, `(1,0,0)`; required nondimensionalization |
| sleeve length `L/a` | `37/20` | `[CALIBRATED]`, `1`; frozen U1 branch ansatz |
| wall/healing width `ℓ/a=ε_r` | `1/20` | `[CALIBRATED]`, `1`; frozen handoff scale |
| brane ordered-slab width `W_B/a` | `1/5=4ℓ/a` | `[POSTULATE]`, `1`; smallest separate slab scale used to regularize the two-sided wall |
| collar rounding `δ/a` | `1/20` | `[CONVENTION]`, `1`; reuse `ℓ`, so no new scale |
| frozen `Σ_i` | formula above, translated with `X_i` | `[POSTULATE]`, surface dimension `(3,0,0)`; minimal smooth fixed-geometry MVP |

The far brane is the ordered slab `|w|<W_B/2`; the sleeve is its body-attached tubular continuation on side `s_i`. The prescribed wall mid-surface `Γ_Σ` is the boundary of that slab-plus-sleeve ordered region. It is an imposed level set, not a lab-pinned object.

For every occurrence below, the fixed reference amplitude is explicit: if `d_Σ(y)` is the signed distance (positive inside the ordered slab-plus-sleeve region) to the rounded surface above, then

`r_Σ(y)=[1+exp(−d_Σ(y)/ℓ)]⁻¹`.

The far-brane reference used in subtractions is

`r_bg(w)={tanh[(w+W_B/2)/ℓ]−tanh[(w−W_B/2)/ℓ]}/{2tanh[W_B/(2ℓ)]}`.

It is even, normalized to `r_bg(0)=1`, and decays on both sides. These profile definitions are `[POSTULATE]`, dimension `1`; provenance: the minimal smooth logistic representative of the already frozen level-set geometry. `S_hold`, rather than an unstated force, maintains this finite slab/sleeve against kink–antikink or cap motion.

### 1.2 Where the frozen-shape reaction is booked

The shared conservative holding functional is

`S_hold=−∫dt ∫_{Γ_Σ} d³A λ_Σ(r_B−1/2)`.

`λ_Σ` has dimension `E/L³=(-1,-2,1)`. Variation in `λ_Σ` fixes the wall mid-surface; variation in `r_B` supplies its conormal jump. Split `Γ_Σ=Γ_B∪(∪_iΓ_Σ,i)`: the far background `Γ_B` is part of the ambient subtraction, while each sleeve/collar piece `Γ_Σ,i(X_i)` moves with its body. When `X_i` varies, `−δS_hold/δX_i` from that local piece is included in **`F_var`**. The equal-and-opposite shape-maintaining reaction belongs to the body scaffold/holding ledger. It is not an unreported laboratory support and not `F_𝔅`.

This assignment is `[CONVENTION]`, force dimension `(1,-2,1)`; provenance: the fixed-`Σ` requirement and the U1 rule that passive holonomic reactions must not be counted twice. Alternatives—finite bending/anchoring or a free `Σ`—are flagged in REVIEW NOTES.

## 2. Shared conservative parent functional

The conservative action is

`S_cons=S_bulk+S_χ+S_scalar+S_mouth+S_hold+S_geon,const`,

where `S_scalar` is implemented in **one** of two mathematically matched forms: the full parent `H` form or its normalized zero-mode reduction in §2.4. They are never summed.

The drain/return in §6 is deliberately non-variational. A background subtraction is made inside every infinite-volume static integral.

For avoidance of a sign ambiguity, `S_mouth:=−∫dt Ω_mouth`; every other static energy displayed below likewise enters the Lagrangian with a minus sign. This is `[CONVENTION]`, action dimension `(2,-1,1)`; provenance: standard source-versus-energy Legendre bookkeeping.

### 2.1 Bulk carrier and the quantum convention package

Use the number-density Madelung action

```text
S_bulk = ∫dt d⁴y [
  −ℏρ ∂tθ
  −(ℏ²ρ/2m)|∇₄θ|²
  −U(ρ)
  −(ℏ²/8mρ)|∇₄ρ|² ],
U(ρ)=(K/4)ρ⁵,       P(ρ)=Kρ⁵.
```

This entry is `[CALIBRATED]`, action dimension `(2,-1,1)`; provenance: the carrier, EOS, and quantum-pressure package already fixed in handoff §2.1. It is not a new electric-sector fit.

The exact quantum convention retained by this card is

```text
Π_Q,ij = −(ℏ²/4m) ρ ∂i∂j lnρ,
ε_Q    =  (ℏ²/8mρ)|∇ρ|²,
j_ε^Q  = [ρQ_B−ε_Q]v − (ℏ²/4m)(∂tρ/ρ)∇ρ,
Q_B    = −(ℏ²/2m)(∇²√ρ/√ρ),
ε      = U+½mρ|v|²+ε_Q,
j_ε    = (ε+P)v+j_ε^Q.
```

`Π_Q`, `ε_Q`, and `j_ε^Q` have dimensions `E/L⁴=(-2,-2,1)`, `E/L⁴`, and `E/(L³T)=(-1,-3,1)`, respectively. The package is `[CONVENTION]`; provenance: the canonical Schrödinger energy current paired with the handoff's fixed `Π_Q`. It obeys

`∂jΠ_Q,ij=−(ℏ²/2m)ρ ∂i(∇²√ρ/√ρ)`

and the source-free local energy identity. Improvement terms are fixed to zero.

### 2.2 Regular wall-amplitude action

The amplitude-only wall action is

```text
S_χ = ∫dt d⁴y [
  ½ Z_χ (∂t r_B)²
  −½ κ_χ |∇₄r_B|²
  −(λ_χ/4) r_B²(1−r_B)² ].
```

It is polynomial and nonsingular at `r_B→0`; no bare phase variable occurs there. The two minima are `r_B=0,1`. With `λ_χ=2κ_χ/ℓ²`, a planar unconstrained wall has logistic width `ℓ` and surface tension `σ_χ=κ_χ/(6ℓ)`.

| Entry | Dimension | Star value | Tag and one-line provenance |
|---|---:|---:|---|
| `Z_χ` | `(-2,0,1)` | `1` | `[POSTULATE]`; constant minimal kinetic normalization |
| `κ_χ` | `(0,-2,1)` | `1` | `[POSTULATE]`; constant isotropic gradient stiffness |
| `λ_χ` | `(-2,-2,1)` | `800` | `[POSTULATE]`; fixed by `2κ_χ/ℓ²`, not by a force sign |
| `V_χ(r_B)` | `(-2,-2,1)` | `(800/4)r_B²(1−r_B)²` | `[POSTULATE]`; lowest regular double well with minima 0 and 1 |
| `κ_bend`, `κ_anchor`, collar two-surface tension | `E/L=(1,-2,1)`; `E/L³=(-1,-2,1)`; `E/L²=(0,-2,1)` | `0,0,0` | `[POSTULATE]`; fixed-`Σ` MVP removes free-shape dynamics; `S_hold` carries the reaction explicitly |

Here and below a star means the units in §7. The Euler–Lagrange equation away from `Γ_Σ` is

`Z_χ ∂t²r_B−κ_χ∇₄²r_B+(λ_χ/2)r_B(1−r_B)(1−2r_B)=0`.

At a smooth wall piece `r_B` is continuous and `[κ_χ∂_n r_B]_-^+=λ_Σ`. At the rounded collar, `r_B` is continuous and the sum of incoming conormal fluxes equals the same holding multiplier; no extra collar functional is hidden.

### 2.3 Localized parent `H` and its reduced `h` mode

The distinct electric scalar is represented first as a parent field `H(x,w,t)` with dimension `L⁻¹`. Define

```text
O_⊥ = −∂w² + [4−6 sech²(w/ℓ)]/ℓ²,
S_H = ½∫dt d⁴y {2M₄(∂tH)²−2K₄[|∇xH|²+H O_⊥H]}.
```

This entry is `[POSTULATE]`, action dimension `(2,-1,1)`; provenance: the lowest factorable reflection-even Pöschl–Teller operator with one gapless localized ground mode and a gapped continuum. It was selected without using a two-body sign.

`O_⊥=A†A`, `A=∂w+2tanh(w/ℓ)/ℓ`, so it is nonnegative. Its zero mode and norm are

`f₀(w)=1/[ℓ cosh²(w/ℓ)]`, `N₀=∫dw 2f₀²=8/(3ℓ)`.

Set `H=f₀h+H_⊥`, with `∫dw 2f₀H_⊥=0`. The zero-mode reduction gives `M_h=N₀M₄`, `K_h=N₀K₄`, `K_h=M_hc_E²`.

| Entry | Dimension | Star value | Tag and one-line provenance |
|---|---:|---:|---|
| `M₄` | `(0,0,1)` | `3ℓ*/8=3/160` | `[POSTULATE]`; chosen so the normalized reduced `M_h*=1` |
| `K₄` | `(2,-2,1)` | `3/160` | `[POSTULATE]`; `K₄=M₄c_E²` with `c_E*=1` |
| `M_h` | `(-1,0,1)` | `1` | `[POSTULATE]`; positive scalar inertia |
| `K_h` | `(1,-2,1)` | `1` | `[POSTULATE]`; positive scalar stiffness |
| `c_E` | `(1,-1,0)` | `1` | `[CONVENTION]`; time unit, not identified with `c_s` or `c_γ` |
| `H³`, `H⁴`, a nonlocal addition to `O_⊥`, and a bulk `H²` mass | coefficient/operator dimensions `E/L=(1,-2,1)`, `E=(2,-2,1)`, `L⁻²=(-2,0,0)`, `E/L²=(0,-2,1)` | `0,0,0,0` | `[POSTULATE]`; minimal linear gapless zero-mode sector |

In a parent-field solve, `S_scalar=S_H+S_u+S_mix`, where `h=P₀H:=N₀⁻¹∫dw 2f₀H` is used in `S_mix` below. In a far-field zero-mode solve, `H_⊥` is Schur-eliminated and `S_scalar` is replaced by the reduced action below. This matching rule is `[CONVENTION]`, action dimension `(2,-1,1)`; provenance: the required Tier-2 embedding and no-double-count discipline.

### 2.4 Phase representation and the `(u_L,h)` scalar block

**Chosen representation:** use `Q_s(u_L,h)` after eliminating the brane phase `θ_B`, and put the Schur inertia shift into `A_eff`; do not use bare `ρ_br` and do not also solve `(u_L,θ_B)`.

```text
A_eff = ρ_br + C_J²/κ_phase,
S_scalar^reduced = S_Lh = ∫dt d³x [
  ½A_eff (∂tu_L)² + ½M_h(∂th)²
  −½B_eff|∇u_L|² −½K_h|∇h|²
  −C_hu ∇u_L·∇h ].
```

The equivalent parent pieces are `S_u=∫dt d³x[½A_eff(∂tu_L)²−½B_eff|∇u_L|²]` and `S_mix=−∫dt d³x C_hu∇u_L·∇(P₀H)`. Thus the `M_h,K_h` terms occur in `S_H` **or** `S_Lh`, never both.

The primitive-sign bookkeeping is `K_θ=−κ_phase<0`. The numerical decomposition `ρ_br*=3/4`, `(C_J²/κ_phase)*=1/4` gives `A_eff*=1`; it records the pathA36 Schur shift without retaining another field.

| Entry | Dimension | Star value | Tag and one-line provenance |
|---|---:|---:|---|
| `A_eff` | `(-3,0,1)` | `1` | `[POSTULATE]`; exact eliminated-phase inertia used instead of bare `ρ_br` |
| `ρ_br/A_eff` | `1` | `3/4` | `[POSTULATE]`; exposes rather than hides the Schur increment |
| `(C_J²/κ_phase)/A_eff` | `1` | `1/4` | `[POSTULATE]`; nonzero pathA36-compatible eliminated-phase contribution |
| `B_eff` | `(-1,-2,1)` | `2` | `[POSTULATE]`; positive longitudinal stiffness; primitive `B=0` remains curl-only |
| `C_hu` | `(0,-2,1)` | `1/2` | `[POSTULATE]`; nonzero minimal sensitivity point, not the decoupling ablation |
| `g_θθB` | `E/L³=(-1,-2,1)` | `0` | `[POSTULATE]`; coefficient of a possible brane `θ−θ_B` interface energy, zero because `θ_B` is eliminated |

`θ_B` therefore has no independent background rotation, chemical potential, initial datum, or ninth R49 tangent. E3's “phase texture” in this card is the **bulk carrier phase `θ`**, never `θ_B`. The brane-order phase's only retained effect is `A_eff−ρ_br=1/4`.

This phase choice is load-bearing. The alternatives `(u_L,θ_B)` or a bare-`ρ_br` `Q_s` are listed in REVIEW NOTES; neither is silently identified with this block.

### 2.5 Geon constant in the electric MVP

The trapped geon contributes only `S_geon,const=−Σ_i∫dt E_g,i`, with `E_g,i` dimension `E`. It is a separation-independent additive constant and is subtracted from the static interaction/force functional. Its value is held out with the magnitude problem; every field-dependent geon coupling to `{r_B,H,h,u_L,ρ,θ}` is zero in G0-v0. This is `[CONVENTION]`; provenance: the geon supplies body rest content but no additional electric closure knob in the fixed-body static MVP.

## 3. Independent fixed-source mouth functional: the `χ↔h` mechanism

### 3.1 Annular support and orientation functional

The mouth annulus is the spherical brane shell

`A_i={x:a≤r_i≤a+ℓ}`, `η_i(x)=3·1_A/[4π((a+ℓ)³−a³)]`, so `∫d³x η_i=1` and `[η_i]=L⁻³`.

Let `o_ℓ(w)=w exp[−w²/(2ℓ²)]/(√(2π)ℓ³)`, an odd kernel of dimension `L⁻²`. For the frozen reference wall `r_Σ,+`, define

```text
N_χ = ∫d³x η_i ∫dw o_ℓ(w)[r_Σ,+²−r_bg²],
Q_χ[r_Σ,s] = N_χ⁻¹ ∫d³x η_i ∫dw o_ℓ(w)[r_Σ,s²−r_bg²].
```

Body conjugation gives `Q_χ[r_Σ,s]=s`; this is a normalized geometric orientation projection, not a winding and not a point charge. `N_χ` has dimension `L⁻¹`. This definition is `[CONVENTION]`; provenance: a compact odd projection is the minimal map from the fixed sleeve amplitude to a signed scalar source.

The orientation convention takes `o_ℓ(w)>0` for `w>0` and therefore `N_χ>0` for the `+w` reference sleeve. This fixes the source sign before any propagator is applied.

### 3.2 Functional, sign, and conjugate natural condition

The fixed-source mouth functional is

```text
Ω_mouth = Σ_i ∫d³x η_i(x)
  [ ½K_m H(x,0)² − J_m Q_χ[r_Σ,s_i] H(x,0) ].
```

The convention is: **positive `Q_χ` and positive `J_m` put `−J_mH` in the static energy**. This fixes only the source sign; it says nothing about the eventual force sign.

| Entry | Dimension | Star value | Tag and one-line provenance |
|---|---:|---:|---|
| `K_m` | `(4,-2,1)=E·L²` | `(ℓ/a)²=1/400` | `[POSTULATE]`; positive minimal Robin impedance |
| `J_m` | `(3,-2,1)=E·L` | `ℓ/a=1/20` | `[POSTULATE]`; unit reduced source chosen before regression |
| reduced `k_m=K_m/ℓ²` | `E` | `1` | `[POSTULATE]`; localized Robin coefficient |
| reduced `g_χh=j_m=J_m/ℓ` | `E` | `1` | `[POSTULATE]`; linear, local-in-`x`, orientation-odd fixed-source coupling |

Variation in `H` gives the independently defined Robin/Neumann jump

`2K₄[∂wH]_-^+ = η_i(K_mH−J_mQ_χ)` on `w=0`, supported only on `A_i`,

and `[H]_-^+=0`. Outside the annuli, `[∂wH]_-^+=0`. In the reduced zero-mode problem this becomes the source term

`η_i(k_m h−g_χh s_i)`

in the `h` Euler equation. Thus `h`'s conjugate natural datum is `k_mh−g_χh s_i`, with the displayed sign and normalization.

For the fixed-source ensemble, `Q_χ=s_i` is held fixed during `h` variation and source-strength variation. Translating the body translates `η_i` and is varied when computing `F_var`. Modulation of `J_m` by dynamical `δr_B`, `h`, pressure, or the neighbor is set to zero at G0; only the body-attached orientation projection selects its sign. This is the chosen **one-way fixed-source character** of `χ↔h` and is a load-bearing postulate, not a neutral simplification.

Eliminating the linear quadratic response gives the fixed-source character `Ω_on-shell∝g_χh²·O_h⁻¹`, i.e. inverse stiffness scaling. This identifies the ensemble without assigning the sign of a two-body force.

The later fixed-defect diagnostic exchanges this Legendre member for `h|_{A_i}=s_i h_def η_i/⟨η_i²⟩` with the conjugate reduced source (dimension `E`) left natural; `h_def` has dimension `1` and its stored energy scales directly with stiffness. It is not part of the two E2/E3 runs on this card and may not be selected after seeing their sign.

## 4. Shared equations and complete common interface/collar conditions

For either branch, away from source supports the stationary reduced scalar equations are

```text
−B_eff ∇²u_L − C_hu ∇²h = 0,
−K_h   ∇²h   − C_hu ∇²u_L + Σ_i η_i(k_mh−g_χh s_i) = 0.
```

The time-dependent quadratic operator is

```text
Q_s(ω,k) = [[A_eff ω²−B_eff k², −C_hu k²],
            [−C_hu k²,          M_h ω²−K_h k²]].
```

Common conditions, unchanged between E2 and E3, are:

1. `[H]=0` and the mouth jump of §3.2 at `w=0`; elsewhere the `H` conormal is continuous. `[POSTULATE]`, `H` dimension `L⁻¹`; provenance: the chosen sheet functional.
2. `r_B` is continuous across every diffuse wall piece and obeys the `S_hold` conormal law in §2.2. `[POSTULATE]`, dimension `1`; provenance: one amplitude field and a fixed level set.
3. `u_L` and its conormal `B_eff∂_nu_L+C_hu∂_nh` are continuous across artificial mesh interfaces; `h` and `K_h∂_nh+C_hu∂_nu_L` are likewise continuous except at the stated mouth source. `[CONVENTION]`, field dimensions `L,1`; provenance: natural conditions of `S_Lh`.
4. At `C_i`, the lumen and ambient traces of `ρ,φ` (dimensions `L⁻⁴,1`) are continuous, their mass fluxes (dimension `L⁻³T⁻¹`) obey Kirchhoff conservation, and there is no collar two-surface energy or source. `[POSTULATE]`; provenance: minimal open-mouth matching.
5. `κ_bend=κ_anchor=` collar two-surface tension `=0`, with dimensions `E/L`, `E/L³`, `E/L²`; no unlisted collar jump is permitted. `[POSTULATE]`; provenance: fixed-geometry minimality.

## 5. Phase and collective zero-mode quotient

The fields are quotiented as follows.

1. Bulk phase: `φ∼φ+constant`; impose `∫_{Ω_IR}W_IR φ d⁴y=0`, where `W_IR=C_IR·1_{3R_out/4<|y|<7R_out/8}` and `C_IR` is fixed by `∫W_IRd⁴y=1`. `[CONVENTION]`, `W_IR` dimension `L⁻⁴`; provenance: an explicit positive quotient weight outside every throat and return source that removes only the carrier-phase constant mode.
2. Frozen-body translations: define `W_i=c_Wi exp(−d_Σ,i²/ℓ²)·1_(|y−(X_i,0)|<2L)`, with `c_Wi` fixed by `∫W_id⁴y=1`, and impose `∫d⁴y W_i δr_B ∂_{x_a}r_Σ,i=0` for `a=1,2,3`. This defines `X_i` and removes the translational Goldstone from the field solve. `[CONVENTION]`, `W_i` dimension `L⁻⁴` and integral dimension `L⁻¹`; provenance: an explicit nondegenerate wall-supported version of U1's moduli-fixing requirement.
3. `h`: static decay at the IR boundary fixes its constant mode. `u_L`: `u_L→0` at the IR boundary fixes rigid longitudinal translation. `[CONVENTION]`, field dimensions `1,L`; provenance: unique elliptic inverse.
4. `θ_B`: absent as an independent variable, so there is no brane-phase zero mode to quotient and no double count. `[POSTULATE]`, dimension `1`; provenance: §2.4 representation choice.

The quotient removes only exact zero modes. It does not delete the localized massless `h` pole or the physical longitudinal mode.

## 6. ACTIVE drain and global return

### 6.1 Spatial source laws

Define the compact normalized polynomial bump

```text
B_n(q;σ) = Γ(n/2+3)/[2π^(n/2)σ^n] · (1−|q|²/σ²)² · 1_(|q|<σ),
D_i(y)   = B_3(x−X_i;a/4) B_1(w−s_i a/2;a/4),
∫d⁴yD_i = 1.
```

The global return kernel is

```text
R_0(y)=η_ret(|x−X_ret|)·½[B_1(w−2L;a)+B_1(w+2L;a)],
η_ret(r)=3·1_[R_ret,R_ret+a]/[4π((R_ret+a)³−R_ret³)],
∫d⁴yR_0=1.
```

`X_ret` is the fixed center of the outer return apparatus (the pair barycenter in a held two-body run), and the kernel is exactly `R_w` even.

The wall-to-drain map is the even frozen-profile projection

```text
N_d = ∫d⁴y D_+(y)[r_Σ,+²−r_bg²],
Q_d[r_Σ,s_i] = N_d⁻¹∫d⁴y D_i(y)[r_Σ,s_i²−r_bg²] = 1.
```

Let `Γ_i=Γ_0Q_d=Γ_0>0` be fixed throughput. The mass sources are

`S_drain=−Σ_iΓ_0D_i`, `S_leakage=(Σ_iΓ_0)R_0`.

They have dimension `L⁻⁴T⁻¹=(-4,-1,0)`, are `[POSTULATE]`, and have provenance “minimal smooth finite-size sink plus one two-sided global return.” Their integrals cancel exactly. `g_χdrain:=Γ_0` has dimension `T⁻¹`, star value `10⁻³`, tag `[POSTULATE]`; provenance: smallest explicitly nonzero weak-drain rate, coupled to the frozen sleeve amplitude through `Q_d` but independent of the electric-source normalization. Dynamical `δQ_d` feedback is zero in this fixed-`Σ` card, as recorded in §9.

### 6.2 Momentum and energy sources—all six source slots fixed

Define the carrier enthalpy per particle

`e_c=(ε+P)/ρ`, dimension `E`, on the positive-density admissible branch. First form the comoving source pieces

```text
f_drain    = m S_drain v,
f_leak,0   = m S_leakage v,
w_drain    = e_c S_drain,
w_leak,0   = e_c S_leakage.
```

Then define the return-controller totals

```text
I_ret = −∫d⁴y (f_drain+f_leak,0),
P_ret = −∫d⁴y (w_drain+w_leak,0),
f_leakage = f_leak,0 + R_0 I_ret,
w_leakage = w_leak,0 + R_0 P_ret.
```

Thus

`∫(S_drain+S_leakage)=0`, `∫(f_drain+f_leakage)=0`, and `∫(w_drain+w_leakage)=0`

identically. `f` has dimension momentum/(four-volume·time) `(-3,-2,1)`; `w` has dimension energy/(four-volume·time) `(-2,-3,1)`. All four momentum/energy functionals are `[POSTULATE]`; provenance: comoving removal/injection plus the unique constant finite-rank correction that closes the globally declared return ledger. The correction is applied only on `R_0`, outside every admissible body control surface.

`I_ret` has force dimension `(1,-2,1)` and `P_ret` has power dimension `(2,-3,1)`; both are `[POSTULATE]` auxiliary controller variables with provenance “exact integrated momentum/energy closure,” not tunable response coefficients.

**Return-response choice, flagged:** `δΓ_i/δh_inc=δΓ_i/δP_inc=0` and `δR_0/δh_inc=δR_0/δP_inc=0`. The return **does not constitutively change its mass rate, location, or partition in response to a neighbor's incident `h` or pressure field**. `I_ret` copies removed momentum, and `P_ret` algebraically copies removed energy (so the energy ledger contains the pressure already present in `e_c`); neither changes `Γ_i` or creates a direct `h/P→F_flux` susceptibility. This nonresponsive fixed-throughput choice is `[POSTULATE]`, dimension `1`; provenance: the minimal S2 fixed-`ṁ` closure. A responsive return can change the orientation-pair flux coefficient and is a required alternative in REVIEW NOTES.

The complete weak balance system is

```text
∂tρ+∇·(ρv) = S_drain+S_leakage,
∂t(mρv)+∇·[mρv⊗v+P I+Π_Q] = f_drain+f_leakage,
∂tε+∇·j_ε = w_drain+w_leakage.
```

The source cells `supp(D_i)` and `supp(R_0)` are declared non-Hamiltonian controller/core cells on which these weak balances, rather than an extra variation of `S_cons`, are imposed. The exterior carrier is phase-potential and obeys the conservative Madelung equations. This interface assignment is `[POSTULATE]`, balance dimensions as above; provenance: an explicit open-system closure without pretending the sink or return pump is conservative.

### 6.3 Control surface and source/force partition

For body `i`, use the comoving four-dimensional control surface `C_i(R_c): |y-(X_i,0)|=R_c`, with `R_c=4a`. It intersects the brane transversely on its induced two-sphere; both bulk and brane stress/flux on that cut are included. It may be deformed smoothly away from a collar or source support, but never so as to omit a field channel. It encloses `Σ_i` and `D_i`, excludes every other `D_j` and excludes `R_0`. Any homologous surface in the source-free window `2L<R_c<R/2` is admissible.

The partition convention is:

```text
F_var  = −δΩ_cons/δX_i at fixed {Γ_i,J_m,Q_χ},
         including all passive wall/H/u/bulk traction and S_hold reaction;
F_flux = −∫d⁴y f_drain,i
         (equivalently the active-source momentum read-out on C_i);
F_𝔅    = multiplier/Rayleigh reaction not generated by S_cons;
F_rad  = outgoing wave momentum at the radiative outer boundary.
```

This assignment is `[CONVENTION]`, force dimension `(1,-2,1)`; provenance: U1's open-body four-channel law. The remote `f_leakage` reaction is reported separately in the global momentum ledger and is not silently assigned to either body. Expanding a surface through the return support requires adding that return share; within the declared source-free family the sum `F_var+F_flux+F_𝔅` is surface-independent by the integrated Noether/momentum identity. Any `∝V̇` intake term is kept in `F_flux` and excluded from the `L_eff` inertia.

The flux channel is live because `Γ_0=10⁻³≠0`; setting it to zero is the named drain ablation, not the card.

## 7. IR closure and nondimensionalization

### 7.1 Two-sided return/outer conditions

Use `Ω_IR=B³_{R_out}×[−R_out,R_out]`, with the two-sided fields related by the declared `R_w` action. At the outer boundary:

```text
ρ=ρ_0,       n·∇φ=0 plus the phase quotient,
r_B=r_bg(w), H=0 (static), u_L=0, h=0.
```

For a moving check, replace `H=0` by the outgoing condition for `O_⊥−(M₄/K₄)∂t²`; the uniform-translation frequency is `ω=k·V`, and the card's stationary solution is its subsonic `V→0` limit. The static order of limits is `R_out→∞`, then `R_ret→∞` with `R/R_ret→0` and fixed `Γ_0`. This entry is `[POSTULATE]`; boundary-data dimensions are `ρ:L⁻⁴`, `φ:1`, `r_B:1`, `H:L⁻¹`, `u_L:L`, `h:1`; provenance: minimal two-sided, zero-net-source return with no static absorbing layer.

For this check, every body-attached object is translated as `X_i(t)=X_i(0)+V_it`: `Σ_i`, `η_i`, `D_i`, and the fixed source `Q_χ=s_i`; the return apparatus remains in the medium-rest frame. The same six source functionals and coefficients are used. “Steady” values in this card mean the continuous `|V_i|<min(c_s,c_E)` outgoing solution followed by `V_i→0`, not an independently frozen snapshot. This is `[CONVENTION]`, velocity dimension `(1,-1,0)`; provenance: scope §5's V→0 guard.

### 7.2 Units, values, independent groups, and hierarchy

Take

`L_0=a`, `T_0=a/c_E`, `E_0=A_eff c_E²a³`.

The wall coefficient units are `Z_0=E_0T_0²/a⁴`, `κ_0=E_0/a²`, `λ_0=E_0/a⁴`; brane coefficients use the dimensions in the tables. The independent dimensionless groups fixed by this card are

```text
geometry:  L/a=37/20, ℓ/a=1/20, W_B/a=1/5,
scalar:    B_eff/(A_eff c_E²)=2,
           M_h/(A_eff a²)=1,
           K_h/(A_eff c_E²a²)=1,
           C_hu/(A_eff c_E²a)=1/2,
phase:     ρ_br/A_eff=3/4, (C_J²/κ_phase)/A_eff=1/4,
wall:      Z_χc_E²a²/E_0=1, κ_χa²/E_0=1, λ_χa⁴/E_0=800,
mouth:     k_m/E_0=1, g_χh/E_0=1,
drain:     Γ_0a/(ρ_0a⁴c_E)=10⁻³,
bulk:      c_s/c_E=2, ξ_Q/a=1/20, ρ_0a⁴=1, m c_E²/E_0=1,
           ℏ/(E_0T_0)=sqrt(2)/10, K/(E_0a^16)=4/5,
IR:        R_c/a=4, R_min/a=20, R_ret/a=100, R_out/a=200.
```

Here `ξ_Q:=ℏ/(√2mc_s)`. The four independent bulk ratios `{c_s/c_E,ξ_Q/a,ρ_0a⁴,mc_E²/E_0}` are fixed witness values; the displayed `ℏ`, `K`, and `μ_∞/E_0=1` then follow from `ξ_Q`, `c_s²=5Kρ_0⁴/m`, and `μ_∞=(5K/4)ρ_0⁴`. No relation ties these bulk groups to `A_eff,M_h,K_h,C_hu`, `c_γ`, or the mouth strength. Their values are `[POSTULATE]`, dimension `1`; provenance: a positive finite-density, resolved-quantum, subsonic witness point. The choice `c_s/c_E=2` expressly avoids imposing a cone lock.

The working hierarchy is

`ℓ=0.05a < W_B=0.2a < a < L=1.85a < R_c=4a ≪ R (R≥20a) ≪ R_ret=100a < R_out=200a`.

Pair results must be extrapolated in `a/R`, `R/R_ret`, and `R_ret/R_out`; the quoted finite numbers are a reproducible first grid, not claimed asymptotic equalities.

## 8. The only branch slot: `𝔅`

All common mouth, collar, `H`, `r_B`, `u_L`, source, return, and IR conditions remain exactly those above.

### 8.1 E3 lead — permeable bulk-phase texture

On `Σ_i`, treat the wall as a diffuse internal texture, not a material boundary for the bulk carrier:

```text
[ρ]_-^+ = 0,       [φ]_-^+ = 0,
[ρ(v−V_i)·n]_-^+ = 0,
[Π n]_-^+ = 0.
```

There is no permeability resistance, surface mass storage, phase jump, or surface Rayleigh term. The “phase texture” is the smooth **bulk `θ`** field. The local drain-source match for any pillbox enclosing the lumen/sink is

`∮_{∂P}ρ(v−V_i)·n d³A=∫_P S_drain,i d⁴y`;

in particular, lateral/cap crossing plus mouth flux equals `−Γ_0` once the full sink is enclosed. This is the required E3 drain-source matching, with the sign inherited from `S_drain<0`.

This completion is `[POSTULATE]`, normal mass-flux dimension `L⁻³T⁻¹`; provenance: the zero-resistance endpoint of the permeability axis. Alternatives are finite Robin permeability or a brane-`θ_B` texture; both add structure and are not this branch.

**Why `F_𝔅=0`:** E3 imposes no multiplier constraint and has no Rayleigh functional. Its smooth conservative interface traction is already in `δS_cons/δX_i`, while its non-variational drain impulse is in `F_flux`. Therefore, by the declared channel definition, `F_𝔅^E3=0`—this follows from the empty multiplier/Rayleigh slot, not from ignoring traction.

### 8.2 E2 contrast — impermeable, action-derived free slip

Treat the two faces of `Σ_i` as a holonomic material wall for the bulk carrier while leaving the mouth open:

```text
(v_±−V_i)·n = 0,                  essential no-penetration,
n·∇ρ_± = 0,                       natural quantum-density condition,
P_∥(Π_±n)=0,  P_∥=I−n⊗n,         natural free-slip traction.
```

The first condition restricts admissible phase variations; the latter two follow from varying `S_bulk` and the surface location with unrestricted tangential slip. No tangential velocity is imposed. The lumen drain is matched only through the open mouth:

Explicitly, the E2 surface variation of the background-subtracted bulk functional contains

```text
δΩ_bulk|_Σ = ∫_Σ d³A [
  mρ(v−V_i)·n δΦ
  +(ℏ²/4mρ)(n·∇ρ)δρ
  +(Πn)·δξ ],
```

where `v=∇Φ`, normal displacement is excluded by the holonomic wall, and tangential `δξ` is arbitrary. The three displayed coefficients therefore give no penetration, `n·∇ρ=0`, and `P_∥Πn=0`, respectively. This is the promised action derivation of free slip; it does not introduce an extra traction law.

`∫_{M_i}ρ(v−V_i)·n_M d³A=∫_lumen S_drain,i d⁴y=−Γ_0`,

where `n_M` points from the lumen toward the ambient, so physical inflow into the lumen is negative by convention.

This completion is `[POSTULATE]`, normal mass-flux dimension `L⁻³T⁻¹`; provenance: the minimal inviscid impermeable/free-slip endpoint. There is no no-slip layer, multiplier, or phenomenological traction.

**Why `F_𝔅=0`:** no-penetration is a holonomic restriction inside the variational problem, and free-slip traction is its action-derived natural condition. Their passive normal traction and the common frozen-shape reaction are both counted once in `F_var`. There is no nonholonomic multiplier and no Rayleigh term, so `F_𝔅^E2=0`. Reclassifying this passive traction as `F_𝔅` would violate the card's stated partition.

### 8.3 Stationary radiation channel

For both E2 and E3, the stationary solution has `∂th=∂tu_L=0`, `∂tρ=0`, and only the uniform chemical-potential rotation in `θ`. The scalar wave energy currents are proportional to time derivatives, and the static outer condition has no outgoing `1/R` wave flux. The steady carrier throughput is assigned to `F_flux`, not to radiation. Hence `F_rad=0` for the stationary E2/E3 solutions. The outgoing slightly-moving problem gives no radiation for uniform subsonic translation and tends continuously to this result as `V→0`.

Thus the passive-branch observable made computable by this card is

`F_total=F_var+F_flux`, with `F_𝔅=F_rad=0`.

Neither `F_var` nor `F_flux` may be dropped; only their total (with the remote return ledger when testing global momentum) is physically invariant.

## 9. Complete declared-zero ledger

Within the retained fields, through quadratic order and at most two derivatives, every symmetry-allowed coupling not displayed above is fixed to zero:

| Declared-zero entry | Tag, dimension, provenance |
|---|---|
| bulk `r_BH`, `r_B²H²`, `Hρ`, `Hδρ`, `H∂tθ`, and `H∇θ` couplings outside `Ω_mouth` | `[POSTULATE]`; coefficient dimensions in listed order `E/L³=(-1,-2,1)`, `E/L²=(0,-2,1)`, `E·L=(3,-2,1)`, `E·L`, `(-1,-1,1)`, `(0,-2,1)`; fixed-source mouth is the sole `χ↔h` mechanism |
| dynamical modulation `δJ_m/δr_B`, `δJ_m/δh`, `δJ_m/δP`, neighbor-source response | `[POSTULATE]`; response dimensions `(3,-2,1)`, `(3,-2,1)`, `(5,0,0)`, `1`; fixed-source means the normalized body source is held fixed |
| `r_Bu_L`, `r_B divu`, `r_Bu_T`, `Hu_T`, `u_Lu_T`, and two-gradient scalar–transverse mixing | `[POSTULATE]`; coefficient dimensions `E/L⁴=(-2,-2,1)`, `E/L³=(-1,-2,1)`, `E/L⁴`, `E/L³`, `E/L⁵=(-3,-2,1)`, `E/L³`; electric MVP retains only `C_hu` scalar mixing |
| cross-kinetic `∂tu_L∂th`, one-time-derivative/Berry `u_L∂th` and `h∂tu_L` | `[POSTULATE]`; coefficient dimensions `M/L²=(-2,0,1)`, `(-2,-1,1)`, `(-2,-1,1)`; time-reversal-even minimal scalar block |
| `u_L²`, reduced `h²`, `(∇²u_L)²`, `(∇²h)²`, `h³,h⁴`, and `(∇u_L)³` | `[POSTULATE]`; coefficient dimensions `E/L⁵=(-3,-2,1)`, `E/L³=(-1,-2,1)`, `E/L=(1,-2,1)`, `E·L=(3,-2,1)`, `E/L³`, `E/L³`; keep both far-field channels gapless and minimal |
| independent primitive `B(divu)²` in addition to `B_eff` | `[POSTULATE]`, `B=(-1,-2,1)`, value `0`; prevents double counting the density-induced modulus |
| independent `θ_B`, its amplitude-weighted kinetic/gradient terms, Josephson `θ−θ_B`, and brane-phase drain | `[POSTULATE]`; representative coefficient/source dimensions `(-2,0,1)`, `(0,-2,1)`, `E/L³=(-1,-2,1)`, `L⁻³T⁻¹=(-3,-1,0)`; eliminated-phase representation, with its effect retained only in `A_eff` |
| wall bending, anchoring, collar two-surface tension, surface number storage, and surface dissipation | `[POSTULATE]`; dimensions `E/L=(1,-2,1)`, `E/L³=(-1,-2,1)`, `E/L²=(0,-2,1)`, `L⁻³=(-3,0,0)`, `E/(L³T)=(-1,-3,1)`; fixed-`Σ` passive MVP |
| dynamic `δQ_d/δr_B`, `Γ(h_inc)`, `Γ(P_inc)`, return-location responses to `h,P`, return-kernel responses to `h,P`, and drain-rate orientation dependence | `[POSTULATE]`; response dimensions `1`, `T⁻¹=(0,-1,0)`, `(2,1,-1)`, `L=(1,0,0)`, `(3,2,-1)`, `L⁻⁴=(-4,0,0)`, `E⁻¹=(-2,2,-1)`, `T⁻¹`; fixed-throughput, charge-even mass drain |
| direct drain sources in the `h` and `u_L` equations, and direct `h` contribution to `e_c` | `[POSTULATE]`; source/coefficient dimensions `E/L³=(-1,-2,1)`, `E/L⁴=(-2,-2,1)`, `E=(2,-2,1)`; drain acts on the one bulk carrier only |
| field-dependent geon derivatives `δE_g/δ{r_B,h,θ,H,u_L,ρ}` | `[POSTULATE]`; dimensions `{E,E,E,E·L,E/L,E·L⁴}`; the held-out geon rest constant cannot source the electric MVP |
| bulk viscosity, tangential drag, E2 no-slip traction, E3 permeability resistance, E3 phase jump | `[POSTULATE]`; dimensions `(-2,-1,1)`, `(-3,-1,1)`, `E/L⁴=(-2,-2,1)`, `(5,-1,1)`, `1`; inviscid endpoint definitions |
| all E4 multipliers, E5 Rayleigh kernels, E1 reactions, and mixture terms | `[POSTULATE]`; reaction-traction dimension `E/L⁴=(-2,-2,1)` and Rayleigh surface-power dimension `E/(L³T)=(-1,-3,1)`; outside the two registered branches |
| Maxwell/gauge fields, point sources, native current law, Coulomb potential prior | `[CONVENTION]`, dimension `∅` because no such field or coefficient exists; prohibited ancestry |

This is an exhaustive zero ledger for the card's stated derivative/field truncation. Higher-order operators are not “implicitly small”; they are zero in G0-v0 and require a new version to activate.

## 10. Instantiability, gates, and checks

### 10.1 Elliptic problem supplied by the card

For either `𝔅`, solve the wall-amplitude Euler equation with its fixed level set, the sourced stationary Madelung continuity/Bernoulli system, the parent `H` equation (or its normalized zero-mode reduction), and the coupled `(u_L,h)` equations, subject to §§3–8. The permanent throat source plus an incident long-wavelength `(h,u_L,ρ,φ)` field defines an affine one-body response map. Its multipoles, local drain impulse, passive Noether force, and return ledger are separate read-outs. No reciprocity is assumed for the active source.

The two-body response is therefore computable in principle by connecting two affine one-body maps with the card's Green operators **and** solving the explicitly displayed finite-rank `{I_ret,P_ret}` return constraint. If that shared constraint prevents one-body factorization, the response-first pass records the specified analytic no-go and escalates to the two-center solve; it is not a missing closure. This statement is about poseability, not about the sign or analytic factorization succeeding.

### 10.2 Exclusion gates wired to the solve

The committed `g0_closure_card_v0_checks.{py,wl}` suite verifies **static algebraic prerequisites only**. Before any force language, require the following gates, with their present dispositions:

1. **Class (1), verified:** `PASS_LOCALIZED_H_NORM` and `PASS_LOCALIZED_H_PHYSICAL_NORM` establish a normalizable localized `H` zero mode with `N₀>0` and `M₄N₀>0`; delocalized support is `FAIL_DELOCALIZED_BULK_1_OVER_R3`.
2. **Class (1), verified:** `PASS_REDUCED_H_INERTIA`, `PASS_STABILITY`, `PASS_TRANSVERSE_FACTORIZATION`, `PASS_PLANAR_WALL_FACTORIZATION`, `PASS_BULK_U_SECOND_DERIVATIVE`, `PASS_BULK_FOURIER_DENSITY_HESSIAN`, and `PASS_PHASE_CONSTANT_MODE_QUOTIENT` establish `M_h>0`, positive reduced scalar stiffness, the nonnegative parent transverse operator, the nonnegative planar unconstrained wall Hessian, the positive homogeneous bulk-density Fourier Hessian, and the explicit bulk-phase constant quotient. **Class (2), deferred:** the curved/held sleeve-slab claim still requires the lowest eigenvalue of the assembled constrained one-body `r_B` Hessian with `S_hold`, boundary conditions, and translation quotient; the inhomogeneous quantum/bulk claim still requires the lowest eigenvalue of the assembled one-body `(δρ,δφ)` Hessian about the sourced solution with its boundary conditions and quotient. Either negative eigenvalue is `FAIL_GHOST_INSTABILITY`.
3. **Class (1), verified:** `PASS_BARE_MOUTH_SOURCE` establishes a nonzero bare projected source from the independently fixed mouth functional. **Class (3), deferred:** `FAIL_NO_MONOPOLE` cannot be cleared until the solved dressed two-body/far-field `h` response supplies a nonzero `1/R` monopole coefficient `q_{h,dressed}`.
4. **Class (1), verified:** `PASS_REDUCED_H_MASSLESSNESS` and `PASS_PARENT_TRANSVERSE_ZERO_EIGENVALUE` establish that the reduced `h` block has no `k⁰h²` term, `Q_s(0,0)=0`, and the parent transverse operator has a zero eigenvalue. These are the static prerequisites for, not a solved assertion of, the assembled Green mode; failure is `FAIL_PINNED_BRANON` or `FAIL_YUKAWA`.
5. **Class (3), deferred:** the `(1,0)` and `(0,1)` Hadamard claims require the solved four-orientation, channel-resolved pair force and its `H_10,H_01` projections (equivalently the solved two-body derivatives `F^(1,0),F^(0,1)`) on the `R_w` background before either cell can be declared zero or ambient contamination.
6. **Class (1), verified:** `PASS_DRAIN_KERNEL_NORMALIZATIONS`, `PASS_DRAIN_MASS_CONTROLLER`, `PASS_DRAIN_MOMENTUM_CONTROLLER`, and `PASS_DRAIN_ENERGY_CONTROLLER` establish normalized `D_i,R_0`, positive `Γ_0`, and the global controller identities. **Class (2), deferred:** one-body finite-volume closure requires the assembled solved `ρ_i,φ_i`, mass flux through `C_i`, integrated stress/source momentum residual, and solved `ε_i,j_{ε,i}` energy residual; the drain ablation requires paired `Γ_0` and `Γ_0=0` assembled one-body solutions and the `F_flux,i` readout. **Class (3), deferred:** pair mass/momentum/energy closure requires `R_ρ^(2),R_p^(2),R_ε^(2)` on the shared-controller two-center solution, and the outer-return ablation requires nominal and return-ablated two-center solutions plus the remote-return momentum/energy ledger.
7. **Class (1), verified:** `PASS_CONSERVATIVE_HESSIAN_SYMMETRY` establishes the action-derived symmetric-Hessian prerequisite. **Class (3), deferred:** pair-force integrability still requires the configuration-space curl/cross-derivatives of the solved `F_var,11^(2)`, and the physical total still requires solved channel readouts, including `F_var` and `F_flux`, before any `U_11` is written.

### 10.3 `pathA_38` regression—barred from calibration

Only after the preceding parent operator, `f₀`, normalization, mouth support, source sign, and coefficients are frozen, the zero-mode facts `N₀=8/(3ℓ)`, three-dimensional `1/R` static Green behavior, and relative `±w` source projection can be compared with `pathA_38_throat_body_electric_localization.md`. Agreement is a pipeline regression; disagreement is reported. `J_m`, `g_χh`, `K_m`, and their sign were **not** fitted to reproduce that report's `U₊₊>0`, and that variational result cannot determine the full physical sign.

### 10.4 Committed static verification record

The committed `g0_closure_card_v0_checks.{py,wl}` suite verifies **static algebraic prerequisites only**; it does not discharge class-(2) assembled one-body claims or class-(3) two-body/far-field claims.

The card values give

```text
M = Matrix([[A_eff,0],[0,M_h]]) = Matrix([[1,0],[0,1]])
K = Matrix([[B_eff,C_hu],[C_hu,K_h]]) = Matrix([[2,1/2],[1/2,1]])
M_h = 1 > 0
B_eff*K_h − C_hu² = 2*1 − (1/2)² = 7/4 > 0
det(K−c²M) = c⁴−3c²+7/4
c_±² = (3±sqrt(2))/2 = {0.7928932188..., 2.2071067812...}
eig(K) = {(3−sqrt(2))/2,(3+sqrt(2))/2} > 0
O_⊥ f₀ = 0
O_⊥ = A†A >= 0
L_χ^(2) r₀' = 0 for r₀=[1+exp(−x/ℓ)]⁻¹
L_χ^(2)/κ_χ = A_χ†A_χ,
  A_χ=∂x+tanh[x/(2ℓ)]/ℓ
bulk witness: K*=4/5, ℏ*=sqrt(2)/10, μ_∞*=1, U''(ρ_0)*=4>0
```

The explicit SymPy check also verifies term dimensions:

```text
wall kinetic = wall gradient = wall potential = E/L⁴,
parent-H kinetic = parent-H gradient = parent-H potential = E/L⁴,
u kinetic = u stiffness = h kinetic = h stiffness = C_hu mix = E/L³,
mouth Robin = mouth source = E after integration,
S = 1/(L⁴T), f = momentum/(L⁴T), w = E/(L⁴T).
```

Class-(1) result, verified by `g0_closure_card_v0_checks.{py,wl}`: `PASS_STABILITY`, `PASS_POSITIVE_GENERALIZED_WAVE_SPEEDS`, `PASS_TRANSVERSE_FACTORIZATION`, `PASS_PLANAR_WALL_FACTORIZATION`, `PASS_POSITIVE_BULK_HESSIAN`, and `PASS_DIMENSIONAL_HOMOGENEITY`. `PASS_PLANAR_WALL_FACTORIZATION` covers only the planar unconstrained kink; full curved/held sleeve-slab wall positivity is a class-(2) claim deferred to the lowest eigenvalue of the assembled constrained one-body wall Hessian with `S_hold`, boundary conditions, and translation quotient.

For the nullspace and Fredholm statements, the scripts verify only the static algebraic prerequisites just listed, together with the explicit phase quotient and reduced `h` masslessness. The full assembled claims are deferred as follows. **Class (2):** identifying the bulk-phase, constrained-wall, coupled `(H,u_L,h)`, and longitudinal kernel bases; proving that the phase, hold/moduli, and outer-boundary conditions remove the unwanted assembled nulls; and establishing compatibility of `{I_ret,P_ret}` with the one-body operator all require the respective restricted kernels/inverses, the assembled constraint Jacobian/Schur block, and the principal symbol of the controller-augmented system. The finite-`R_out` elliptic/Fredholm and nullspace-fixed claims additionally require the assembled operator domain and boundary map, Fredholm index, kernel/cokernel bases, and one-body compatibility residuals. **Class (3):** preservation of the physical `h` far-field Green mode requires the constrained zero-frequency far-field Green function `G_hh^(2)` and its pole residue; convergence of the IR extrapolation requires the two-center solution family and the `R_out,R_ret→∞`, `R/R_ret→0` limit of the pair observables.

## REVIEW NOTES

This section is deliberately redundant: choices not marked user-ratified remain live for review; the marked choices are frozen for v0 while their listed alternatives remain visible.

### Load-bearing choices and live alternatives

1. **`χ↔h` character:** `[USER-RATIFIED 2026-07-20]` chosen as a linear, orientation-odd, fixed-source Robin sheet functional built from a normalized odd projection of the frozen wall amplitude. Alternatives: a fixed defect/Dirichlet datum; a reciprocal dynamical `r_BH` coupling; `r_B²H²`; derivative or nonlocal coupling; a source on the sleeve rather than mouth; or no monopole. The choice is minimal and finite-size, but it is the electric mechanism itself and could change range, residue, or sign.
2. **Drain/return character:** `[USER-RATIFIED 2026-07-20]` chosen active at fixed nonzero throughput, with an `R_w`-even remote return and exact global mass/momentum/energy correction. The sign-channel constitutive answer is **NO direct neighbor response**: `δΓ/δh_inc=δΓ/δP_inc=0` and the return kernel is fixed. The auxiliary controller nevertheless copies the momentum/energy actually removed—including pressure in the energy ledger—and is solved as a shared finite-rank constraint. Alternatives: `Γ=Γ(h,P)`, a pressure-controlled return, a deformable pair-dependent return kernel, a local recirculation loop, or an externally supported state-blind injector. State-blind return is the v0 base, but the responsive `Γ(h,P)` susceptibility check is a **MANDATORY in-grind sign-gate**: no electric force sign may be promoted until the linear responsive perturbation is computed and shown not to flip sign or range. It is not deferred to a later card.
3. **Phase representation:** chosen `Q_s(u_L,h)` with `θ_B` eliminated and `A_eff=ρ_br+C_J²/κ_phase`, not bare `ρ_br`. Alternatives: retain the regular amplitude-weighted `(u_L,θ_B)` block; use a bulk-to-brane pullback; or deliberately use bare `ρ_br` as a separate comparator. E3 texture is the bulk `θ`, not `θ_B`. Changing representation changes poles and the field count.
4. **`C_hu`:** `[USER-RATIFIED 2026-07-20]` chosen nonzero, `C_hu*=1/2`, with stability margin `7/4`. Alternatives: a preregistered nonzero scan or the explicit `C_hu=0` **ablation**. Zero is not neutral: pathA39 pole residues and scalar admixture depend on it.
5. **Localized `H` operator:** chosen factorable Pöschl–Teller `A†A` with one gapless localized zero mode. Alternatives: a wall-amplitude-dependent localization operator, another reflection-even well, a delocalized four-dimensional scalar, a pinned/Yukawa mode, or no scalar. The analytic match to the old zero-mode profile is only a regression check.
6. **Frozen shape:** chosen cylinder-plus-cap `Σ`, logistic wall width, and a body-attached holonomic mid-surface constraint whose reaction is in `F_var`. Alternatives: finite bend/anchor stiffness or free-`Σ`. A lab-pinned surface is not an allowed reinterpretation.
7. **E3:** chosen zero-resistance bulk-`θ` texture with smooth transmission. Alternatives: finite permeability, phase jump, surface storage, or a brane-`θ_B` texture; those would be different `𝔅` cards.
8. **E2:** chosen no-penetration with action-derived quantum/free-slip natural traction. Alternatives: E1 no-slip, finite slip, or a multiplier implementation. Those alter the channel partition and are not this contrast.
9. **Mouth Robin impedance:** chosen positive unit reduced `k_m`. Alternatives: pure Neumann `k_m=0`, stronger Robin impedance, or the mandatory later fixed-defect exchange. It was not fit to the pathA38 energy sign.
10. **Bulk/source energy convention:** chosen the canonical Madelung stress/current and enthalpy-per-particle source bookkeeping. Alternatives differ by explicit improvement terms or reservoir thermodynamics; they must be versioned and recheck finite-volume balance.
11. **IR return:** chosen remote, two-sided, state-independent geometry with a finite-rank conservation controller. Alternatives include the one-sided pathA29 slab or a genuinely nonlocal pair return. The electric census remains conditional on the `R_w` ambient.
12. **Scalar numerical point:** `A_eff*=M_h*=K_h*=1`, `B_eff*=2`, `C_hu*=1/2`. These are postulated witness values, not fits. A future sensitivity card must vary them before any sign is promoted.
13. **Force partition:** passive action-derived traction and shape holding are in `F_var`; local sink impulse is in `F_flux`; only multiplier/Rayleigh reactions could enter `F_𝔅`. Alternatives can repartition passive traction, but must transform all channels together and leave the total invariant. The displayed convention was chosen to make `F_𝔅=0` demonstrable for both passive branches.
14. **Wall functional:** chosen constant `Z_χ,κ_χ` and the lowest polynomial `r_B²(1−r_B)²` double well. Alternatives include field-dependent stiffness, polar/smectic invariants, or additional slab-stabilizing interactions. The holonomic fixed-shape card is minimal but cannot be promoted to free-`Σ` physics.
15. **`R_w` representation of the mixed scalars:** `[USER-RATIFIED 2026-07-20]` chosen `H,h,u_L` all odd in the retained orientation channel so nonzero constant `C_hu` is symmetry-compatible. Alternatives are a charge-even density `u_L` with `C_hu=0`, an orientation-dependent spurion coefficient, or a larger even/odd longitudinal block. This parity assignment is as load-bearing as the numerical value of `C_hu`.

### Every declared-zero coupling

The complete list is the §9 table: bulk wall–`H` terms outside the mouth; source modulation; wall/longitudinal/transverse cross terms; cross kinetic and Berry terms; scalar masses, nonlinearities, and higher gradients; primitive `B`; independent `θ_B` and Josephson/interface terms; bending/anchoring/line/storage/dissipation terms; return response and orientation-dependent drain rate; direct drain-to-`h/u_L`; field-dependent geon couplings; viscosity/drag/no-slip/finite permeability; E1/E4/E5/mixture reactions; and all Maxwell/gauge/point/current imports. None is “understood nonzero.”

### What this card does not fix

It does **not** fix free-`Σ` dynamics (deferred beyond R1-MVP), E4 or E5 (deferred reduced static solves), the electric/gravity hierarchy magnitude, or magnetism/current structure (R2+). It also does not select the physical mouth ensemble; the fixed-defect exchange remains mandatory before a sign can be promoted.

### Committed static-check result for reviewers

The checked numerical stability margin is `7/4`; the generalized squared wave speeds are `(3−sqrt(2))/2≈0.7928932188` and `(3+sqrt(2))/2≈2.2071067812`; `M_h=1`; the transverse operator annihilates `f₀=sech²(w/ℓ)/ℓ` and factorizes as `A†A`; the planar wall Hessian annihilates `r₀'` and factorizes as `A_χ†A_χ`; the bulk witness gives `U''(ρ_0)*=4`; all action and source terms pass their `(L,T,M)` dimension equalities. The committed `g0_closure_card_v0_checks.{py,wl}` suite verifies these class-(1) static algebraic prerequisites only; class-(2) assembled one-body and class-(3) two-body/far-field claims remain unproved.
