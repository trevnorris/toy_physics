# S11c-a — SHARED PHYSICS (background & geometry: the tilted-face shape-derivative)

**Orchestrator-written; folded once after two legs (Codex + Grok), 2026-08-23.** The physics authority both
engines read for **S11c-a**, the first sub-step of the S11c family (`directives/S11c_decisions.md`,
`N1`/`N2`). ⛔ Not a build directive (per-engine wiring is separate). ⚠ The first draft's defects are folded
below and marked ⭐**[folded]**; each was verified against source by the orchestrator (rule 13).

## 0 · What this sub-step is, and what it is NOT

S11c-a **un-freezes the in-plane background** and computes the **first-order geometry** the rest of the
family needs: the tilted-face kinematics and the shape-derivative of **every** interface law S11b wrote for a
flat wall. ⭐ It exports **geometric boundary operators and shape-derivative residuals**, ⛔ not the
transverse↔thickness coupling itself.

⛔ **NOT in S11c-a** (each a later sub-step, `directives/S11c_decisions.md:N2`): the off-diagonal coupling
kernel and the variable-coefficient energy/constitutive operator + its new invariants (S11c-b, `N15`);
**solving** the outgoing curved-bulk problem and constructing the **nonlocal** DtN/self-energy (S11c-c) — ⭐
S11c-a exports the *geometric* boundary shape-operators (the trace and normal-derivative/conormal), S11c-c
**solves**; any spectrum/scattering/Bloch/WKB object (S11c-d); the leakage observable, the confinement
verdict, the bench-optics falsification (S11c-e). ⛔ **NOT anywhere in S11c**: the nonlinear
DC/harmonic/sideband program (`G14`). ⛔ **No terminal `VERDICT`.**

## 1 · Inherited setup — SUPPLIED, imported from the closed S11b, ⛔ not re-derived

The geometry, fields, conventions and interface laws are S11b's (`directives/S11b_SHARED_PHYSICS.md`); point
at them. The load-bearing objects S11c-a takes the shape-derivative of, verbatim:

- **Fields** (`S11b_SHARED_PHYSICS.md:72-88`): in-plane displacement `u` (**3 components**, ⛔ no
  `w`-component); thickness `δW`, `e_W ≡ δW/W₀`; **Eulerian** densification `θ` (`ρ_4D = ρ_4D⁰(1+θ)`);
  background inertia `ρ_br⁰ ≡ ρ_4D⁰ W₀`; the twist modulus `μ_R`. ⚠ Use the **harmonic-in-time** convention
  only; ⛔ the plane-wave `∝ exp[i(k·x−ωt)]`, `k` real table (`:87`) is for the uniform problem and ⛔ must
  **not** be applied to the background profile of §2 (`N5`).
- **Faces** (flat, in S11b, `:368`): `w = ±W₀/2`, outward normals `n̂_± = ±ŵ`, `ζ_± = ±δW/2`, outward face
  velocity `V_± ≡ (∂_tζ_±)(±1)`.
- **Relative flux** (`:194-201`): `J_± ≡ ρ_m(v_w − ∂_tζ_±)·(±1)`, the mass flux leaving the slab **along the
  face's outward normal** — the outward sign appears **once**.
- **Face closure** (`:212-227`): `J_± = Λ_A(ω)𝒜_± + Λ_V(ω)V_±`, `Λ_I = Λ_I⁰/(1−iωτ_I)`; affinity
  `𝒜_± ≡ μ_s − δp_±/ρ_m`, `μ_s ≡ μ_θ/ρ_br⁰`, `μ_θ ≡ δU/δθ`. ⛔ Use exactly.
- **Virtual constraint** (`:337`, on the uniform background): `δ_vθ + δ_v e_W + ∇_x·δ_v u = 0`, binding from
  the **material** law `δ_v Σ_mat = 0`, `Σ_E ≡ ρ_4D W`, `Σ_mat ≡ Σ_E(x(X,t),t)𝒥_x` (`:320-335`). ⛔ `δ_vΣ_E
  = 0` is **not** the constraint (`:341`).
- **Physical (sourced) mass balance** (`:640`, kinematics — a **separate EVOLUTION** equation, ⛔ not the
  same object as `δ_vΣ_mat=0`, `:356`): `∂_t Σ + ∇_x·(Σ v) = −(J₊ + J₋)`. ⭐ Supplied here because T-h
  linearises it.
- **Virtual work & traction** (`:370-372`): `δ_vx_± = n̂_± δ_v(δW)/2`, `t_± = −(δp_± + Λ_X(ω)𝒜_±)n̂_±`,
  `δ_v𝒲_bulk = Σ_s t_s·δ_vx_s = −½[δp_+ + δp_− + Λ_X(𝒜_+ + 𝒜_−)]δ_v(δW)`. ⭐ Supplied so T-d can
  shape-differentiate it; note on a flat wall `t_∥ = 0` so this reduces to thickness-only work.
- **Bulk face response** (`:593-637`): the projection identity B0a (⭐ including its **dynamic window**
  `Ω(w;x,t)` term enumeration, `:600-606`), and the interfacial mass balance `v_bulk,± = V_± + J_±/ρ_m`
  (`:622`, kinematics).
- **The drain** `v_dr ≡ v_bulk_normal_0` (`:99-111`): the bulk **normal** drain; a scope limit only.
- **Energy & method** (`:255-380`): the stored-energy basis `U` (⭐ its `½κ_W W₀²|∇(δW)|²` gradient term,
  `:260`); balance laws (⛔ not an action principle); the virtual-displacement rule; **variational** (⛔ not
  ordinary-partial) derivatives; the sign conventions; the **exit semantics** (`:518`, `:543`).

⚠ **Everything in §1 is SUPPLIED** and **unfalsifiable within this build** — a passing S11c-a does not
re-verify S11b. What S11c-a **tests** is the shape-derivative of these objects under §2's un-freezing.

## 2 · The non-uniform background — un-freezing, the varying moduli, admissibility, anchoring, power counting

**2a · Un-freeze the wall AND the twist modulus in-plane (`N14`, §5).** Introduce two `O(1)` dimensionless
profiles and reserve **new** names: the background thickness `W_bg(x) = W̄₀(1 + η w₁(x))` and — ⭐**[folded:
this was omitted]** — the **twist modulus profile** `μ_R,bg(x) = μ̄_R(1 + η m₁(x))`. ⭐ `∇μ_R ≠ 0` is the
**principal** transverse↔longitudinal mixing source (`V3_STEP_PLAN.md:1179`); S11c-a **states that μ_R
varies, at order η, with its anchoring (2c)** — the operator it produces is S11c-b's (`N15`), ⛔ but S11c-a
may not leave μ_R constant. ⛔ Do not reuse the imported constant keys `W_0`, `mu_R`, `rho_br`, `e_W`, `v_0`.

**2b · The freeze reconciliation — STATE which density is constant, then COMPUTE the consequence (`N11b`).**
S11b froze `ρ_br⁰ = ρ_4D⁰ W₀` to the constant `rho_br` (`:75`). State the modeling choice — which of
`ρ_4D⁰`, `W_bg`, `ρ_br⁰` is spatially constant — as an explicit premise, then **compute** `∇Σ_E⁰ =
∇(ρ_4D⁰ W_bg)`. ⭐ Report both admissible representatives: (i) `ρ_4D⁰` constant ⇒ `ρ_br⁰(x)` varies; (ii)
`ρ_br⁰` constant ⇒ `ρ_4D⁰(x)` varies, `∇Σ_E⁰ = 0`. Reserve a new name for any varying `ρ_br⁰(x)`.

**2c · Profile anchoring — BOTH branches, separately named (`N4`). ⭐[folded: was engine-local]** The
inhomogeneous background is either **advected with the material** or **held fixed in lab/Eulerian space**;
the two give different `O(εη)` terms and the difference is the channel being sought — ⛔ so the choice may
**not** be engine-local. Supply **both** as separately named physical branches and compute each: the
**lab-held** branch has the profile a fixed function of Eulerian `x` (e.g. `W_bg(x)`); the **material-
advected** branch has it attached to material points `X`, i.e. `W_bg(x−u)` to first order, or the
parametric face map `R_s(X,t)`. ⚠ Eulerian/material *variables* are the same physics after `Δρ = δρ_E +
u·∇ρ⁰`; **anchoring** is the genuine physical input, ⛔ not the representation.

**2d · Admissibility — STATE the premise and name the support; the failing check is computed with S11c-b's
operator. ⭐[folded: the drafted residual could not fail]** ⛔ Do **not** insert `W_bg − W̄₀` as a static
`δW` into S11b's perturbation equations (they are linear about constant `W₀`; on a zero-field background the
residual is identically zero — a control that cannot fail). Instead: **name the background state**
(`θ⁰, p_s⁰, J_s⁰` and any boundary load) and **declare** whether the background is **self-supporting** (a
stationary point of the background-order energy) or **externally held** (name the support-force density). ⭐
The background-order stationary balance uses the **variable-coefficient** energy that N15 assigns to S11c-b;
S11c-a fixes the premise and the named support, and the **residual that can fail is computed in S11c-b**.
⇒ `S11CA_BACKGROUND_STATE`, `S11CA_ADMISSIBILITY_PREMISE`.

**2e · Power counting — CONTRAST vs SLOPE, on EVERY object (`N12`). ⭐[folded: contrast≠slope]** Carry
**three** small parameters, kept distinct: wave amplitude `ε` (`u, δW, θ ∼ O(ε)`); background **contrast**
`η` (`W_bg = W̄₀(1+η w₁)`); and the geometric **slope** `|∇W_bg| ≪ 1` (the small parameter the shape
expansion is in — ⚠ a slit edge, `N7`, has order-unity slope at small or large contrast, so `η` and the
slope are **not** the same). ⭐ **Every emitted object carries its computed order.** The corrections S11c-a
exports are `O(η)`/`O(slope)` (background) and `O(ε·η)`/`O(ε·slope)` (wave × tilt); ⛔ never discard an
`O(εη)` term as "second order."

## 3 · What to compute — the tilted-face shape-derivative of every interface law (⛔ name the object, ⛔ do not pre-state its form)

⭐ **Supply (the geometry only):** each face is the oriented level set with the **outward** unit normal
fixed by orientation, ⛔ not by a bare gradient sign. For faces `F_s(x,w,t) = 0` (`s=±`; lab-held
`F_s = w − ½ s W_bg(x) − ζ_s`; material-advected `F_s = w − ½ s W_bg(x−u) − ζ_s`), the outward unit normal
`n̂_s` is the unit normal with `s(n̂_s·ŵ) > 0` (pointing **out of** the slab `|w|<W_bg/2`). ⭐**[folded: the
draft's `F_±` oriented the lower face inward]**. Compute and emit, **each to its `(ε,η,slope)` order**, for
**both** anchoring branches (2c):

- **T-a · Outward unit normal** `n̂_s`. ⇒ `S11CA_FACE_NORMAL`.
- **T-a′ · Normal-derivative / conormal operator** `n̂_s·∇₄` on the tilted face — the geometric boundary
  operator S11c-c consumes (⛔ do not solve the bulk problem; export the operator). ⇒ `S11CA_CONORMAL_DERIV`.
- **T-b · Outward face velocity** — the face's outward normal velocity (the S11b `V_s=(∂_tζ_s)(±1)` is the
  flat limit); the material branch carries the extra tilt term. ⇒ `S11CA_FACE_VELOCITY`.
- **T-c · Relative flux — the FULL first-order `J_s`.** ⭐**[folded: was a component filter]** Emit
  `J_s = ρ_m(v − v_face)·n̂_s` as an **equation**, with **no** component filter — the in-plane velocity
  (transverse **included**) projects onto the tilted normal at `O(εη)`. ⛔ Do **not** decide which components
  enter by a divergence-free argument (`:821-823`); that in-plane projection is the kinematic image of the
  coupling S11b proved zero on a flat wall. ⇒ `S11CA_RELATIVE_FLUX`.
- **T-c′ · Tilted interfacial kinematic mass balance** — `n̂_s·v_bulk,s = V_{n,s} + J_s/ρ_m`, the permeable
  boundary datum S11c-c inherits (⛔ leaving it as the flat Cartesian `v_w = V + J/ρ_m` is inconsistent at
  `O(η)`). ⇒ `S11CA_KINEMATIC_BALANCE`.
- **T-d · Traction and its FULL virtual work.** ⭐**[folded: the in-plane pairing was missing]** The traction
  `t_s = −(δp_s + Λ_X𝒜_s)n̂_s` with the tilted `n̂_s`, **and** the shape-derivative of `δ_v𝒲_bulk` — which
  now carries, beside the thickness pairing `∝ δ_v(δW)`, the **in-plane** pairing `t_∥·δ_v u` (`t_∥ = O(εη)`
  once `n_∥ = O(η)`), a second `O(εη)` channel. ⇒ `S11CA_TRACTION`, `S11CA_TRACTION_INPLANE_WORK`.
- **T-e · Bulk-field evaluation at the shifted face** — the `O(η)` Taylor shift from `w=±W̄₀/2` to the
  tilted `w=±W_bg(x)/2`. ⇒ `S11CA_FACE_SHIFT`.
- **T-f · Projection current with the DYNAMIC window** `Ω(w;x,t)`. ⭐**[folded: a static `Ω(w;x)` kills the
  `∂_tΩ` term]** Report every term present then and absent from the flat/static case (`:600-606`), including
  `∂_tΩ` (an `O(εη)` term for a material-following window); state the window's anchoring (2c). ⇒
  `S11CA_PROJECTION_EXTRA_TERMS`.
- **T-g · The VIRTUAL constraint on the non-uniform background.** ⭐**[folded: the draft pre-stated a
  dimensionally-illegal `u·∇Σ_E⁰` and used the physical `u`]** Linearise the **material** law
  `δ_v Σ_mat = 0` (properly normalised — ⛔ do not add a `Σ`-dimension term to the dimensionless flat
  constraint) about the non-uniform background, under each anchoring, and report **all** `O(εη)` terms it
  produces — ⛔ do **not** pre-state any. ⚠ The pullback samples the **virtual** displacement `δ_v u`, ⛔ not
  the physical `u`; and unfreezing `W₀→W_bg(x)` re-touches `e_W ≡ δW/W₀` (a `W̄₀` vs `W_bg(x)` factor that
  survives **even** where `∇Σ_E⁰=0`) — report it as its own term or introduce a fresh `e_{W,bg}` with an
  explicit map to imported `e_W`. ⇒ `S11CA_VIRTUAL_CONSTRAINT`.
- **T-h · The PHYSICAL (sourced) mass balance.** ⭐**[folded: was missing]** Linearise `∂_t Σ + ∇_x·(Σ v) =
  −(J₊+J₋)` (`:640`) about the **non-uniform** background — a **separate** object from T-g's virtual
  constraint (`:356`), carrying the flux source. ⇒ `S11CA_EVOLUTION_MASS_BALANCE`.
- **T-i · The assembled-closure shape-derivative.** ⭐**[folded: was missing]** The shape-derivative of the
  assembled face law `J_s − Λ_A𝒜_s − Λ_V V_s = 0` about the non-uniform background — including the
  **affinity-normalisation** `O(εη)` term from a varying `ρ_br⁰(x)` in `μ_s = μ_θ/ρ_br⁰` (representative (i)),
  a slab-side term T-e's bulk-field shift does not generate. ⇒ `S11CA_CLOSURE_SHAPE_DERIV`.

⭐ Each `S11CA_*` reduces to its S11b object as `(η, slope) → 0` — ⚠ report that reduction as a **regression
smoke test** (`N6`), ⛔ not the primary control (`η→0` is known-vacuous for the coupling; §8 is the real
control).

## 4 · The representation-invariance control (`N4`/`N6`) — two routes, EMIT the residual, one-sided corruption

Derive T-g (and one of T-c/T-i) **twice** — once in **Eulerian** variables, once after flattening the faces
into material coordinates (`x = X + u`, `w′ = w/W_bg`) — **using the same explicit anchoring map (2c) on both
routes**. ⭐ **Emit both routes' operators AND their symbolic difference** (print, then guard); ⛔ do not
assert what the difference equals (rule 5) — a nonzero residual is the finding. ⚠ **Independence by
one-sided corruption:** break **one** route only (omit its advective-density contribution, or flip one
face's slope term), and report that the residual moves; ⛔ a zero residual proves nothing if the second
route is derived from the first. ⇒ `S11CA_REP_INVARIANCE_RESIDUAL`, `S11CA_REP_INVARIANCE_CORRUPTION`.

## 5 · Name reservations (`N14`) — ⛔ never reuse an imported constant key for a varying object

Every spatially-varying object gets a **fresh** name, distinct from the imported constant keys `W_0`, `e_W`,
`rho_br`, `mu_R`, `v_0` (live in `scripts/S11b_exports.py`): the profiles `W_bg(x)`, `w₁(x)`, `μ_R,bg(x)`,
`m₁(x)`; any varying `ρ_br⁰(x)`, `ρ_4D⁰(x)`; the contrast `η`; a fresh `e_{W,bg}` if introduced (T-g). ⛔
Reusing a constant key is an `F9` false-equal (`G3`), and for `rho_br` it silently freezes `∇Σ_E⁰=0` and
drops a channel. Reserve `v_bulk_normal_0` for the drain; ⛔ never the glyph `v₀`/key `v_0`.

## 6 · Method, and exit semantics — ⛔ no VERDICT, ⛔ no nonzero exit on a physics disagreement

Inherit S11b's derivation method verbatim (`:301-380`): balance laws; the virtual-displacement rule;
variational (⛔ not ordinary-partial) derivatives; the sign conventions; the **three script clauses**
(`.claude/skills/build/SKILL.md`) and the structural rule — the physical symbols are combined by hand
**only** in §1's setup and §2's background ansatz; every other expression is **reached by computation**;
every control re-enters at the setup, ⛔ never at a result. ⛔ **No terminal `VERDICT`/`PASS`/`FAIL`** (`:543`).
⭐**[folded]** Inherit S11b's exact exit rule (`:518`): a **physics disagreement EMITS and CONTINUES** (exit
0); a nonzero exit is for **operational failure only**. ⛔ A guard may not exit nonzero on a nonzero
representation/shape residual — that residual is a computed object to print, not a crash.

## 7 · Dimensions & homogeneity

Restore units and report the dimension of every emitted object; `η`, `w₁`, `m₁` are dimensionless, `∇W_bg`
is length/length-scale. ⛔ A shape-derivative not dimensionally homogeneous with its `(η,slope)→0` limit is a
defect — report the check able to fail. ⇒ `S11CA_DIMENSIONS`.

## 8 · Controls — ⛔ FORM controls (rule 14); a coefficient rescale tests nothing

- **C-1 (form).** ⭐**[folded: the weak branch is removed]** Drop the **off-diagonal** `∇W_bg` component of
  the tilt structure in T-a/T-c across all `D_brane = 3` in-plane directions, re-run, and report the literal
  diff. ⛔ Setting two in-plane normal components equal is **only a profile restriction** (a wrong
  contraction can agree there) — ⛔ **not** an accepted form control; a coefficient rescale of `η` stays in
  the family; `η→0` is the known-vacuous uniform limit (`N6`).
- **C-2 (independence).** The one-sided corruption of §4.
⇒ `S11CA_CONTROL_FORM`, `S11CA_CONTROL_INDEPENDENCE`.

## 9 · Tag grammar — ⛔ both engines emit PARALLEL tag sets

One tag per named object (the `S11CA_*` above). ⛔ Do not let the two engines choose two names for one
object. Declare any engine-local tag so parity is meaningful. The SymPy engine imports the S11b LEDGER and
carries it forward; the Wolfram engine is **blind** and re-derives §§1–5 from this spec.

## 10 · Supplied vs. tested; report back (⛔ under 25 lines)

**Supplied** (unfalsifiable here): all of §1 (the S11b setup); the background ansatz forms of §2a; the
oriented level-set + orientation law of §3. **Tested** (the build's output): every `S11CA_*` shape-
derivative and residual, the admissibility premise/state (2d), the representation-invariance residual (§4),
the form/independence controls (§8). ⛔ State in the report that the §1 objects were supplied, so a passing
build does not read as if it verified them. Report back the emitted tags, their orders, and the two control
residuals — ⛔ no conclusions.
