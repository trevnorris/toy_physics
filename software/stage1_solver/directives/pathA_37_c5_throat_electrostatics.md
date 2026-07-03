# pathA_37 — C5-with-a-throat: does the electric field emerge as a throat-sourced COUNTERFLOW (Gauss + Coulomb + charge⊥mass + signed force)?

**Status:** v5 — folded the v4 flow design-review (`NOT_SOUND`, 7 blockers); re-review pending. **Type:** focused, flat-brane,
linearized consistency test on a **finite-throat punctured domain**, with a **cheap pre-check phase (P-1) run BEFORE the full
dual-engine derivation**. **Frozen baseline:** pathA_35 G0 light sector `T0_SHEAR_FROZEN(d9520d3819c3)` + the **established
gravity-flow sector** (pathA_29 drain; PN ladder — LEVERAGED, not re-derived). **Author:** orchestrator. **Date:** 2026-07-03.

---

## §0. Provenance and scope — from "deformation" to "flow" to "counterflow"

pathA_36: the flat brane carries **2 transverse photons at `c_γ²=μ_R/ρ_br`** (light), but the longitudinal `∇·u` mode could not be
gauge-removed. Reframe: it is the **`c_s` density/flow/gravity channel**, and its **static response to a charge is the electric field**.
- **v1–v3 (deformation):** electric = elastic brane deformation. GLM found a structural energy no-go — generic first-gradient
  elasticity forces `u_L~1/r²`, so the gradient energy `(∂u)²~1/r⁶ ≠` Coulomb `~1/r⁴` (self-energy `~1/a³ ≠ ~1/a`).
- **v4 (flow):** gravity works **because it is a flow** (kinetic energy `½ρv²~1/r⁴`), so electric = a throat-sourced flow. The v4
  design-review confirmed this is the right class of object, but a **one-fluid** flow `v_E` has `ρv_E =` the mass current → charge is
  mass by construction (`FAIL_CHARGE_EQUALS_MASS`).
- **v5 (counterflow — this directive):** the electric field is a **relative current / counterflow** between two nondissipative
  carriers of the one medium:
  `J_m = Σ_A ρ_A v_A = 0` (zero NET mass transport), `J_E = Σ_A σ_A ρ_A v_A ≠ 0` (electric flux), `U = ½Σ_A ρ_A|v_A|² > 0`
  (field-squared energy). This carries kinetic energy while moving **no net mass** → **charge⊥mass falls out**, and the `±w`
  orientation plausibly supplies the two counter-flowing carriers (charge-sign as a mechanism, not an assertion).

**What the counterflow buys (and why "form-only" is NOT a win).** A conserved current in 3D gives `∮J·dA=Q` (Gauss) and `1/r²`
(Coulomb) **almost automatically**, and `½ρv²` gives `~1/r⁴` energy. **These are tautological for any flow** — a headline built on
them alone is `FORM_ONLY_TAUTOLOGY`, NOT a pass. The gate's real content is the three cruxes below.

**THE THREE CRUXES (each first-class, able-to-fail):**
1. **charge ⊥ mass via a NATIVE counterflow.** The electric carrier must be a conserved current with a nonzero component in
   `ker(J_m)` (zero net mass transport) and positive energy — AND the two carriers must be **identified with native model structure**
   (the `±w` orientation, brane-vs-bulk flow, or existing DOF), not a bolted-on second fluid (a viscous/thermal normal component is
   inadmissible at `T=0`). One-fluid → `FAIL_CHARGE_EQUALS_MASS`; non-native carriers → `FAIL_CARRIER_NOT_NATIVE`/drift.
2. **static closure.** A stable static charge loses nothing → the steady counterflow must **close** (an internal relative-current loop
   or return path; `J_m=0`; no secular accumulation in the throat; no static radiation). Gravity has its explicit drain/return
   closure (pathA_29); the electric analog must be exhibited. Else `FAIL_STATIC_FLOW_NO_CLOSURE` / `FAIL_DC_LIMIT_NOT_FLOW`.
3. **the signed force (make-or-break).** The `±w` orientation must give **like-charges-repel / unlike-attract** via the **field
   interaction energy** — *distinct from* gravity's universal attraction (which is **advective**, a kinematic mechanism). The electric
   branch must be governed by the signed energy-gradient force with **zero net mass/advection** channel; if the throats move only by
   advection in a net-mass flow → `FAIL_SIGN_STRUCTURE` / `FAIL_ADVECTIVE_ELECTRIC_FORCE`.

**If all three pass with provenance:** gravity (mass current, advective, one sign) and electricity (native counterflow, signed
energy-force, `±w`) are the **same throat-sourced-flow mechanism**, differing only in the conserved quantity and its sign → **the
forces are connected.** **Any failure is a first-class, informative no-go** about why charge and mass differ.

**Scope.** Flat-brane, linearized, exterior of a finite throat `Ω=ℝ³\H_a`. **In scope:** the FORM (Coulomb/Gauss/energy — necessary but
tautological), the three cruxes, transverse-light non-disturbance, the gravity parallel, longitudinal radiation. **Out of scope
(downstream):** full charge *quantization spectrum* (Gate Q — but **sign inheritance from `±w` is proven here**); magnetism (Gate S);
`λγ` (Gate K); throat interior/mass (Gate T). Magnitude calibrated; the cruxes are the predicted surplus.

---

## §1. Methodology (BINDING)

1. **Analog, not derivation.** A genuine no-go (no native counterflow, no closure, gravity-like signs) is **first-class and welcome**
   ([[feedback-falsification-is-the-goal]]).
2. **Emergence is a DIAGNOSTIC, never an input.** Gauss/Coulomb must be **derived theorems** from the source-free exterior continuity,
   never imposed. A bulk `δ³` source or an imposed mouth flux is a result-emitter → reject.
3. **Freeze the current-family MENU, not "a flow" (v4 Blocker-1/2).** §2.2 freezes a **finite menu** of current families with
   predetermined pass/fail; the derivation may NOT invent or select a current after seeing results (`FAIL_ACTION_NOT_FROZEN`). The
   one-fluid branch is a **mandatory ablation** that must return `FAIL_CHARGE_EQUALS_MASS`.
4. **Pre-checks BEFORE the expensive run (§2.6).** The cheap algebraic pre-checks (counterflow nullspace, static closure,
   energy-vs-advection sign, mouth non-tautology) run FIRST; a pre-check no-go returns the verdict without the full derivation.
5. **`FORM_ONLY_TAUTOLOGY` is not a win.** Coulomb+Gauss+energy without the three cruxes is explicitly a non-winning tag; a positive
   headline REQUIRES charge⊥mass (crux 1) + closure (crux 2) + the signed force (crux 3).
6. **Native provenance.** The counterflow carriers, the current identity, and the throat coupling must trace to native model structure
   (`±w`/brane-bulk/existing DOF); only the **magnitude** is calibrated. Bolted-on structure is drift ([[feedback-calibrate-predict-methodology]]).
7. **Sign force is energy-not-advection (v4 Blocker-4).** The signed force must be the fixed-charge field-energy force `−dU_int/dR`
   with the electric carrier having zero net mass/advection; a purely advective throat motion is `FAIL_ADVECTIVE_ELECTRIC_FORCE`.
8. **Able-to-fail is mandatory** ([[feedback-decisive-test-not-tautological]]); every crux + FORM check must be reachable-as-FAIL.
9. **Dual-engine REQUIRED** (SymPy + Mathematica; `ENGINE_AGREE` on the headline — the nullspace result, closure, flux theorem, power
   law, energy + self-energy, the interaction sign, the verdict — each engine assembling the headline itself, no `x−x`).
10. **No result-emitters; negative verdicts get HARDER scrutiny** ([[feedback-negative-verdict-short-circuit]]).
11. **Dimensional firewall** ([[feedback-dimensional-consistency-check]]): units-restored, ≥2 able-to-fail dim ablations.

---

## §2. Frozen inputs, the current-family menu, and the pre-checks

### §2.1 Frozen sectors (LEVERAGE, do not re-derive)
- **Light (transverse):** MacCullagh `½μ_R(∇×u)²`; inertia `½ρ_br(∂_t u)²`; `c_γ²=μ_R/ρ_br`; 2 photons. Non-disturbance checked.
- **Gravity (the flow template):** a drain → static inflow `v_cm ~ 1/r²` (`∇·(ρv_cm)=Ṁδ³`), kinetic energy `½ρv_cm²`, Newtonian
  `1/r²` with **universal attraction via advection**; `research/4d_*pn*`, `pathA_29`. **Do NOT re-derive** — the electric branch is
  built as the DELTA (native counterflow, `J_m=0`, signed energy-force).

### §2.2 The FROZEN current-family menu (the primitive — pick nothing after the fact)

Freeze a carrier space `{v_A}` (A = 1…N components of the one medium) with `J_m = Σ_A ρ_A v_A`, an electric current
`J_E = Σ_A σ_A ρ_A v_A` (weights `σ_A` from native structure), and positive energy `U = ½Σ_A ρ_A|v_A|²`. The **menu of branches**, each
with predetermined pass/fail, is frozen HERE (not chosen by the executor):

- **Branch (i) — one-fluid (`N=1`, mandatory ablation).** `J_E ∝ J_m` → charge is mass. MUST return `FAIL_CHARGE_EQUALS_MASS`.
- **Branch (ii) — two-carrier counterflow (`N=2`, the candidate).** Requires a nonzero `(v_1,v_2) ∈ ker(J_m)` (i.e. `ρ_1v_1+ρ_2v_2=0`)
  with `J_E≠0` and `U>0`, AND the carriers + weights `σ_A` **identified with native structure** (`±w` orientation / brane-vs-bulk /
  existing DOF). A non-native second fluid → `FAIL_CARRIER_NOT_NATIVE`.
- **Branch (iii) — polarization / circulation current (optional).** Admitted only with a **concrete continuity law and kinetic
  energy** tracing to native structure; otherwise excluded.

**Naming (v4 Suggestion-1):** reserve `v_cm` for the mass/center-of-mass flow (gravity); use `v_rel` / `J_E` / `χ` for the electric
counterflow carrier. No relabeling `v_cm` as `J_E`.

**Relation to pathA_36 (consistency).** The static electric field is a **sourced particular solution** of pathA_36's second-class
`c_s` longitudinal sector (like Newtonian gravity), NOT a Maxwell first-class constraint. Forbidden: a stray free **`c_γ`** photon;
allowed/expected: the `c_s` sound mode (`LONGITUDINAL_SOUND_PRESENT`).

### §2.3 The frozen throat source (P0) — mechanical mouth datum; derived flux, not imposed
- **Finite-throat punctured domain** `Ω=ℝ³\H_a`; frozen action in `Ω`; a **mechanical/geometric mouth datum** on `∂H_a` (its `±w`
  orientation sign `s` + one amplitude `A_mouth`). **No bulk `δ³`.**
- **Derived flux is REQUIRED, imposed flux is FORBIDDEN (v4 Blocker-6 fix).** A mouth condition may fix a conjugate
  geometric/head/polarization variable, and the solved exterior BVP may then map it to a derived capacity/admittance
  `Q = C(a)·A_mouth` — **that is allowed** (Gauss remains a theorem). It may **NOT** directly fix `∮J_E·dA`, `n·J_E`, or `Q` as
  Neumann/flux data → `FAIL_MOUTH_FLUX_IS_GAUSS_INPUT`. Distinguish *derived capacity* from *imposed flux*.
- **Point limit only after the freeze;** core/rim affects only the calibrated magnitude (else `FAIL_CORE_DEPENDENT`). No-throat
  control is **`A_mouth → 0`** (with `Q` then derived = 0), not `Q→0` as a primitive.

### §2.4 Steady-state closure (crux 2, a frozen requirement)
A stable static charge consumes/loses nothing. The steady counterflow must close: an **internal relative-current loop / return path**
with `J_m = 0`, **no secular accumulation** of any carrier in the throat, and **no static radiation**. The closure must be exhibited
(the electric analog of the gravity drain/return closure); if it cannot be written, `FAIL_STATIC_FLOW_NO_CLOSURE`.

### §2.5 Ledger
Count new fields/constants (the carrier space `{v_A}`, weights `σ_A`, `A_mouth`, `s`, closure structure) as drift; **flag each as
native (traced to `±w`/brane-bulk/existing DOF) vs new structural postulate vs calibration input.** Report `DRIFT(n)`.

### §2.6 PRE-CHECK phase (P-1 — cheap, dual-engine, runs BEFORE the full derivation)
Run these algebraically first; a no-go here returns the verdict without the expensive Green's-function/interaction derivation:
1. **Counterflow nullspace check.** For each frozen branch, form the linear maps `{v_A} → (J_m, J_E)`; is there `(v_A) ∈ ker(J_m)`
   with `J_E≠0` and `U>0`? Branch (i) fails by construction; branch (ii) passes iff the electric weights are non-proportional AND
   native. Output: pass/fail + the native-carrier identification.
2. **Static closure check.** Integrate continuity through the finite throat; admissible iff a return path / internal loop exists with
   `J_m=0`, no secular depletion. Else `FAIL_STATIC_FLOW_NO_CLOSURE`.
3. **Energy-vs-advection sign check.** Two fixed throats at separation `R`, signs `s_i=±1`: compute (a) `F_energy=−dU_int/dR` at fixed
   derived `Q_i`; (b) the advective throat velocity in the other's flow. Admissible iff the electric branch has zero advection and
   `F_energy` gives like-repel/unlike-attract; else `FAIL_SIGN_STRUCTURE`/`FAIL_ADVECTIVE_ELECTRIC_FORCE`.
4. **Mouth non-tautology check.** Freeze one non-flux mechanical mouth datum; verify (symbolically, capacity-level) that `Q` is an
   OUTPUT and not secretly `n·J_E`/`∮J_E·dA`; else `FAIL_MOUTH_FLUX_IS_GAUSS_INPUT`.

---

## §3. The full computation, and the verdict grammar

If the P-1 pre-checks pass, the full dual-engine derivation computes: the exterior continuity/EOM + the **derived flux theorem**; the
Green's function → power law + arbitrary-surface flux; the field energy density + self-energy scaling; the two-throat interaction
energy + **signed force** (vs `±w`) with advection explicitly separated; the transverse non-disturbance; the provenance status; all §4
controls.

**Verdict grammar.**
- **`C5_THROAT_ELECTRIC_COUNTERFLOW_WITH_PROVENANCE`** — all three cruxes pass with native provenance: a native counterflow (`J_m=0`,
  `J_E≠0`, positive energy) gives Coulomb + Gauss (theorem) + `~1/r⁴` energy + `~1/a` self-energy; static closure holds; the `±w`
  sign structure yields like-repel/unlike-attract via the energy-force (zero advection); transverse light undisturbed; magnitude
  calibrated. **THE WIN** — electricity is a throat-sourced counterflow, the same mechanism as gravity, differing in carrier + sign →
  **the forces are connected** → Gate Q.
- **`C5_THROAT_ELECTRIC_COUNTERFLOW_BY_TUNING`** — the cruxes hold only with a tuned/non-native current or coupling. Report the condition.
- **`FORM_ONLY_TAUTOLOGY`** — Coulomb+Gauss+energy hold but a crux is unmet or untested. **NOT a win**; not green for Gate Q.
- **`FAIL_CHARGE_EQUALS_MASS`** — the admissible current is `∝` the mass drain (one-fluid, or no native `ker(J_m)` carrier). The
  primary physical no-go.
- **`FAIL_CARRIER_NOT_NATIVE`** — charge⊥mass only via a bolted-on second fluid not traceable to native structure (drift).
- **`FAIL_STATIC_FLOW_NO_CLOSURE` / `FAIL_DC_LIMIT_NOT_FLOW`** — no steady closed counterflow for a static charge (secular
  accumulation / a static radiative leak / the `ω=0` limit is not a genuine steady current).
- **`FAIL_SIGN_STRUCTURE` / `FAIL_ADVECTIVE_ELECTRIC_FORCE`** — signs are gravity-like, or the only motion is advective (not the
  signed energy-force).
- **`FAIL_FORCE_LAW_UNDEFINED`** — the two-throat force cannot be given an unambiguous energy-vs-advection prescription.
- **`FAIL_NO_COULOMB_FIELD` / `FAIL_FIELD_ENERGY_MISMATCH`** — not `1/r²` / surface-dependent flux / energy not field-squared `~1/r⁴`.
- **`FAIL_SOURCE_NOT_FROZEN` / `FAIL_MOUTH_FLUX_IS_GAUSS_INPUT` / `FAIL_ACTION_NOT_FROZEN` / `FAIL_COUNTERFLOW_NOT_IN_FAMILY`** — the
  source/current was invented, imposed as a mouth flux, cherry-picked, or a needed branch was excluded from the frozen menu.
- **`FAIL_TRANSVERSE_DISTURBED` / `FAIL_CORE_DEPENDENT` / `FAIL_CHARGE_NONCONSERVATION`** — as in v4.
- **`INCONCLUSIVE_SIGN_NOT_TESTED`** — the sign crux was not decided (explicitly **not green**, not a positive result; replaces the
  removed soft `ELECTRIC_FLOW_PARTIAL`).
- Sub-tags: `LONGITUDINAL_SOUND_PRESENT`; `GRAVITY_PARALLEL_HOLDS|NOTE_ONLY` (structural comparison, **not** a win condition);
  `LONGITUDINAL_RADIATION(...)`; `DRIFT(n)`; `AXIS_RE_ADMITTED`; `U_W_COLLISION`.

---

## §4. Able-to-fail controls (each a real computation; each must fire correctly)

1. **One-fluid ablation (mandatory).** Branch (i) MUST return `FAIL_CHARGE_EQUALS_MASS`. *Proves the gate can see charge=mass.*
2. **Counterflow nullspace (crux 1).** Branch (ii) passes iff a native `ker(J_m)` carrier with `J_E≠0`, `U>0` exists; a
   proportional-weights ablation (`σ_A ∝ 1`) must fail. *Proves charge⊥mass is derived, not asserted.*
3. **Native-carrier control.** A bolted-on non-native second fluid must fire `FAIL_CARRIER_NOT_NATIVE`. *Guards drift.*
4. **Static closure (crux 2).** A configuration with secular throat accumulation / static leak must fire `FAIL_STATIC_FLOW_NO_CLOSURE`;
   a closed relative-current loop passes. *Proves the static charge is steady and loses nothing.*
5. **Energy-vs-advection sign (crux 3).** `(+,+)`/`(−,−)` must **repel**, `(+,−)` **attract**, via `−dU_int/dR` at fixed `Q_i` with
   zero advection; the all-attract (advective/gravity-like) ablation must fire `FAIL_SIGN_STRUCTURE`/`FAIL_ADVECTIVE_ELECTRIC_FORCE`.
6. **FORM-only tautology guard.** A run that yields Coulomb+Gauss+energy but skips a crux MUST be tagged `FORM_ONLY_TAUTOLOGY`, not a
   win. *Prevents a false pass on the easy part.*
7. **Deformation contrast (demoted).** Recompute the elastic `1/r⁶`/`1/a³` scaling as a contrast (verifying the GLM no-go); it
   motivates the pivot but does **not** make the flow result non-tautological. *Documentation, not a discriminator.*
8. **Source-freeze: capacity vs imposed flux.** A derived `Q=C(a)A_mouth` from a non-flux mechanical datum PASSES; a directly imposed
   `n·J_E`/`∮J_E·dA` FIRES `FAIL_MOUTH_FLUX_IS_GAUSS_INPUT`. Mutants: bulk `δ³`; cherry-picked current → `BY_TUNING`/`FAIL_*`.
9. **Arbitrary-surface Gauss-flux.** `∮J_E·dA` equal on two radii AND a non-spherical surface (`=Q`); surface-dependent under a
   non-conserved-current ablation.
10. **Transverse non-disturbance.** `c_γ` unchanged by the counterflow sector + throat coupling; static charge excites no photon.
11. **Radiation / charge-conservation.** Static `Q` → no radiation; moving fixed-`Q` → report `LONGITUDINAL_RADIATION(...)`; `Q(t)`
    primitive → `FAIL_CHARGE_NONCONSERVATION`.
12. **Core-dependence.** Vary rim/core at fixed `Q`; far-field exponent unchanged (else `FAIL_CORE_DEPENDENT`).

---

## §5. Dual-engine + dimensional firewall (REQUIRED)
- **SymPy (.py) and Mathematica (.wl)** independently run the §2.6 pre-checks and (if passed) the full derivation; `ENGINE_AGREE` on
  the nullspace result, closure, flux theorem, power law, energy + self-energy, interaction sign, and verdict — each engine assembling
  the headline itself (no `x−x`).
- **Dimensional firewall** — every term (inertia, MacCullagh potential, carrier kinetic energies `½ρ_A v_A²`, `J_m`, `J_E`, `c_γ²`,
  `c_s²`) units-restored; **≥2 able-to-fail dim ablations that MUST fire**; carrier-free check per the Gate-4 lesson.
- **Timeouts** — every runner `timeout 600`; a 124 is a failure → reformulate, never raise the cap.

---

## §6. Discipline / scope
- Flat-brane exterior of a finite throat; magnetism, quantization spectrum, `λγ`, throat interior out of scope except as
  non-disturbance / sign-inheritance checks (sign inheritance IS proven here).
- **Do not rescue.** Report any crux FAIL plainly with the physical reason — that reason is the deliverable's most valuable content.
- **Do not bolt on an abstract gauge field, or a non-native second fluid.** Emergence must be from the native throat-sourced
  counterflow.
- Orchestrator reviews + arbiter-re-runs; Codex codes/runs/iterates; no orchestrator code mutation ([[feedback-claude-reviews-codex-codes]]).
- **Deliverables:** `tools/pathA_37_throat_{sympy.py,.wl}`; `reports/pathA_37_c5_throat_electrostatics.md` + `_results.yaml`; the P-1
  pre-check results + (if passed) the verdict + DERIVED nullspace / closure / flux theorem / power law / energy + self-energy /
  interaction sign + provenance + all §4 controls + dim firewall + the `LONGITUDINAL_RADIATION` report.

---

## §7. Acceptance (what "verified" means)
**Nothing is accepted until:** (a) **derived via BOTH SymPy and Mathematica** with `ENGINE_AGREE` on the headline; AND (b)
**double-reviewed by Claude AND Codex** — this v5 design-reviewed to green (Codex + Claude) before execution (a full re-review, not a
confirm-pass; a GLM tertiary is optional since GLM/Codex authored the flow/counterflow ideas — treat them as non-independent here);
then post-execution the standing tri-review (arbiter re-run + fidelity audit + adversarial-with-ablation on **fresh clean agents** +
Claude's read). **A negative verdict gets harder scrutiny** (confirm the FAIL is DERIVED and **able-to-pass**). **A positive verdict
gets the ablation scrutiny too** — confirm the three cruxes genuinely flip under the control-1/2/4/5 ablations, and that the headline
is NOT `FORM_ONLY_TAUTOLOGY` (a clean "it all works" is suspicious).

---

## §8. Changelog
- v1–v3 (2026-07-03) — electric-as-elastic-**deformation** gate; hardened to `SOUND_AS_IS` before GLM found the energy-scaling no-go.
- v4 (2026-07-03) — **PIVOT to electric-as-FLOW** (GLM + user). Cruxes: charge⊥mass + `±w` sign. Design-review `NOT_SOUND`: a one-fluid
  `½ρv²` is charge=mass by construction; needs an explicit counterflow family + closure + a sharpened sign test.
- v5 (2026-07-03) — folded the v4 design-review (7 blockers) + its 4 cheap pre-checks. **Electric = a native COUNTERFLOW** (relative
  current, `J_m=0`, `J_E≠0`, `U>0` → charge⊥mass). Changes: (1) froze an explicit **current-family menu** (§2.2: one-fluid = mandatory
  FAIL ablation; two-carrier counterflow with **native** carriers; optional polarization/circulation) — no post-hoc current selection
  (Blocker 1/2). (2) added the **P-1 pre-check phase** (§2.6: nullspace, closure, energy-vs-advection sign, mouth non-tautology) run
  BEFORE the full derivation (the 4 cheap decisive checks). (3) **static-closure** requirement + tags (§2.4; Blocker 3). (4) **signed
  force = energy-not-advection**, two parallel computations (Blocker 4). (5) removed `ELECTRIC_FLOW_PARTIAL` as a success →
  `INCONCLUSIVE_SIGN_NOT_TESTED` + added `FORM_ONLY_TAUTOLOGY`, `FAIL_CARRIER_NOT_NATIVE`, `FAIL_STATIC_FLOW_NO_CLOSURE`,
  `FAIL_DC_LIMIT_NOT_FLOW`, `FAIL_ADVECTIVE_ELECTRIC_FORCE`, `FAIL_FORCE_LAW_UNDEFINED`, `FAIL_COUNTERFLOW_NOT_IN_FAMILY` (Blocker 5,
  Suggestion 2). (6) **source-freeze fix** — a *derived* capacity `Q=C(a)A_mouth` is allowed; only a directly imposed flux fails;
  no-throat control `A_mouth→0` (Blocker 6). (7) `FORM_ONLY_TAUTOLOGY` guard so the easy part can't false-pass (Blocker 7); carrier
  renamed `v_rel`/`J_E` vs `v_cm`; native-vs-drift ledger. Next: Codex design-review of v5 → fold → execute (P-1 pre-checks first).
