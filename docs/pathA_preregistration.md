# Path-A (Promoted-Throat-Field) Branch-Realization — Pre-Registration Record

**STATUS: 🟡 DRAFT — for the USER's conceptual GO (stage 1 of a two-stage gate; NOT yet freeze-ready).** This
pre-registration declares a **new, independent model**: the promoted-throat-field (Path-A) branch. It is
**NOT** a continuation of the frozen effective-closure Stage-1 pre-registration (`docs/stage1_preregistration.md`).
Per the non-rescue rule (`docs/branch_realization_parent_status_decision.md` point 5; effective-closure prereg
§A/§K), Path A is allowed only as a separate, independently pre-registered model with its own audit trail — a
Path-A result may not be advertised as a continuation of the effective-closure Stage-1 attempt, and vice versa.

**Two-stage gate (this is stage 1):**
1. **Conceptual GO (USER, now):** authorize the `S_Σ[R]` promotion + commit to the Path-A program, having seen
   the honest prior (§L). This DRAFT is the substantiation for that decision.
2. **Freeze (USER, later):** after the GO, Claude+Codex produce the **Path-A GATE-A freeze sheet** (the
   target-blind `S_Σ` constitutive forms + the closed-background extraction map + the hash definition — the
   REQUIRED pre-freeze deliverable, §E/§M); only then is this record completed and `frozen: YES` signed off,
   and only then may any solve/compute run (brief §3.4). **This draft is intentionally not freeze-ready: the
   new `S_Σ` constitutive packet is not yet pinned (§E).**

**Why this is USER-gated.** Promoting `S_Σ[R]` to an autonomous parent dynamical throat field is the one move
that is *conceptual*, not methodology (`parent_status_decision` point 4: "The only thing that WOULD make this
conceptual is … the autonomous-field/eigenmode lane"). Freezing this record + authorizing the self-consistent
solve is the user's call. The derivation that grounds it is target-blind analysis (done); the numerical solve
is future work, gated behind this freeze (brief §3.4 compute-spend gate).

**Purpose.** Fix — *in advance of any self-consistent solve* — the promoted parent action, the closed
matter/gauge↔wall loop, the gauge/boundary conventions, the frozen `S_Σ` constitutive packet, the derived
return sources, the extraction formulas, the observables, the targets, the primary observable, and the
non-rescue rules (including the sharpened no-fine-tuning rule from the crux resolution). No equation is
reconstructed here: every formula is quoted from the source ledger or the Path-A derivation with a citation.

**Companion documents.** Effective-closure prereg `docs/stage1_preregistration.md` (the model this supersedes
for Path A); scientific brief `docs/branch_realization_brief.md`; execution plan
`docs/branch_realization_execution_plan.md`; parent-status decision
`docs/branch_realization_parent_status_decision.md`; the crux resolution
`software/stage1_solver/decisions/08_pathA_crux_resolution.md`; the grounding derivation
`software/stage1_solver/derivations/pathA_01_return_source_and_balance.md`. Canonical equation source:
`notes/moving_throat_pde_program_compact.md` (cited as *compact*).

---

## §A. Pre-registered parent-status decision (the conceptual change)

**`parent_action_status = promoted_throat_field`** (Path A).

The parent action is promoted from the effective-closure form to the parent-complete moving-throat statement
(*compact* §2.4 L496-505):
```
S_total = S_ψ[ψ,A;Σ] + S_EM[A] + S_Σ[R].
```
`R(Ω,w,t)` becomes an **autonomous parent dynamical field** (no longer a frozen geometry / confinement-only
argument). At quadratic order around a stationary branch, `S_Σ → S_η^(2) + O(η³)` (*compact* L213-216), so the
quadratic wall operator is unchanged in form; what changes is (i) `R0(w)` is now a **derived self-consistent
equilibrium** of the promoted action (not a posited shape), and (ii) the matter/gauge→wall **return loop is
closed** by the derived sources `S_η^(ψ,A)`, replacing the effective-closure `f_ext` surrogate.

This is a deliberate departure from the effective-closure parent-status (`parent_status_decision` §"Decision";
effective-closure prereg §A). Path B (effective closure) remains a valid, separately-recorded falsification of
its own branch; this record does not retroactively reinterpret that result.

---

## §B. Why Path A is the right next test (the localized open question)

The effective-closure Stage-1 run produced a robust, target-blind, dual-engine `R_norm = −10.8` miss (the
forward transfer `P0 = N0/D0 = 6.6e-7` is ~7 orders below the GR target `54/5`; `STAGE1_VERDICT.md`). The crux
resolution (decision-08, Claude+Codex) localized the open question precisely:

- The miss is **not** rescuable by deriving a "missing" wall↔gauge coupling: reciprocity is exact at the action
  level (the return kernels are adjoints/variations of the *same* parent kernels — `S_η^(ψ) = −k1·δρ` with
  `k1 = V_wall'(Σ0/ℓc)/ℓc`; `S_η^(A)` matter-mediated), so closing the loop introduces **no new numerator
  magnitude** in `N0 = (P/Δ)²`, `P = Ω_U²G_W + R·G_U` (derivation D1/D2/D4).
- The effective-closure run **froze `R0(w)` by hand** (the GATE-A cubic smoothstep, `decisions/07`). Path A's
  only legitimate lever on `R_norm` is the **self-consistent background / denominator `D0 = K − B0 − Z0`**:
  promoting `S_Σ[R]` makes `R0(w)` a *derived equilibrium* whose bundle `(K, B0, Z0, N0)` differs from the
  hand-frozen one.
- This is exactly the self-consistent feedback closure the execution plan always intended:
  *"the feedback structure (inflow pressure ↔ wall energy ↔ standing-wave support ↔ brane closure) is what
  determines those coefficients on the realized branch … a failure tells you which sector of the feedback
  didn't close"* (`branch_realization_execution_plan.md` §5). The effective-closure branch ran with the
  feedback **open**; Path A closes it.

So Path A is the decisive `R_norm` re-test: does the *self-consistent* equilibrium land where the
effective-closure shortcut could not?

---

## §C. Branch identity + fields (what is solved)

- **WP1 — stationary isotropic promoted moving-throat branch.** Stationary (drops `t`), isotropic/spherically-
  symmetric (drops `Ω`); fields on the `(r,w)` 2-D mesh (full `r²` measure), same arena as the effective-
  closure branch (`branch_realization_execution_plan.md` §2.1) — but `R0(w)` is now **solved**, not frozen.
- **WP3 — grouped real `P2` weak-axisymmetric tangent**, linearized around the converged closed WP1 branch.
- **Closed coupled fields** (the promotion makes `R` dynamical; the return loop is on):
  - matter `ψ` (or `ψ_real, ψ_imag`), `ρ = |ψ|²`;
  - localized gauge `A_M = (A_0, A_i)`, `F_{MN}`;
  - **autonomous throat field `R(Ω,w,t) = R0(w) + η`** with `R0(w)` the derived static equilibrium;
  - the same gauge-fixing scalar / Lagrange multiplier (H=Z) and continuation parameters.
- Geometry: open finite throat, `Σ(X,t) = r − R(Ω,w,t)`, `R(0)=a`, `R(L)>0` (*compact* §12.2, L310-325).

---

## §D. Geometry, boundary, and gauge conventions (frozen)

Same conventions as the effective-closure prereg §D (open finite-throat; Dirichlet mouth anchoring `R(0)=a`;
open-impedance Robin exit `R_exit>0`, no hard cap; gauge **H=Z**, `ξ_4=ξ`; projection-first Maxwell
`W,Z,H=Z,S=Z/Z_int`; natural point-particle source-map `N_Q = chi_Q⁻¹`, `m̂0² chi_Q N_Q = 1`), with these
Path-A specializations:

| Item | Path-A choice | Source |
|---|---|---|
| Throat shape `R0(w)` | **DERIVED self-consistent equilibrium** of the static balance (§F), NOT a frozen profile | derivation D3; decision-08 |
| Wall closure | **Promoted `S_Σ[R]`** (autonomous), reducing to `S_η^(2)` at quadratic order | *compact* §2.4 L496-528 |
| Return loop | **CLOSED** — driven by the derived `S_η^(ψ,A)`; the `f_ext` surrogate is removed | derivation D1/D2 |
| EM localization weight | `Z(w)` (geometry-independent in the declared action — **no direct η→Z source**) | *compact* L590-620, 678-685; derivation D2 |
| Moving-EM-boundary terms | **Out of scope** for this record (the declared action uses fixed-domain `Z(w)`; a shape-calculus moving-boundary EM action is not declared and is not added here) | derivation D2 open-item; *compact* L310-325 |

**Source-physics scope (frozen).** Active: the closed promoted WP1/WP3 PDE source/current terms, projection-
first Maxwell, the natural source-map export, and the WP3 `P0/P2` response as extracted from the frozen tangent
before residual evaluation. **Excluded** (as in the effective-closure prereg §D): the post-242 relaxed/open-
system §12.12 companion lanes, and any moving-EM-boundary shape-calculus action. They may not be added after
residuals or used to reinterpret a miss.

---

## §E. Promoted wall action + constitutive packet (frozen — the new GATE-A analog)

The promoted wall action is (*compact* §2.4 L510-528):
```
S_Σ[R] = ∫ dt dw dΩ L_Σ,
L_Σ = ½ μ_Σ R_t² − ½ T_{w,Σ} R_w² − ½ T_{Ω,Σ}|∇_Ω R|² − U_Σ.
```
Its quadratic reduction around `R0(w)` reproduces the effective wall coefficients (derivation D3, verified
against *compact* L523-528, L1269-1305):
```
μ_η = μ_Σ(R0,w),  T_w = T_{w,Σ}(R0,w),  T_Ω = T_{Ω,Σ}(R0,w),
K_η = U_{Σ,RR}(R0,w) − ∂_w(T_{w,Σ,R}(R0,w) R0') + ½ T_{w,Σ,RR}(R0,w)(R0')².
```

**Posits (red-team `free_choice`):** the constitutive functions `μ_Σ, T_{w,Σ}, T_{Ω,Σ}, U_Σ` (and hence their
reductions `μ_η, T_w, T_Ω, K_η`). **The promotion adds throat dynamics; it does NOT derive these from a deeper
material law** (derivation D3; `parent_status_decision` point 3 — a promoted nonlinear `S_Σ[R]` is a *fresh*
parent-level posit). Under the calibrate-and-predict methodology (§H), these split into two parts that are
treated differently:

- **Functional FORMS — fixed target-blind, by physics, with a BOUNDED parameter count.** The shapes of `μ_Σ,
  T_{w,Σ}, T_{Ω,Σ}, U_Σ` are chosen on physical/structural grounds (smoothness, positivity, natural wall
  physics), NOT by looking at any target. This is what bounds the degrees of freedom and keeps the DOF count
  honestly small (no hidden parameters — §K).
- **PARAMETERS within those forms — may be calibrated to the stated anchor** (`R_norm` = the GR quadrupole),
  then frozen (§H). This is the legitimate, openly-declared calibration; the held-out surplus, not the anchor,
  is the test.

⚠ **REQUIRED PRE-FREEZE DELIVERABLE (this DRAFT does NOT yet pin them).** Unlike the effective-closure
constitutive functions (which had existing red-team provenance), the Path-A `S_Σ` forms are new and are **not
specified in this draft**. Before this pre-registration can be frozen (`frozen: YES`), a **Path-A GATE-A freeze
sheet** must be produced (Claude+Codex math determination, in the spirit of `decisions/07`) giving: the explicit
target-blind functional forms + domains + units + rationale + **the bounded calibrated-parameter list (the DOF
count)**; the closed-background field→coefficient extraction map; and the fixed `candidate_freeze_hash_def`.
Produced **after** the user's conceptual GO and **before** the freeze (§M, two-stage gate). The user's gate is
the conceptual GO (the `S_Σ` promotion); the specific forms + the calibration plan are then a Claude+Codex
determination, recorded and frozen before any held-out claim.

---

## §F. Derived return sources + extraction (frozen — from the Path-A derivation, applied to the closed system)

**Return sources** (derivation D1/D2; reciprocity ⇒ no new free magnitude):
```
S_η^(ψ) = −k1 δρ,            k1 = V_wall'(Σ0/ℓc)/ℓc      (matter return; reciprocal of forward δV_conf=−k1 η)
          (+ k2 ρ0 η self-stiffness folded into K_η,eff = K_η − k2 ρ0, k2 = V_wall''(Σ0/ℓc)/ℓc²)
S_η^(A) := δΓ_gauge,med[R]/δη   (gauge return; matter-mediated adjoint of η→δV_conf→δJ_ψ→δA — named,
          determined by the same covariant matter/current action; NO direct η→Z term)
```
The matter return is fixed to quadratic order in the declared confinement interaction; cubic+ source terms are
not claimed. The gauge return is the matter-mediated on-shell variation — a **closed local kernel requires the
full closed matter/gauge solve** (it is realized numerically by the closed solve, not as a posited formula).

**Static self-consistent throat balance** (derivation D3 — this is what determines `R0(w)`):
```
−∂_w(T_{w,Σ}(R0,w) ∂_w R0) + ½ T_{w,Σ,R}(R0,w)(∂_w R0)² + U_{Σ,R}(R0,w) = S_η^(ψ,A)|static[ρ0(R0), A0(R0)].
```
`R0(w)` solves this jointly with the static matter/gauge equations on the same geometry — the self-consistent
equilibrium.

**Extraction formulas (frozen, verbatim from the ledger — applied AFTER the closed solve).** Identical to the
effective-closure prereg §F (conservative bundle `D0=K−B0−Z0`, `D2,D4,u2,u4,P0,P2,P4`; outgoing scalar
`chi_Q, N_Q`, `Δ_norm`, `P0^target = 54Gc_s⁵/(5a⁵c⁵)`; prefactor slope `Ξ_1`; rigid-mouth `(U,V)`; placement
`varrho_phys`), now evaluated on the **derived self-consistent background**, with the mixed-Schur identities
(*compact* L3946-4010, verbatim; `R_mix` = the off-diagonal U↔W coupling, NOT the throat radius):
```
K_* := K − C²/ϖ²
Δ   := Ω_U²Ω_W² − R_mix²
Q   := G_U²Ω_W² + 2 G_U G_W R_mix + G_W²Ω_U²
P   := Ω_U²G_W + R_mix G_U
D0  := K_* − Q/Δ        (= K − B0 − Z0 packet form: B0 ⊃ C²/ϖ², Z0 = Q/Δ)
Λ   := P/Δ,   N0 = Λ²,   P0 = N0/D0
```

---

## §G. Observables to extract (target-blind)

Same as the effective-closure prereg §G: the bundle `(D0,D2,D4,N0,N2,N4,P0,P2,P4)`; source-map scalars
`(m̂0, S_port, chi_Q, N_Q)`; prefactor slope `Ξ_1`; orbit-lock chart `(U,V)`; placement `varrho_phys`. Each
carries a stated uncertainty (§J). Grouped by their role under the calibrate-and-predict frame (§H):

- **Calibration-feasibility readout (PRIMARY RESULT):** the **derived `R0(w)` profile**, the **realized `D0` and
  its margin to the Schur boundary** (`K − (B0+Z0)`), and whether a **physical, stable `S_Σ`** reached the GR
  near-pole `D0 → 0` (the calibration), with the achieved `R_norm` and its uncertainty. This is the decisive
  Stage-1 output (§L). `R_norm`/`chi_Q` here are the *anchor*, not held-out predictions.
- **Held-out surplus (the predictive test):** the observables NOT used in calibration — at Stage 1 primarily
  `R_pole`, `P2`, `P4` — each with uncertainty, reported alongside the calibrated-parameter (DOF) count so the
  surplus is honestly weighable (§K).

---

## §H. Targets (used ONLY for comparison, AFTER extraction) + the mechanism statement

**Do not insert into the solve, tune any parameter to, choose a branch by, or set a stopping criterion from
these.** Same V2 target card as the effective-closure prereg §H:
```
R_pole = D0(B4+Z4) − 3(M+B2+Z2)² = 0
R_norm = m̂0² S_port N0/D0 − 54 G c_s⁵/(5 a⁵ c⁵) = 0
P2 = P4 = 0
R_tail = Θ_tail (c/c_s)³ − 1 = 0   (tail_sector_active = false; quoted for fidelity, NOT evaluated in this run)
```
WP3 tangent targets: `d ln R_tr = d ln R_target = d ln ε_η = 0`, `N_Q = 1`.

**Test methodology: CALIBRATE-AND-PREDICT (Standard-Model style — `[[feedback-calibrate-predict-methodology]]`),
NOT blind-freeze-and-hope.** This is the project's own bootstrap (`n=5 → β=3 → κ_PV → …`, from
`research/1pn_orbital_dynamics/` onward) and how the Standard Model is validated: calibrate a small,
honestly-counted parameter set against a TRUSTED known anchor, freeze, then judge by the **held-out predictive
surplus** = (independent observables matched) − (free parameters tuned). The failure mode to avoid is the
string-theory one (`#knobs ≥ #predictions`), NOT calibration itself.

- **Calibration anchor (stated openly):** `R_norm = 0`, i.e. the GR quadrupole `P0 = N0/D0 = 54/5` — the only
  externally-benchmarked Stage-1 constant (Peters 1964 / `benchmarks.yaml`). The `S_Σ` parameters (within the
  physically-fixed forms of §E) **may be calibrated** so the self-consistent equilibrium reproduces it.
  **`R_norm` is therefore the calibration anchor, NOT a held-out prediction** — we say so plainly.
- **Mechanism (decision-08 — binding):** Path A reaches the anchor **only through the denominator `D0` /
  self-consistent background** — never through a new numerator / free-`P` magnitude (reciprocity ⇒ no added term
  in `P`; `N0` may differ only via the re-solved bundle). So calibrating to `R_norm` = driving the
  self-consistent `D0` toward the Schur softening edge `D0 → 0` (`K → B0+Z0`).

**PRIMARY RESULT — the reframed central question:** *can a **physical, stable** `S_Σ` reach the GR near-pole
`D0 → 0` at all?* — i.e. does a physically admissible (stable, `D0>0`, passivity-respecting; §I) promoted wall
action exist whose self-consistent equilibrium calibrates the GR quadrupole? **For a toy analog, reaching it at
all is a genuine positive** (the self-consistent feedback CAN carry the GR-scale outgoing transfer); failing to
reach it with any physical `S_Σ` is a clean, strong negative (the model structurally cannot). Both outcomes are
informative — see §L.

**HELD-OUT surplus — the predictive test (judged by surplus-over-DOF).** The observables NOT used in
calibration. **Independence map (methodology review — calibrating `R_norm` fixes `D0 = N0/T`, `T :=
54Gc_s⁵/(5a⁵c⁵)/(m̂0² S_port)`, so every bundle observable that contains `D0` is partially contaminated by the
anchor):**
- `chi_Q` / `Δ_Q`: **NOT independent** — same ω⁰ transfer as the anchor (`Δ_norm = P0_target(chi_Q⁻¹−1)`,
  `compact:6436-6454`). Not a held-out prediction.
- `R_pole = D0(B4+Z4) − 3(M+B2+Z2)²`: **partially held out** — the `D0` part is anchor-fixed; only the
  one-pole *shape* (`B4,Z4,M,B2,Z2`) is genuinely held out. Internal-consistency, not a second external target.
- `P2 = (D0 N2 − 2 D2 N0)/D0²`, `P4` (∝ `D0⁻²`, `D0⁻³`): **partially held out, strongly sharing the calibrated
  `D0`** — they test higher-moment bundle *correlations*, not an external benchmark.

**Honest Stage-1 scope:** report Stage 1 as **"anchor calibration within an N-parameter principled `S_Σ` family,
plus internal-consistency held-outs"** — the surplus is THIN (the GR quadrupole is the only external Stage-1
benchmark and it is the anchor; the held-outs share `D0`). The substantial *external* surplus is DOWNSTREAM
(the WP3 response sector `d ln R_tr/R_target/ε_η`, `N_Q`; Stage-2 divergence extraction; the EM-side
predictions), collected as the prediction set widens — exactly the bootstrap. The calibrated-parameter (DOF)
count and the realized surplus (with full calibration-covariance propagation, §J) are reported (§K).

---

## §I. Stability gates (frozen)

Same as the effective-closure prereg §I (wall positivity; stable BdG/Krein certificate; mixed-sector
positivity; `D0 > 0`; non-dark `N0`), **with the Schur boundary made explicit** (*compact*/v2_09): the branch
must remain on the **positive side** of `Δ0>0`, `K − Q0/Δ0 > 0` (i.e. `D0>0`). A branch failing any gate is
reported as a non-convergence/instability diagnostic, not retuned to pass. **The GR target corresponds to the
asymptotic approach to (not crossing of) this boundary** — a passive, stable near-pole, not an instability.
Every exported packet attaches `parent_action_status=promoted_throat_field`, `boundary_protocol`,
`stability_certificate`, `source_hashes`, `freeze_hash`.

---

## §J. Validation requirements (non-negotiable — brief §5; mandatory before any physics claim)

The five effective-closure prereg §J requirements carry over unchanged (validate on a known-answer limit;
≥3-level convergence study with per-observable order; boundary control below target signal; conservation
diagnostics; explicit error budget + noise floor per observable). **Plus two Path-A additions:**

6. **Self-consistent-balance validation.** The derived `R0(w)` must come from a validated solve of the §F
   static balance (existence/convergence demonstrated; the equilibrium is a genuine solution, not a partial
   iterate). The promoted closed Newton solve must converge with a documented residual.
7. **Margin uncertainty.** The realized `D0` and its margin to the Schur boundary (`K − (B0+Z0)`) carry an
   explicit uncertainty; whether the equilibrium "sits near `D0→0`" is meaningful only relative to that
   margin and the §J error budget. (Note the stiff `P=Kρ⁵` conditioning failure mode, plan §6 — a poorly
   conditioned near-pole must be distinguished from a genuine physical one.)
8. **Calibration-covariance propagation into the held-out residuals (methodology review).** Because the held-out
   observables share the calibrated `D0` and carry inverse powers of it (`P2 ∝ D0⁻²`, `P4 ∝ D0⁻³`,
   `compact:6809-6828`), a small calibrated `D0` AMPLIFIES the calibration uncertainty into the surplus. Every
   held-out residual must be reported with the full calibration covariance propagated through (not just its own
   discretization error); a "surplus match" inside the amplified covariance is NOT a match. Margin-to-Schur-
   boundary error bars are mandatory.

A result without §J does not count as an answer (brief §6).

---

## §K. Freeze boundary + non-rescue rules (binding after freeze)

Standard freeze flags (as effective-closure prereg §K): `pre_target_freeze=true`, `target_blind=true`,
`no_post_residual_refit=true`, `candidate_freeze_hash` over {promoted action + S_Σ constitutive packet,
boundary protocol, gauge convention, source map, derived return sources, extraction formulas, solver
tolerances, mesh}; no target residuals / pass flags / target values in the frozen packet.

Standard non-rescue rules carry over (no post-residual mutation of source support / boundary class / gauge /
port normalization / extraction; no algebraically-projected zero-residual packet reported as a physical hit).
**Plus the calibrate-and-predict discipline (decision-08, refined by `[[feedback-calibrate-predict-methodology]]`
+ the methodology review) — this REPLACES the earlier blunt "no tuning to `D0→0`" rule, which wrongly forbade
legitimate calibration:**

**The TWO-FREEZE protocol (binding sequence; do these in order, each step logged):**
1. **Freeze the model class — target-blind.** Fix the `S_Σ` functional FORMS, domains, units, passivity/
   boundary class, the calibration objective, the optimizer + tie-breaker, and EVERY discrete family choice
   (§E). Hash this. No reference to any target. This bounds the DOF.
2. **Calibrate — only to the stated anchor.** Optimize the declared parameters (only) so the self-consistent
   equilibrium reproduces `R_norm = 0` (the GR quadrupole). Nothing else is tuned.
3. **Freeze the calibrated values + hash.** Record the calibrated parameters; no further change.
4. **Evaluate held-out observables — no residual refit.** Compute every held-out observable (§H) with full
   calibration-covariance propagation (§J). Whatever comes out, stands.

**The accounting rules:**
- **Calibrate-vs-predict separation (load-bearing).** No observable may be both calibrated-to AND claimed as a
  prediction. `R_norm`/`chi_Q` are the anchor → never reported as held-out hits. Every claimed prediction is
  held out of step 2.
- **Count ALL DOF honestly (methodology review).** The "knob count" includes: every calibrated coefficient,
  spline knot / basis amplitude, regularization weight, branch-selection switch, normalization choice, AND any
  post-hoc candidate-family selection (trying several form families and keeping the best counts as DOF). Report
  the count + the held-out surplus = (independent matches) − (DOF). Meaningful only with
  `#independent-predictions > #knobs` (the string-theory guardrail); a thin/negative surplus is reported as such.
- **No post-residual refit (unchanged).** After step 4, parameters are not re-tuned. A held-out miss stands.
- **No hidden freedom.** Forms are fixed in step 1 to bound DOF; hidden parameters that secretly let `S_Σ` also
  fit a held-out observable are forbidden (they convert surplus into fitting).
- **Non-continuation.** This promoted-throat-field model is distinct from the effective-closure Stage-1
  pre-registration. A Path-A miss/pass falsifies/supports **this** model; it may not be reported as a
  continuation of, or a rescue of, the effective-closure attempt (and vice versa). Each stands on its own
  frozen record.
- **No moving-boundary EM bailout.** If Path A misses, the excluded moving-EM-boundary shape-calculus terms
  (§D) may not be added post-hoc to reinterpret the miss; that would require its own new pre-registration.

---

## §L. Honest prior — what a pass and a miss actually mean (decision-08; brief §7 analog)

Recorded **before** any solve, so the verdict is not re-spun afterward. The question is **calibration
feasibility**, not a blind hit:

- **Feasibility is only informative relative to a NARROW, predeclared family (the load-bearing caveat —
  methodology review).** With a *broad/arbitrary* `S_Σ`, reaching `D0 → 0` is **trivial and meaningless**: one
  can simply prescribe `U_{Σ,R}` along `R0(w)` to satisfy the static balance and `U_{Σ,RR}` / derivative terms
  to set `K_η`, hitting any `D0 = ε > 0` (`compact:3946-3978`; derivation D3). So the question "can a physical,
  stable `S_Σ` reach the near-pole?" has TEETH only when "physical" = a *small, physically-principled, frozen
  material family* with a bounded, honestly-counted DOF (§E/§K). Feasibility within that narrow family is the
  meaningful question; feasibility under broad forms is not a result.
- **Feasibility alone is NOT surplus.** Even reaching `D0→0` within the narrow family is at most a *mild*
  positive (the constrained self-consistent feedback *can* carry the GR-scale transfer — operationally
  demonstrable analog, `branch_realization_execution_plan.md` §9 / `project-analog-framework-goal`, not an
  ontology claim). The scientific content is the **held-out surplus** (§H), judged by surplus-over-DOF.
- **Outcome — feasibility fails within the narrow family:** no small principled `S_Σ` reaches the near-pole
  without instability. A clean, strong, publishable negative — the *fully closed* model cannot carry the GR
  quadrupole on this branch under a principled material law, with the `D0`-margin + obstruction pinpointing why.
  Much stronger than the effective-closure miss.
- **Outcome — feasibility succeeds within the narrow family:** report it as **"anchor calibration achieved
  within an N-parameter principled family,"** then report the **held-out surplus**. Honest Stage-1 scope: that
  surplus is THIN and mostly internal-consistency (the GR quadrupole is the only external Stage-1 benchmark and
  it is the anchor; `R_pole/P2/P4` share the calibrated `D0` — §H). The real *external* surplus is DOWNSTREAM
  (WP3 response, Stage 2, EM side), collected as the prediction set widens — the bootstrap.
- Either way the result is **conditional on the calibrated `S_Σ` (forms target-blind + narrow, parameters
  openly calibrated to the GR anchor)** and is judged by the **held-out surplus**, never by the anchor it was
  calibrated to.

---

## §M. Two-stage gate + Freeze Block

**Stage 1 — Conceptual GO (USER) — ✅ GIVEN 2026-06-17.** The user authorized the `S_Σ[R]` promotion + the
Path-A program, having read §L (the honest prior) and the calibrate-and-predict reframe (§H, decision-09).
**This UNBLOCKS the target-blind solver build + validation** (the execution-plan §7 sequencing: extend the
self-consistent closed solver with a dynamical `R0(w)` + the derived return loop; MMS per operator; convergence;
conservation; error budget on TARGET-BLIND surrogate observables) — done with **parameterized/placeholder `S_Σ`,
NO calibration, NO `R_norm`/target comparison, NO physical export.** It does NOT unblock the frozen calibration
run (Stage 2).

**Stage 2 — Freeze (USER, after the GO + the required pre-freeze deliverables).** The freeze block is signed off
ONLY after the §E pre-freeze deliverables exist (the narrow `S_Σ` GATE-A freeze sheet + DOF count, the
closed-background extraction map, the validation plan). **The CALIBRATION RUN** (two-freeze protocol §K: calibrate
`S_Σ` to `R_norm`, freeze values, evaluate held-out) — the only compute that touches the anchor/targets — may run
ONLY after `frozen: YES` (brief §3.4). Target-blind tooling/validation (Stage 1) is not gated by this.

```
frozen: YES           # ✅ user-authorized 2026-06-18 (stage-2 freeze gate passed)
conceptual_go: YES    # ✅ user-authorized 2026-06-17 (stage-1 conceptual gate passed)
conceptual_go_by: user
conceptual_go_date: 2026-06-17
freeze_date: 2026-06-18
freeze_commit: TBD   # this commit sets frozen: YES; its SHA is recorded by the immediately following commit (a self-hash here is impossible)
source_revision: 8bd82b9             # repo HEAD at freeze stamp
candidate_freeze_hash: ed358569393fed5fc29c0c13286a07cd438db467da6c1bc663a09bb04b1691c9
candidate_freeze_hash_def: SHA-256 over the decision-11 §8 canonical spec (family forms+ties+g=0, geometry, source/port mhat0=S_port=1, calibration objective incl. the GR anchor constant 54Gc_s^5/(5a^5c^5), §5 extraction-formula strings, mode/branch selection, tolerances+mesh ladder, source_file_sha256); EXCLUDES any target residual / computed D0/R_norm / pass-flag. Stamped by src/stage1_solver/patha_gate_a_freeze.py → frozen/pathA_gate_a/<hash>/.
parent_action_status: promoted_throat_field
signed_off_by: user
```

**REQUIRED pre-freeze deliverables (must exist before stage 2; math determinations — Claude+Codex; user gates
GO + conceptual + freeze):**
1. the **Path-A GATE-A freeze sheet** (`decisions/`-recorded, decision-07 analog) — the **narrow,
   finite-dimensional, target-blind** `S_Σ` constitutive forms `μ_Σ, T_{w,Σ}, T_{Ω,Σ}, U_Σ` + domains + units +
   passivity/boundary class + physical rationale; the **calibration objective + optimizer + tie-breaker**; and a
   **complete DOF count** — every calibrated coefficient, spline knot / basis amplitude, regularization weight,
   branch-selection switch, normalization choice, and any post-hoc candidate-family selection (§E/§K). The forms
   must be narrow enough that feasibility is non-trivial (§L) and the surplus is honestly weighable;
2. the **closed-background field→coefficient extraction map** (how `K, B_n, Z_n, N0, D0` are read off the
   self-consistent solution);
3. the **engine decision** for the self-consistent closed solve (extend the torch WP1 to a dynamical `R` +
   closed return loop, vs a Mathematica build);
4. the **§J validation plan** for the self-consistent balance (existence/convergence of `R0(w)`).

None may be set to rescue a target after residuals (§K, incl. the no-fine-tuning rule). Items 2–4 are
build-design that may follow the conceptual GO; item 1 (+ a fixed `candidate_freeze_hash_def`) is the gating
artifact that converts this DRAFT into a freeze-ready record.
