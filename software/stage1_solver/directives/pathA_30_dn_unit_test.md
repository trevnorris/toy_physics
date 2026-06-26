# Directive pathA_30 — Phase-1 Frozen-Wall D/N Unit Test

**STATUS:** ✅ DONE — executed 2026-06-25, verdict **`DN_UNITTEST_BC_DEPENDENT`** (tri-review CLEAN: arbiter
re-run reproduced both engines + fidelity CLEAN + adversarial CLEAN; user banked the D/N calibration). The
geometry lift's frozen reduction gives the const-coeff Helmholtz resonator (operator/speed `c_s`/domain `L₀`/
DtN `−(ω/c_s)tan`/half-shifted ladder/Robin counterfactual all genuine, dual-engine agreed); D/N recorded as
an **imposed calibration input** (`bc_derivation_emitted=False`) — earning it (→ `PASS`) is an optional later
upgrade. Two non-blocking NITs logged in the completion ladder. Review trail: Codex (8) → re-confirm → GLM
(Bogoliubov `k⁴` + 5) → Codex post-GLM (2) → close-out `SOUND`, xhigh throughout.
**Push:** the full nonlinear moving-throat PDE closure (memory `project-moving-throat-pde-push`). This is the
**first, cheapest able-to-fail compute gate** of that push (roadmap Phase 1; scaffold §9 + §11 item 1).

**One-line goal.** Show that the **distributed geometry lift** `Σ=r−R(Ω,w,t)` with `V_conf(Σ/ℓ_c)`,
linearized about the stationary straight reference throat and **frozen** (`η=0`), **reduces** the
longitudinal-along-throat matter support channel to the exact finite-throat D/N benchmark — mouth DtN
`Z₀₀(ω)=−(ω/c_s)·tan(ωL₀/c_s)`, pole ladder `ω_n=(πc_s/L₀)(n+½)` — and that a *wrong* lift/BC **fails** this.

**VERDICT LADDER (the script must compute, not assert, which rung fires; `BC_DEPENDENT` is the DEFAULT —
`PASS` must be *earned* per §2.4):**
- `DN_UNITTEST_PASS` — the frozen reduction of the lift yields const-coeff Helmholtz `d²/ds²+(ω/c_s)²` on
  `s∈[0,L₀]` with speed = bulk `c_s` and domain = throat depth `L₀` (emitted **reduction certificate**, §1);
  the mouth DtN equals the `tan` form (dual-engine, **independent routes**, `FullSimplify` agree); poles at
  the half-shifted ladder; round-trip closure `R_rt=1, φ₀≡0`; static limit `Z₀₀(ω→0)=0` *emerges* from the
  small-ω expansion; the Robin-family counterfactual (§3) is **rejected**; **AND** the script emits a genuine
  variational/asymptotic derivation that the D-at-mouth / N-at-cap assignment **follows from the lift's
  mouth/cap physics**.
- `DN_UNITTEST_BC_DEPENDENT` *(default)* — everything above holds **except** the BC derivation: the D/N
  assignment is **imposed**, not shown to follow from the lift. Honest toy-model outcome: the BC is a labeled
  modeling input (calibration), not a derivation. Records provenance; does **not** block the push. **This is
  the verdict unless the BC-provenance derivation is actually emitted.**
- `DN_UNITTEST_FAIL_<reason>` — `FAIL_OPERATOR_INTRUSION` (an unsuppressed variable-coefficient / `Q` /
  `δV_conf` term survives in the reduced operator `L_s`, §1/§2.1); `FAIL_WRONG_SPEED` / `FAIL_WRONG_DOMAIN`;
  `FAIL_POLE_LADDER` (poles off the `(n+½)` ladder); `FAIL_COUNTERFACTUAL` (one Robin determinant does **not**
  reproduce D/N at `α→0` and D/D at `α→∞`, or the half-shift survives `α→∞` ⇒ hardcoded lookup);
  `FAIL_ENGINE_DISAGREE` (the two engines' `FullSimplify` do not agree).
- `DN_UNITTEST_ENV_BLOCKED` — missing Mathematica / license / timeout / infrastructure failure prevented the
  dual-engine check. **Not a physics verdict** — re-run; never report infrastructure failure as a `FAIL`.

---

## §0. Scope, framing, and what makes this a real test (not a tautology)

**Toy-model contract (standing).** Calibrations are allowed but every claim must be *earned*, not asserted;
we may later derive from first principles if we find a way. The goal of the whole push is to demonstrate a
**legitimate brane/bulk/soliton structure that supports EM+GR** — this gate is one necessary plumbing check
on the way there. (memories `project-calibrated-pde-goal`, `feedback-calibrate-predict-methodology`.)

**Nothing is static (keystone).** "Frozen wall" means the **geometry background** `R₀(w)` is the *stationary*
reference (`η=0`), **not** a static-field reduction. The matter perturbation `δψ ~ e^{-iωt}` is genuinely
**dynamic**; `Z₀₀(ω)` is the full ω-dependent mouth response. The static value `Z₀₀(ω→0)=0` must be obtained
as the **ω→0 limit of the one dynamic operator** (the small-ω expansion), never inserted as a separate
static solve. (memory `project-model-mechanics-corrections` §1.)

**Why this is able-to-fail and not "solve a textbook 1D resonator."** Re-deriving `Z₀₀` from *imposed*
constant-coefficient Helmholtz + *imposed* D/N BCs is tautological and is **explicitly forbidden as the
deliverable** (§4). The teeth are in the **reduction** — does the lift's frozen limit, built from the *actual*
GNLS matter sector + `P=Kρ⁵` EOS + `V_conf(Σ)` coupling, genuinely produce that idealized operator? A wrong
lift gives a variable-coefficient operator, a renormalized speed, a shifted domain, or the wrong BC — each a
real failure surface. Decompose and check each independently (§2).

**This gate is necessary, not sufficient.** Passing means the geometry lift's plumbing is internally
consistent with the finite-throat D/N branch before any grouped-`P2` (ℓ=2) work. It does **not** establish the
quadrupole normalization, the leak/return, or the Maxwell coupling — those are later phases (§7).

---

## §1. The object to build (the frozen reduction)

Start from the **already-frozen** parent stack — do not re-postulate it:
- matter: gauged 4D GNLS, `i ħ D_t ψ = [−(ħ²/2m)D_iD_i + V_conf(X;Σ) + h(ρ)]ψ`, `h=5Kρ⁴/4`,
  `c_s²(ρ)=5Kρ⁴/m` (handoff §2.1, §3.1; scaffold §1.1);
- geometry lift: `Σ(X,t)=r−R(Ω,w,t)`, `V_conf(X;Σ)=V_wall(Σ/ℓ_c)`, reference `Σ₀=r−R₀(w)`, perturbation
  `R=R₀+η`, `δV_conf=−V'_wall(Σ₀/ℓ_c)/ℓ_c · η` (handoff §5.3; scaffold §2, §4).

**Reference straight throat (the unit-test branch, scaffold §9.1):** a straight finite throat of depth `L₀`,
mouth at `s=0`, regular cap closure at `s=L₀` (`R₀(L₀)=0`). Use the longitudinal centerline coordinate
`s∈[0,L₀]`. Reduce the **frozen** (`η=0`) linearized matter support channel along `s` to its 1D ODE.

**Deliverable of this section (emit, don't assume):** the **actual reduced longitudinal ODE** obtained from
the linearized GNLS support channel on the frozen straight reference, written as
`L_s ψ̂(s) = 0` with `L_s` produced by the reduction — *then* compare it to the idealized
`d²/ds² + (ω/c_s)²`.

**Reduction certificate (mandatory, emitted).** The reduction is the able-to-fail core, so it must be
*auditable*, not asserted. The script must emit a certificate listing: (a) the stationary background
(`ψ₀(s), ρ₀(s)`; `A_{M0}=0` for this scalar unit test), (b) the projection that collapses the bulk operator
onto the longitudinal `s` channel — **name the measure (e.g. `√γ₀(s)`) and show it is constant, else it
induces a first-derivative term `f(s)d/ds` (§2.1)** — and (c) the **explicit fate of each potentially-intruding
term**, each with the *exact stated condition* under which it vanishes or is retained:
- background density variation `ρ₀(s)` and the induced `c_s(s)` variation — vanish at `O(ρ₀'/ρ₀)` on the
  uniform straight reference;
- background `V_conf(Σ₀)` — absorbed into the reference / vanishes off the wall;
- **the quantum potential, split into TWO distinct terms (this is the sharpest "nothing is static" trap):**
  (i) the **background** `Q`-gradient `∇Q(ρ₀)` — vanishes at `O(ρ₀'/ρ₀)` on a uniform background; **and**
  (ii) the **perturbation** `δQ = −(ħ²/4mρ₀)∇²δρ` — the **Bogoliubov `k⁴` term** (full BdG dispersion
  `ω² = c_s²k² + ħ²k⁴/4m²`), which does **NOT** vanish on a uniform background. The linearized matter sector
  is genuinely **BdG (2×2 coupled / 4th-order in hydrodynamic form)**; the reduction to scalar Helmholtz
  `d²/ds²+(ω/c_s)²` and the exact D/N ladder hold **only in the phonon / long-wavelength limit** `kξ ≪ 1`,
  where **`ξ = ħ/(mc_s)` is the healing length** (the scale at which `ħ²k⁴/4m²` becomes comparable to
  `c_s²k²`) — equivalently `ω ≪ mc_s²/ħ`. **`ξ` is distinct from the confinement length `ℓ_c`** in
  `V_wall(Σ/ℓ_c)`; do not conflate them. The certificate must **state `kξ≪1` explicitly** and record `δQ` /
  the `k⁴` term as an *explicitly-parameterized higher-phase term that vanishes at this order* under it.
  Dropping the `k⁴` term **without** the stated `kξ≪1` condition ⇒ `FAIL_OPERATOR_INTRUSION`.

Any potentially-intruding term that does **not** demonstrably vanish under a stated condition and remains in
`L_s` ⇒ `FAIL_OPERATOR_INTRUSION` — it may **not** be silently dropped or hand-waved into the §7 deferral
ledger.

---

## §2. The five able-to-fail sub-checks (each independently computed)

1. **OPERATOR.** Is the reduced `L_s` the constant-coefficient Helmholtz `d²/ds² + (ω/c_s)²`, or does the
   reduction introduce an `s`-dependent coefficient, a **first-derivative / measure-weight term `f(s)d/ds`**
   (from a non-constant projection measure e.g. `√γ₀`), a **`∇⁴`/`k⁴` Bogoliubov term** (from `δQ`, §1), or an
   extra potential term (from `ρ₀(s)` variation or `δV_conf`)? Emit `L_s` symbolically. **Decidable rule:**
   any such term not shown to vanish under a stated reduction-certificate condition (§1) is an `L_s` intrusion
   ⇒ `DN_UNITTEST_FAIL_OPERATOR_INTRUSION`. Only terms *explicitly parameterized as higher-phase* (Phase-2
   variable-coefficient corrections, or the phonon-limit-deferred `k⁴` term) may be recorded in the §7
   deferral ledger — and only after the certificate shows they vanish at this order.
2. **SPEED.** Trace the propagation speed in `L_s` to the bulk `c_s=√(5Kρ₀⁴/m)` from the EOS — confirm it is
   **not** a healing-length / wall-stiffness renormalized speed.
3. **DOMAIN.** Confirm the interval length is the throat depth `L₀` (cap from `R₀(L₀)=0`), not an effective
   length.
4. **BC TYPE + PROVENANCE.** Apply D-at-mouth / N-at-cap. **Classify, with the humble default:** the verdict
   is `BC_DEPENDENT` **unless** the script emits a genuine variational/asymptotic derivation that yields the
   D/N assignment as an **output, not an input**. To count as a derivation (anti-postulate-to-win, the
   pathA_29 lesson where a postulated BC selected the favorable outcome), it must emit: (a) the `V_conf`
   gradient at the **mouth** and at the **cap** from the wall profile `V_wall(Σ₀/ℓ_c)`; (b) the regularity
   condition at the pinch-off `R₀(L₀)=0`; (c) the BC that **follows** from each — Neumann at the cap from
   regularity, and at the mouth the clamp-to-bulk-reservoir condition. A bare "we impose D/N", or a narrative
   without (a)–(c), ⇒ `BC_DEPENDENT`, not `PASS`. **Framing caveat:** D-at-mouth here is the **fixed-geometry
   (η=0, Born–Oppenheimer) clamp to the quasi-static bulk reservoir**, NOT a radiation BC — the frozen-wall
   DtN is **lossless** (real-valued at all ω); the radiation impedance (imaginary part) only appears when the
   brane↔bulk coupling is restored in later phases (§7). Round-trip reflection coefficients: define
   `r_L ≡ r_N`; with `r_D=−1, r_N=+1` (handoff §7.3) the product gives the trapped closure — but state
   `R_rt=1, φ₀≡0` **only after** substituting the D/N pole ladder, not as a standalone identity.
5. **DtN + POLE LADDER.** From the BVP, derive the mouth DtN in the **outward-mouth convention**
   `Z₀₀ = −ψ'(0)/ψ(0)` (`s` increasing *into* the throat — this is the sign that yields the target) and
   simplify to `Z₀₀(ω)=−(ω/c_s)·tan(ωL₀/c_s)`; independently locate its poles as the zeros of the relevant
   transcendental denominator and confirm `ω_n=(πc_s/L₀)(n+½)`. **Dual-engine, genuinely independent
   algebraic routes** (e.g. SymPy: `dsolve` the BVP + impose BCs + form `−ψ'(0)/ψ(0)`; Mathematica:
   transfer-matrix / resolvent route — NOT a renamed replay of the SymPy algebra; cf. redteam stage_097
   lesson). Agreement via `FullSimplify[a−b==0]`.

**Plus two consistency controls:**
- **STATIC↔DYNAMIC (procedural, not an independent test — honest scope).** `Z₀₀(ω→0)=0` must come out of the
  small-ω series of the dynamic `Z₀₀(ω)` (`Z₀₀ = −(L₀/c_s²)ω² + O(ω⁴)`); emit the series. **Caveat:** this
  only verifies **Taylor-coefficient consistency of the derived DtN** — the `O(ω²)` onset is guaranteed by the
  `tan` form, so it is *not* an independent able-to-fail "one-object" test (unlike pathA_29's
  two-independent-computation `static_dynamic_consistency`). The one-object-two-limits principle is enforced
  **procedurally** here — by *forbidding* a separate static solve (§0) — not verified by an independent
  computation. Do not over-claim it as a decisive check.
- **ROUND-TRIP.** `R_rt = r_D r_L e^{2ikL₀} = 1`, `φ₀≡0 (mod 2π)` on the ladder (handoff §7.3).

---

## §3. Mandatory counterfactual self-test (the pathA_29 lesson)

The script **must** prove a wrong answer **fails its own gate** — this is non-negotiable after pathA_29's
three assert-not-solve rejections (hardcoded `r²` → spectral lookup → hand-written `1/r`).

**Parameterized Robin counterfactual (mandatory — defeats a label-keyed two-branch lookup).** A plain
D/D-vs-D/N swap can be passed by code that hardcodes *each* branch separately
(`if bc=='DN': tan; elif bc=='DD': cot`). To force a genuine parameterized solve, impose a **symbolic Robin
cap BC** `ψ'(L₀)=α·ψ(L₀)` and derive the mouth response from the **one** coefficient-matrix / determinant
**assembled from the `dsolve` general-solution coefficients** (NOT a looked-up Robin DtN formula that happens
to interpolate correctly) as a function of `α`. Then require that the **same** determinant reproduces, as
special cases:
- `α→0` ⇒ Neumann cap ⇒ the D/N mouth DtN `−(ω/c_s)tan(ωL₀/c_s)`, half-shifted ladder `(n+½)πc_s/L₀`;
- `α→∞` ⇒ Dirichlet cap ⇒ the D/D mouth DtN `+(ω/c_s)cot(ωL₀/c_s)` (sign/branch differs), **integer** ladder
  `nπc_s/L₀` with `n≥1` (no `n=0`); [N/N, for completeness, shares the integer ladder but **includes** the
  `n=0` zero mode];
- one **auditor-chosen numeric `α`** ⇒ a Robin ladder that is **neither** family (poles solve the
  `α`-dependent transcendental condition), confirming the `α`-dependence is real.

Emit `counterfactual_guard: {robin_determinant_emitted: true, recovers_DN_at_alpha0: true,
recovers_DD_at_alpha_inf: true, halfshift_destroyed_for_DD: true, numeric_alpha_distinct: true,
dtn_mismatch: true}`. If one fixed determinant cannot reproduce both limiting families, or the `(n+½)` ladder
survives the `α→∞` limit ⇒ `FAIL_COUNTERFACTUAL` (the answer was a lookup, not a solve).

**Relabel (physics):** changing the cap BC does **not** change the differential operator `L_s` (it stays
`d²/ds²+(ω/c_s)²`); it changes the **DtN / pole condition**. Hence the emitted flag is `dtn_mismatch`, not
`operator_mismatch`.

**Tell to watch for (smoking gun):** if the engine uses a genuine `dsolve`/transfer-matrix anywhere but the
headline `Z₀₀` is a typed-in `−(ω/c_s)tan(...)` read back by a solver-named wrapper, that is the same fake one
level deeper. The DtN must be *produced* from the BC-applied solution, with the solution emitted as an
artifact.

---

## §4. Anti-tautology firewall + execution enforcement

- **FORBIDDEN:** writing `Z₀₀=−(ω/c_s)tan(ωL₀/c_s)`, the pole ladder, or the reduced operator **by hand** and
  reading it back. Every headline quantity must be **emitted from a real solve** (the `dsolve` solution, the
  BC-applied determinant/denominator, the transcendental pole equation, the small-ω series).
- **EMIT SOLVE ARTIFACTS** so the audit can verify a real solve: the general ODE solution, the BC-applied
  coefficient relations, the symbolic DtN before simplification, the denominator whose zeros are the poles,
  and the counterfactual's mismatched DtN / pole condition (the differential operator `L_s` is BC-independent
  — it is the boundary response that differs).
- **Controls must be able to fire from computed inputs** — no literal `pass=True`, no `x−x` identities, no
  spectral-class lookups standing in for a solve.
- **Reduction first, comparison second.** `L_s` (§1/§2.1) is produced by the reduction; the idealized
  operator is the *target* it is compared against, not the starting point.

---

## §5. Dimensional consistency (standing pre-numbers step)

Before any pole-ladder numerics, run a **units-restored** homogeneity check (do not trust `c_s=ħ=m=1` pins):
confirm `ωL₀/c_s` dimensionless, `Z₀₀` has dimension `ω/c_s = 1/length`, `c_s²=5Kρ⁴/m` dimensionally consistent
with the EOS, and the pole ladder `ω_n` has dimension `1/time`. Implementation may use `sympy.physics.units`
**or** explicit `{M,L,T}` dimension-tuple arithmetic — either is acceptable provided the check is
non-tautological (it must FAIL if a dimension is perturbed). (memory `feedback-dimensional-consistency-check`.)

---

## §6. Deliverables

- `software/stage1_solver/tools/pathA_30_dn_unit_test_sympy.py` and `..._dn_unit_test.wl` (dual-engine,
  independent routes; each `timeout 600`; exit 0).
- `software/stage1_solver/reports/pathA_30_dn_unit_test.md` — line 1 = the verdict rung; then the emitted
  reduced operator `L_s`, the DtN derivation, the pole ladder, BC provenance classification, the static-limit
  series, the round-trip closure, the counterfactual result, and `engine_agreement`.
- `software/stage1_solver/reports/pathA_30_results.yaml` — machine-readable: `verdict` (one of the §0 ladder
  rungs); `reduction_certificate: {background, projection, Q_fate, rho0_fate, cs_of_s_fate, Vconf0_fate}`;
  `operator_is_helmholtz`; `speed_is_cs`; `domain_is_L0`; `bc_type`; `bc_provenance: derived|imposed` plus
  `bc_derivation_emitted: true|false` (`false` ⇒ forces `BC_DEPENDENT`); `dtn` (string of the *derived*
  expression); `pole_ladder`; `halfshift: true`; `static_limit_series`; `round_trip: {r_L, R_rt, phi0}`;
  `counterfactual_guard` (the Robin block of §3); `engine_agreement`; `dim_check: pass`.
- A short feed note appended/linked under `research/pde_ledger/notes/stages/` recording the result against
  scaffold §9/§11.1 (do **not** edit the audited ledger stages).

**YAML/markdown only for anything an LLM reads** — no JSON (memory `feedback-no-json-for-llm-io`).

---

## §7. Explicitly deferred / out of scope (record, do not silently assume away)

Defer to later phases, but the script must **flag** any that intrude on the reduction rather than dropping
them:
- variable-coefficient corrections from `ρ₀(s)`, `c_s(s)`, the quantum potential `Q`, and `δV_conf` (Phase 2,
  full finite-throat linearized bulk problem) — **but only as explicitly parameterized higher-phase terms
  that the reduction certificate (§1) shows vanish at this order**; an *unsuppressed* intrusion is
  `FAIL_OPERATOR_INTRUSION`, not a deferral;
- coupling to the localized Maxwell mixed channels `A_w, J^w, F_{μw}` (kept alive later; scaffold §1.3 — do
  **not** zero them as ontology, just don't excite them in this frozen scalar unit test);
- the grouped real `P2` (ℓ=2) sector and the quadrupole normalization `m̂₀²P₀=54Gc_s⁵/5a⁵c⁵` (Phases 3–6);
- the brane↔bulk leak/return `S_leak`/`S_return` and the ℓ=0/1 residual (pathA_28/29; reconciliation note) —
  this is where the DtN acquires its **imaginary (radiative) part**; the frozen-wall DtN here is deliberately
  **lossless** (real-valued), the fixed-geometry Born–Oppenheimer response.

---

## §8. Process

- **Codex codes + runs both engines; Claude reviews only** (memory `feedback-claude-reviews-codex-codes`).
  The directive states requirements + acceptance criteria; it does **not** pre-design the script route beyond
  the "independent routes" and "emit artifacts" constraints.
- **Dual-engine required** (Mathematica can independently verify here): Codex needs
  `--sandbox danger-full-access` to run `math`; or the orchestrator runs `math` directly as arbiter.
  ≤2 concurrent `math -script` seats.
- **Review ordering:** iterate Codex (design-review, `model_reasoning_effort=xhigh`) to GREEN → one GLM pass →
  fold → Codex confirm to green → execute. Execution backgrounded; **never** wrap the Codex session in shell
  `timeout`.
- **Post-exec tri-review** (clean agents): orchestrator arbiter re-run (reproduce + refresh committed
  output) + transliteration-fidelity audit (code-vs-equations, term by term) + adversarial audit ("is the
  headline EARNED, or pass-by-construction?"). The adversarial leg is non-optional — fidelity alone cannot
  catch a faked "intended math" (pathA_29 lesson). Then user gate.
