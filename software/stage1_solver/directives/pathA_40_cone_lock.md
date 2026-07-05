# Directive pathA_40 — The speed cone-lock: is `λγ=1` / `c_E=c_γ` DERIVED, CALIBRATED, or a NO-GO? (consistency-knit gate 1)

**Status:** DRAFT v2 (Codex design-review + GLM-5.2 tertiary both folded — see changelog §10; pending Codex confirm-to-green → user gate
before execution). **Scope provenance:** `_scratch/consistency_knit_scope_v3.md` (`SOUND_AS_A_SCOPE`, Codex×3 + GLM folded) +
`_scratch/consistency_knit_sequencing_codex.md` (the Claude↔Codex sequencing conversation → settled option **(iv)**).
**Conceptual source (read first):** `docs/conceptual_foundation.md` §3/§4 ⭐ v8; `notes/consistency_knit_handoff.md`; `STATUS.md`
▶ RESUME HERE.

**This is the FIRST gate of the consistency knit** (do the four earned sectors live in ONE self-consistent parameter set?). It tests
ONLY the speed cone-lock + its provenance. NG5 (`SECOND_MEDIUM_DRIFT`) and the `pde_ledger` assembly are **follow-on directives** that
consume this gate's output (§9).

---

## 0. What this gate is worth — read BEFORE executing (honesty preregistration; GLM Finding 5)

**The cone-lock verdict is EXPECTED to be near-mechanically `CONE_LOCK_CALIBRATED` — it returns the honest prior by construction, not
by discovery.** On the earned ledger, with the moduli free, a rank/codimension increase `Δr=2` and "not derived" both fall out
mechanically, and `NO_GO` cannot fire. This gate is therefore NOT where the program's falsification power lives. Its genuine, non-vacuous
outputs are:
1. the **derivability pre-pass** (§3A) — the ONE cheap place a real surprise (a forced lock ⇒ `CONE_LOCK_DERIVED`) could appear, and,
   failing that, a *graded* honest negative that says WHERE a future derivation would have to come from;
2. the **route checked-negatives** (§3A) — computed, falsifiable, cited-to-report, not "not found";
3. the **coupled-pole off-cone residual** `Δ_pole ∝ C_hu²` (§3C) — a real computation in the gate's own parameter set;
4. an **algebraic-compatibility certificate** + the **parameter-freedom caveat list** feeding NG5.

The substantive falsification pressure lives DOWNSTREAM — NG5's one-medium drift count, the coupled-pole test (task #110), and the
held-out predictions (g−2, 5PN, ringdown, multi-defect). **Do not dress a mechanical `CALIBRATED` as a hard-won result.** Symmetric
scrutiny (binding): a `CONE_LOCK_CALIBRATED, Δr=2` verdict earns the SAME extra scrutiny as a surprising `DERIVED` — the gate must SHOW
the two locks are algebraically independent on the ledger (that `Δr=2` is genuine, not merely presumed), not assume it.

---

## 1. The earned ledger `E(θ)` (imported verbatim — do NOT re-derive; provenance-tagged)

Shared parameter vector `θ = (ρ, ρ_br, ρ_B0, χ_c, μ_R, K, m, M_h, c_E, C_hu, K_h, B_eff, Z, Q_E, b, ℓ, …)`.

**Earned equalities** (each source-tagged):
- `c_s² = 5Kρ⁴/m` — gravity/bulk GNLS EOS `P=Kρ⁵` (∝ρ²).
- `c_γ² = μ_R/ρ_br` — light transverse MacCullagh shear (`pathA_36`; Stage-4 `u_T1,u_T2` block). **`μ_R` "MacCullagh" is a borrowed NAME
  for a POSTULATED native brane stiffness** (`pathA_35` G0 tags `μ_R, ρ_br` = `postulated-ingredient`), not an imported constitutive law.
- `c_L² = B_eff/ρ_br`, `B_eff = ρ_B0²/χ_c` — on-brane density mode (`pathA_36`; Stage-4 `u_L` block).
- `c_h² = K_h/M_h = c_E²`, `K_h = M_h c_E²` — embedding-Goldstone/branon; `c_E` = the `pathA_38` dynamic Green **speed**. The Green FORM
  `exp(iRω/c_E)/(4πR)` is IMPORTED; this gate consumes only the speed `c_E`, so the form does not bind (tag: "speed earned; Green form
  imported").
- Stage-4 scalar 2×2 block over `(u_L,h)`: `M(ω,k)=[[ρ_br ω²−B_eff k², −C_hu k²],[−C_hu k², M_h ω²−K_h k²]]` — actual poles = mixed
  roots of `det M=0`, NOT the bare `c_L², c_h²` (§3C).

**Earned inequalities / signs:** `ρ_br,ρ,K,m>0` ⇒ `c_s²>0`; `μ_R>0` ⇒ `c_γ²>0`; `χ_c>0, ρ_B0≠0` ⇒ `B_eff>0`; `M_h>0, c_E²>0` ⇒
`K_h>0`; **`C_hu² < B_eff K_h`** (Stage-4 scalar stability — STRICT; encode via a positive slack variable, never `≤`); `Z<0` (gravity
drain premise, `pathA_29`); `Q_E≠0, b>0, ℓ>0` only where a field-content claim uses `q_h=2Q_E tanh(b/ℓ)/b`.

**Branch / provenance tags (do NOT silently swap):** real `pathA_36` branch = `SECOND_CLASS_PAIR` (one longitudinal DOF), NOT the tuned
Maxwell locus (`B_eff=0`, `K_θ=J²ρ_B0²/ρ_br`); its sign fact `K_θ=−κ_phase<0` (vs tuned `K_θ>0`). Real `pathA_39` field content =
`FIELD_SCALAR_VECTOR_DEPARTURE` (4 DOF). `μ_R,ρ_br,c_E,K` **independent unless a relation is explicitly earned** (§3A). The older
`λγ`/`BETA_GENUINE_GAP` gap (phrased `β_bulk_to_brane`) and the current `μ_R/ρ_br` phrasing are the **SAME** missing
speed-normalization relation — ONE gap (this consolidation prevents a double-count rig).

**Ledger housekeeping (list-or-exclude — nothing silently assumed):** `J`, `K_θ`, `κ_phase` (all appear only in the tuned-Maxwell branch
tag `K_θ=−κ_phase=J²ρ_B0²/ρ_br`) — out of θ / out of scope, since the tuned locus is not used (GLM Finding E1). `P0, g_mhat, g_G` (the
`pathA_22b` `P0·χ_Q·g_mhat²·λγ⁵/g_G=54/5` anchor) — downstream context, NOT cone-lock inputs. `χ_Q` — out of θ by design (the two-`χ_Q`
open item is the pde_ledger step's, §9).

---

## 2. The cone-lock crux

Full Lorentz/Maxwell needs BOTH:
- **(A) `λγ=c_γ/c_s=1`** ⟺ `L_A: μ_R/ρ_br − 5Kρ⁴/m = 0` (denominator-cleared `L_A': m·μ_R − 5Kρ⁴·ρ_br = 0`). Light cone = gravity cone.
- **(B) `c_E=c_γ`** ⟺ `L_B: c_E² − μ_R/ρ_br = 0` (cleared `L_B': c_E²·ρ_br − μ_R = 0`). The **Stage-2 Lorentz-form principal-cone
  diagnostic** (residual `(c_E²ρ_br−μ_R)/ρ_br`), the retired `η_T=1`; locks the **bare/principal** h speed to the transverse speed — NOT
  the coupled poles (§3C).

**Decisive question:** are (A),(B) **derivable** (forced once all four sectors share one microscopic medium action), **calibration
inputs**, or **incompatible** with `E(θ)` (a NO-GO)?

---

## 3. The gate (settled option (iv)): pre-pass → lean cone-lock → (NG5 follow-on)

### 3A. Derivability / freedom PRE-PASS (able-to-fail, graded — where a surprise WOULD appear if the ledger contained a throat integral)

Run FIRST. It inventories the earned reports for the two forcing routes and the freedom assumptions, and is able to return a genuine
`WELL_POSED` (⇒ a real derivation exists and this gate must HALT and promote it) — so it is not baked to a negative. (Honesty, GLM
Finding 5: a `DERIVED` surprise requires R1–R5 to be present, which — per the current earned ledger, where deriving `μ_R` as a
bulk-defect integral is not well-posed — they are not; on the current ledger this pre-pass is expected to return the graded negative.)

- **Route A — brane-as-bulk-defect constitutive integral for `μ_R`.** Question: does the earned ledger already contain the object a
  forcing derivation would integrate? Decide by a **concrete source-fact inventory** (each present/absent with a cited report location —
  no "not found" by taste; Codex design-review Finding 3):
  - **R1:** a nonlinear bulk/throat/interior action containing BOTH the GNLS bulk EOS and the in-plane brane shear displacement `u^a`.
  - **R2:** a transverse profile `f_u(w)` for the IN-PLANE shear mode (NOT the `h`-Goldstone profile).
  - **R3:** a common normalization/integral deriving BOTH `ρ_br` and `μ_R` from the same bulk/throat object.
  - **R4:** an algebraic reduction proving `μ_R/ρ_br = 5Kρ⁴/m`.
  - **R5:** no residual free geometric factor (`ℓ`, thickness, profile norm, compactness). **R5 is a NEGATIVE condition gated on R4 (GLM
    Finding 2 / R5-ambiguity):** R5 = `present` ONLY if R4 is `present` AND the R4 reduction leaves no free factor; if R4 is absent, R5 =
    `not_applicable` (a missing derivation has no residual factor to audit — do NOT mark R5 `present` vacuously, which would inflate the
    inventory toward `WELL_POSED`).

  Classify by the **strongest inventory state**:
  - `ROUTE_A_WELL_POSED` — R1–R5 all present ⇒ **HALT the cone verdict** (this gate does not classify a ledger a real derivation would
    modify); spawn/require a dedicated derivation branch (fresh directive), then re-run cone-lock on the updated `E(θ)`.
  - `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` — a candidate throat/defect bridge exists or is explicitly deferred for
    brane/shear derivation, but ≥1 of R1–R5 is missing ⇒ emit the **missing-object list** (records WHERE a future `λγ` derivation must
    come from; ties to `project-simulation-deferred-complete-pde-strategy`).
  - `ROUTE_A_ABSENT` — there is NO candidate bridge at all for deriving `μ_R/ρ_br` from the bulk/throat; the relevant reports only tag
    `μ_R, ρ_br` as independent/postulated.

  **Convention (binding — pick one so the labels are exclusive, per Finding 3):** `pathA_38`'s `sech²/ℓ` h-throat + the deferred
  nonlinear-throat program **ARE treated as a candidate bridge** for a future brane/shear derivation. Therefore the expected landing on
  the current files is `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` with missing `{R1,R2,R3,R4,R5}` as applicable (R2 missing: the
  `sech²/ℓ` profile is the `h`-Goldstone, not `f_u`; R1/R3/R4 missing: `μ_R,ρ_br` are `postulated-ingredient` in `pathA_35` G0, no bulk
  integral). `ROUTE_A_ABSENT` is reserved for the counterfactual where the inventory is restricted to `pathA_35/36` and finds no bridge
  attempt at all. The pre-pass MUST be constructed so `ROUTE_A_WELL_POSED` is reachable if R1–R5 were present (able-to-fail — the §5
  WELL_POSED control injects a synthetic ledger carrying R1–R5 as SOURCE FACTS and REQUIRES `WELL_POSED`/`HALT`, recomputed from those
  facts by the SAME inventory predicates). **Provenance of the `sech²/ℓ` form (GLM Finding E2):** tag it — the `sech²` kink profile is a
  standard soliton FORM (imported), not native-derived structure; crediting it as a "candidate bridge" component counts only the
  *existence* of a throat program, not earned in-plane-shear structure. **`UNDERDETERMINED` does NOT count as earned structure for NG5
  (GLM Finding 2a):** the `UNDERDETERMINED_MISSING_NONLINEAR_THROAT` label names a DEFERRED (not-yet-executed) program — "throat absent,
  program deferred," not "throat present, incomplete." It may NOT be counted as earned structure in NG5's one-medium drift count; only
  the §3A parameter-freedom caveat list feeds NG5.
- **Route B — shared brane elastic tensor / field-distinctness.** Treated as a **closed checked-negative + over-import guard**: `h`
  (embedding Goldstone SCALAR, in the `(u_L,h)` block) and `u_T` (transverse VECTOR, separate block) are **different DOF** (settled by
  `pathA_36`/`pathA_39` Stage-4) — so `c_E=c_γ` is NOT near-trivially forced by field identity. A thin-plate "one 2D elastic tensor,
  `K_h∝μ_R·thickness²`" relation is **over-import** (standard thin-plate dressed as native) and may NOT be credited as a derivation;
  the guard flags any attempt to introduce it.
- **Freedom check (the sharpest rig risk — GLM Finding 3).** Inventory whether `C_hu` and `ρ_br` are FREE on the earned ledger or already
  tied by an earned relation: is `C_hu` earned-derived (e.g. from the charge residue `q_h=2Q_E tanh(b/ℓ)/b`) rather than free? is `ρ_br`
  = `ρ×(brane thickness)` rather than independent of `ρ`? Emit `FREEDOM_UNCONSTRAINED{C_hu, ρ_br}` (expected) or, if a tie is found,
  `FREEDOM_TIED{param, relation}` — a tie can flip the Layer-1 NO_GO result (§3B) and MUST be carried into the SAT domain.

### 3B. Lean cone-lock — two-layer decision (run only if the pre-pass returns ABSENT/UNDERDETERMINED, NOT `WELL_POSED`)

**Layer 1 — NO-GO (real-algebraic satisfiability).** Baseline sanity: `E(θ)` alone UNSAT → `NO_GO(sector-ledger)` (stop; not the lock's
fault). Then: is `E(θ) ∧ L_A ∧ L_B` satisfiable over the positivity/stability semialgebraic domain (with any §3A `FREEDOM_TIED`
relations appended)? UNSAT → `NO_GO(cone-lock)` + the minimal conflicting earned constraints (unsat core).

**Layer 2 — provenance, per lock, genericity-free primitives:**
- **Entailment (DERIVED?) via the REAL radical.** For each `L_i`, test whether the earned target-blind relations force it over the REAL
  semialgebraic set — **real-radical membership (RealNullstellensatz) / Positivstellensatz certificate**, NOT algebraic-radical over ℂ
  (positivity can block complex-only zeroes; complex radical = cheap necessary pre-filter only). Pass = `L_i` DERIVED. **If a genuine
  real-radical certificate is not mechanizable with available tooling, emit `ENTAILMENT_INCONCLUSIVE(L_i)` — do NOT let the second engine
  reassert the first engine's boolean** (Codex dual-engine realism note).
- **Calibration count via solution-variety CODIMENSION (Krull dim), not generic Jacobian rank.** `Δr = dim V(E_=) − dim
  V(E_= ∧ L_A ∧ L_B)` on the REAL locus. Generic Jacobian rank equals this only if `E_=` is a regular sequence — NOT guaranteed
  (syzygies / non-transverse intersections over 16+ params); use the Krull dimension of `ℝ[θ]/I` before vs after quotienting by
  `(L_A,L_B)`; generic-rank is a sufficient-condition proxy only. **Name both dropped assumptions** (regular-sequence; complex-vs-real)
  in the report. Guard that the dimension is the REAL-locus dimension, not a complex Gröbner dimension counting non-real components.

### 3C. Field-content overlay + the coupled-pole residual (first-class output)

Attach `FIELD_CONTENT = FIELD_SCALAR_VECTOR_DEPARTURE` (imported from `pathA_39` Stage-4, UNCHANGED — 4 physical DOF; the Maxwell
counterfactual reaches 2 only by REMOVING the h block, which equal speeds do NOT do). **Never** label the cone lock `MAXWELL_RESTORED`.

Compute + emit the coupled-pole residual. Evaluate `det M` at `ω²=c_γ²k²` under `c_E²=c_γ²=K_h/M_h`: the second diagonal entry vanishes,
giving `det M|cone = −C_hu² k⁴` — the coupled scalar poles are generically OFF the light cone, residual `∝ C_hu²`. Because the scalar
block has **two** mixed roots, emit BOTH (Codex design-review Finding 1): (i) the **determinant residual** `det M|_{ω²=c_γ²k²} = −C_hu²k⁴`
(the unambiguous primary diagnostic), and (ii) the **root pair** `{Δ_pole_−, Δ_pole_+}` where `Δ_pole_± = (v²_±(θ) − c_γ²)` and `v²_±` are
the two mixed-root phase speeds² from `det M=0` (solve for `ω²/k²`; report both, with the root-selection/normalization rule stated
explicitly). It lives entirely in this gate's parameter set and does NOT pull the #110 verdict in; it prevents a clean `CALIBRATED` from
being misread as "scalar sector is Lorentz-consistent." (B) locks only the bare/principal h speed; the full coupled-pole Lorentz question
(needs constraining `C_hu, q_L, B_eff/ρ_br`) is task #110.

---

## 4. Verdict grammar — a STRICT ordered decision table (mutually exclusive, complete; Codex design-review Finding 5)

Two facts drive the primary verdict: `Δr` (§3B codimension increase) AND, per lock, a **provenance status** ∈ {real-entailment PROOF
(`ENTAILED`) / real non-entailment WITNESS (`WITNESSED`) / `INCONCLUSIVE`}. `Δr` alone does NOT identify which lock is derived (a `Δr=1`
can be one *shared* calibration with neither lock entailed). Evaluate in THIS order; the first matching case is the verdict:

1. `NO_GO(sector-ledger)` — baseline `E(θ)` UNSAT (deeper than the lock; carries the unsat core).
2. `HALT_ROUTE_A_WELL_POSED` — the §3A pre-pass returns `ROUTE_A_WELL_POSED` (R1–R5 present); cone verdict deferred to a dedicated
   derivation branch.
3. `NO_GO(cone-lock)` — `E ∧ L_A ∧ L_B` UNSAT (proven); carries the minimal conflicting constraints (unsat core) **and, if a §3A
   `FREEDOM_TIED` relation caused the failure, that relation** (so the freedom diagnostic appears on the NO_GO, not only on successful
   verdicts). **SAT computational-failure fallback (GLM Finding 4):** if neither SAT nor UNSAT is *proven* over the semialgebraic domain
   (CAD blowup / timeout / unknown on either engine), this is **NOT** `NO_GO` — proceed to Layer 2 carrying a `SAT_UNCERTIFIED` rider; if
   Layer 2 also cannot certify, the verdict is `CONE_LOCK_PROVENANCE_INCONCLUSIVE` (step 9).
4. **PROVENANCE GUARD (closes the `SHARED_CALIBRATION` leak — GLM BLOCKER Finding 3):** if EITHER lock's provenance is `INCONCLUSIVE` (no
   entailment proof AND no non-entailment witness), OR the engines DISAGREE on a lock's status, OR `Δr` real-dimension certification
   failed, OR a `SAT_UNCERTIFIED`/dimension-uncertified state persists → the verdict is `CONE_LOCK_PROVENANCE_INCONCLUSIVE` (step 9). **No
   calibration verdict (`CONE_LOCK_CALIBRATED` or `CONE_LOCK_SHARED_CALIBRATION`) may absorb an `INCONCLUSIVE` lock.** Only once BOTH locks
   are definitively `ENTAILED` or `WITNESSED` do the `DERIVED`/`PARTIAL`/`CALIBRATED`/`SHARED_CALIBRATION` cases (steps 5–8) apply.
5. `CONE_LOCK_DERIVED` — BOTH `L_A,L_B` `ENTAILED` AND `Δr=0`. Too-clean → **extra scrutiny** (Gate-K prior) before banking.
6. `CONE_LOCK_PARTIAL(derived=A, calibrated=B)` — `L_A` `ENTAILED`, `L_B` `WITNESSED`, `Δr=1`; symmetric `(derived=B, calibrated=A)`.
7. `CONE_LOCK_CALIBRATED` — BOTH locks `WITNESSED` AND `Δr=2`. Expected (§0). **Triggers the SAME extra scrutiny as DERIVED:** SHOW the
   two locks are algebraically independent on the ledger (that `Δr=2` is genuine, not presumed) — else it is a generic-default artifact.
8. `CONE_LOCK_SHARED_CALIBRATION(Δr=1, derived=none)` — `Δr=1` with BOTH locks `WITNESSED` (neither entailed; `E` relates them so the two
   locks cost one shared dimension). **Requires both non-entailment witnesses** — an `INCONCLUSIVE` lock is caught by the step-4
   PROVENANCE GUARD and can never land here.
9. `CONE_LOCK_PROVENANCE_INCONCLUSIVE` — the terminal fallback for any state not resolved above (redundant with the step-4 PROVENANCE
   GUARD, kept as a catch-all).

**Rider attachment (GLM Finding E3 — two rider classes, do not conflate):** every verdict carries the §3A pre-pass result (Route-A grade
+ missing-object list; Route-B checked-negative; freedom result). Every non-`NO_GO`/non-`HALT` verdict additionally carries, atomically:
- **Computation-cited riders** (each MUST cite the producing computation; no decorative constants):
  - `coupled_scalar_poles: det_residual=−C_hu²k⁴, {Δ_pole_−,Δ_pole_+}; OFF_CONE_under_AB (∝C_hu²); status: OPEN(#110)` (← §3C);
  - `field_content: FIELD_SCALAR_VECTOR_DEPARTURE` (← imported `pathA_39` Stage-4);
  - `SAT_UNCERTIFIED` / `ENTAILMENT_INCONCLUSIVE(L_i)` where the certification actually failed.
- **Scope-boundary caveats** (fixed statements, no producing computation — labelled as such, not cited):
  - `[Δr is the LOCK-CONSTRAINT SLICE ONLY; full one-medium drift = NG5, expected larger]`;
  - `NO_GO_non_firing_conditional_on: {C_hu, ρ_br} freedom — NG5 must certify; re-run trigger if revoked`.

---

## 5. Controls (able-to-fail; `feedback-negative-verdict-short-circuit`)

**Discipline (BINDING — anti-rig; Codex design-review Finding 2):** every control fixture MUST mutate the input *equations / provenance
source-facts* and then **recompute** the verdict through the SAME predicates the production run uses. The result emitter MUST NOT accept
a direct verdict/flag input — no `route_a_verdict`, `forced_lock=True`, `freedom_tied=True`, `route_a_well_posed=True`, typed boolean, or
expected-label short-circuit anywhere in a fixture. A control that sets its own answer is a rig (the Gate-3/4/pathA_38 pattern).

**The verdict-flipping counterfactuals (each MUST flip the verdict, recomputed from mutated source facts):**
- **WELL_POSED control (Route-A yes-path):** inject a synthetic ledger carrying R1–R5 (§3A) as SOURCE FACTS (the throat/interior action,
  the `f_u` profile, the common normalization, the `μ_R/ρ_br=5Kρ⁴/m` reduction, no residual factor) → the route inventory MUST recompute
  `ROUTE_A_WELL_POSED` / `HALT` from those facts.
- **ABSENT control (Route-A no-bridge path; GLM Finding 1A):** restrict the inventory to a synthetic ledger containing ONLY the
  `pathA_35/36` postulated-modulus facts (no throat/defect bridge of any kind) → MUST recompute `ROUTE_A_ABSENT`. Without this, an
  implementation that always returns `UNDERDETERMINED` would pass every other control — this control proves the discriminator actually
  distinguishes "candidate bridge present" from "no bridge."
- **PARTIAL-inventory control (missing-object list; GLM Finding 1C):** inject R1,R2 present + R3,R4,R5 absent → MUST recompute
  `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` with the SPECIFIC missing set `{R3,R4,R5}` (verifies the missing-object list is
  computed from the inventory, not hardcoded).
- **Forced-lock counterfactual (GLM Finding 1B — must exercise the real entailment path, not pass-by-construction):** append a
  target-blind FAKE earned constitutive relation that is **syntactically DISTINCT from `L_A,L_B`** and **IMPLIES both** (the implication
  verified by the entailment engine, and the relation kept simple enough that a real-radical/Positivstellensatz certificate is
  mechanizable) → MUST recompute `Δr=0` / `CONE_LOCK_DERIVED`. **The fixture may NOT append `L_A,L_B` (or `L_A',L_B'`) themselves** — that
  would make `Δr=0` and entailment trivial by construction and would test nothing.
- **A-only partial counterfactual:** append a target-blind earned relation (syntactically distinct from `L_A`, implying it) forcing `L_A`
  while leaving `L_B` free → MUST recompute `Δr=1` / `CONE_LOCK_PARTIAL(derived=A, calibrated=B)` (exercises the partial grammar so the
  implementation cannot only recognize DERIVED/CALIBRATED/NO_GO). Add the symmetric **B-only** counterfactual (→ `derived=B`) if cheap.
- **Over-constrained counterfactual:** append a real semialgebraic contradiction (e.g. an earned relation pushing `C_hu` above
  `√(B_eff K_h)`) → MUST recompute `NO_GO` (Layer-1 can fail).
- **Freedom-tie control:** append `C_hu = f(q_h)` at a value that, with the lock, violates `C_hu²<B_eff K_h`, into the SAT domain → MUST
  recompute `NO_GO` with the tie relation in the diagnostic (proves the §3A freedom result actually feeds Layer-1).

**Baseline PRODUCTION check (NOT a control — it does not flip anything):** the real earned `θ` run. It MUST compute and report the
verdict + the codimension increase `Δr` (not merely "SAT"). If the computed result differs from the §0 prior, **preserve the result and
trigger extra scrutiny** — it is NOT an expected-landing assertion and a non-`CALIBRATED` computed result is NOT a failure of this check.

---

## 6. Dual-engine (BINDING; `feedback-dual-engine-required`, Stage-4 lesson)

Every computed claim (SAT satisfiability, real-radical entailment, Krull-dimension `Δr`, the `Δ_pole` residual, each control) is
**dual-engine**: SymPy + a genuinely INDEPENDENT Mathematica route, engine-agreement asserted per quantity. **The second engine must
DERIVE, not echo the first engine's boolean** (the Stage-4 vacuous-reassertion pattern). The independence must live in the METHOD (Codex
design-review Finding 4):
- **SAT / NO_GO:** SymPy builds the polynomial/slack system and searches exact witnesses or algebraic contradictions; Mathematica uses
  `Resolve`/`Reduce` over `Reals` on the semialgebraic formula. Mathematica MUST NOT consume SymPy's SAT boolean.
- **Non-entailment (NOT_DERIVED):** either engine proves it by exhibiting a REAL witness/cell for `E ∧ (L_i ≠ 0)`. Lack of a real-radical
  certificate is NOT sufficient to call a lock non-derived.
- **Positive entailment (DERIVED):** requires a real UNIVERSAL proof/certificate — `ForAll[θ, E ⟹ L_i==0]` via Mathematica
  `Resolve`, or an explicit Positivstellensatz / real-radical certificate. If neither a proof nor a counterexample is obtained, emit
  `ENTAILMENT_INCONCLUSIVE(L_i)` (→ §4 step-4 PROVENANCE GUARD → `CONE_LOCK_PROVENANCE_INCONCLUSIVE`) — do NOT let either engine reassert
  the other.
- **`Δr`:** SymPy computes algebraic dimension via Gröbner/Hilbert machinery + real-feasibility checks on components; Mathematica
  computes semialgebraic/CAD cell dimension or an independent elimination/CAD dimension — NOT a repeated generic rank or a translated
  Gröbner answer. Guard that the dimension reported is the REAL-locus dimension, not a complex Gröbner dimension counting non-real
  components.
- **`Δ_pole`:** both engines independently derive the determinant/root residual from the 2×2 scalar block.

The orchestrator independently re-runs both engines as arbiter + a transliteration-fidelity audit + a clean adversarial-with-ablation
agent (per `docs/development_pipeline.md`).

---

## 7. Deliverables

- `directives/pathA_40_cone_lock.md` (this file, post-gauntlet).
- `tools/pathA_40_cone_lock_{sympy.py,.wl}` — dual-engine (pre-pass inventory logic + the two-layer decision + `Δ_pole` + all §5
  controls; both exit 0; engine-agreement asserted).
- `reports/pathA_40_cone_lock.md` (verdict on line 1) + `reports/pathA_40_cone_lock_results.yaml` — the ledger `E(θ)`, the pre-pass
  graded verdicts + cited report lines, the two-layer result (unsat core / entailment certificates / `Δr` with the two named
  assumptions), the `Δ_pole` residual, the controls, the parameter-freedom caveat list for NG5, provenance.

---

## 8. Review plan (the gauntlet; `docs/development_pipeline.md`)

1. **This directive** → Codex design-review (xhigh) → fold → GLM tertiary → Codex confirm-to-green → **user gate** before execution.
2. **Execution:** Codex codes + runs + iterates to exit 0 + dual-engine; Claude reviews only.
3. **After the run:** orchestrator arbiter re-run (both engines) + transliteration-fidelity audit + clean adversarial-with-ablation
   agent (the §5 verdict-flipping counterfactuals — WELL_POSED, ABSENT, PARTIAL-inventory, forced-lock, A-only/B-only partial,
   over-constrained, freedom-tie, each recomputed from mutated source facts with no direct-boolean short-circuit — are the adversarial
   focus) → verdict → **user gate**.
4. **Never alter the calibrated process unilaterally** (`feedback-never-alter-calibrated-process`): if a step fails, HALT + bring to user.

---

## 9. Scope boundary + honest priors

- **This directive = cone-lock classification + field-content overlay ONLY.** Follow-ons: (2) **NG directive** — the one-medium
  drift/counting gauntlet (NG5 `SECOND_MEDIUM_DRIFT`, count every independent input across the four sectors) over the earned param set,
  consuming this gate's parameter-freedom caveat list; (3) **`pde_ledger` directive** — assemble the accepted four-sector spec +
  cone-lock/NG verdicts into `research/pde_ledger/`, resolving the two-`χ_Q` open item (`pathA_22b ≈0.712` vs `pathA_33 =1` — different
  computations, do not silently merge).
- **Honest priors (recorded before executing):** cone-lock = `CONE_LOCK_CALIBRATED, Δr=2` (mechanical; §0). Pre-pass Route A =
  `ROUTE_A_UNDERDETERMINED_MISSING_NONLINEAR_THROAT` with missing `{R1,R2,R3,R4,R5}` (§3A convention). Route B = closed checked-negative.
  `det M|cone = −C_hu²k⁴ ≠ 0` (scalar sector off-cone, departure persists → #110). A `DERIVED` or `NO_GO` would be a genuine surprise →
  extra scrutiny either way, not celebration. The gate's value is the checked-negatives + the residual + the freedom caveat + the
  compatibility certificate, NOT a surprise. **Calibrate-predict honesty** (`feedback-calibrate-predict-methodology`): a calibrated lock
  ABSORBS the EM/`λγ` anchor (zero surplus from the lock); report the earned-vs-calibrated split plainly, never dress the absorption as a
  prediction; the surplus lives in the held-out predictions.

---

## 10. Changelog

- **v1 → v2 (GLM-5.2 tertiary review folded).** BLOCKER: the §4 `CONE_LOCK_SHARED_CALIBRATION` exclusivity leak — a calibration verdict
  could absorb an `INCONCLUSIVE` lock, the exact silent-calibration failure the table claims to prevent; fixed with an explicit
  **PROVENANCE GUARD** (step 4) that routes ANY inconclusive/engine-disagreement/uncertified state to `PROVENANCE_INCONCLUSIVE` before any
  calibration verdict, and `SHARED_CALIBRATION` now requires BOTH locks `WITNESSED` (§4). NITs: (1A) added an ABSENT counterfactual
  control (restrict inventory to pathA_35/36 → require `ROUTE_A_ABSENT`) so an always-`UNDERDETERMINED` impl can't pass; (1B) the
  forced-lock control must be syntactically DISTINCT from `L_A,L_B` and IMPLY them (entailment-verified, Positivstellensatz-mechanizable)
  so it exercises the real proof path, not pass-by-construction; (1C) added a PARTIAL-inventory control (R1,R2 present + R3–R5 absent →
  specific missing set `{R3,R4,R5}`) so the missing-object list is computed not hardcoded (§5); (2a) `UNDERDETERMINED` does NOT count as
  earned structure for NG5 (deferred program, not present throat) + (2b) R5 is a negative condition gated on R4 (`not_applicable` if R4
  absent, never vacuously `present`) (§3A); (4) SAT computational-failure fallback (`SAT_UNCERTIFIED` → not `NO_GO`, proceed to Layer 2 /
  `CONE_LOCK_PROVENANCE_INCONCLUSIVE`) (§4); (5) §3A header oversell reworded; (E1) `K_θ,κ_phase` list-or-excluded (§1); (E2) `sech²/ℓ` form provenance tagged as
  imported soliton (§3A); (E3) riders split into computation-cited vs scope-boundary-caveat classes (§4).
- **v0 → v1 (Codex design-review xhigh folded).** BLOCKERS: (F3) Route-A grading made operationally decidable via a concrete R1–R5
  source-fact inventory + a binding convention (`pathA_38` h-throat + deferred nonlinear-throat = candidate bridge → expected
  `UNDERDETERMINED` with the missing-object list; `ABSENT` reserved for the pathA_35/36-only counterfactual) (§3A); (F2) controls
  rig-locked — every fixture must mutate source facts and recompute (no direct verdict/flag inputs), the misnamed "baseline" demoted to a
  non-flipping production check, and the missing A-only/B-only PARTIAL counterfactuals added (§5); (F5) verdict grammar replaced by a
  strict ordered decision table (8 cases, mutually exclusive/complete) with `CONE_LOCK_SHARED_CALIBRATION` and
  `CONE_LOCK_PROVENANCE_INCONCLUSIVE` added, `CALIBRATED` disallowed without a real non-entailment witness per lock, and rider-attachment
  rules incl. the freedom-tie/unsat-core on `NO_GO` (§4). NITS: (F1) `Δ_pole` split into the determinant residual + the `{Δ_pole_−,
  Δ_pole_+}` root pair with a stated root-selection rule (§3C); (F4) dual-engine method-divergence pinned per computation (SAT /
  non-entailment / positive-entailment / `Δr` / `Δ_pole`) so Mathematica cannot echo SymPy's booleans (§6); (F6) baseline expected-landing
  phrasing removed, riders required to cite their producing computation (§4/§5).
