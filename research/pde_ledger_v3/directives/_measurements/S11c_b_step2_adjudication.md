# S11c-b step 2 — per-case adjudication (WIP; user chose FULL RIGOR all families 2026-08-29)

Grounds: step-1 grade fingerprint (`scripts/S11c_b_grade_fingerprint.py` over the committed multigrade run;
`~/.s11_build/S11c_b_grade_fingerprint.txt`) + raw operands from the adjudication run
`~/.s11_build/S11c_b_adjud_fullres_run.out`. `(a,b)` = `(eta_bg^a, sigma_W^b)`; A=SymPy(PY), B=blind Wolfram(WL).
⛔ COMPUTE, do not assert; ⛔⛔ never blanket-collapse — the four families adjudicate DIFFERENTLY.

## 1. ADMISSIBILITY_OPERATOR_OPERAND (BODY_FORCE, DOF=THETA) ×4 — genuine WL under-retention (candidate repair)
A(PY)=`{(0,1),(1,1)}`, B(WL)=`∅`. PY carries a background-Laplacian body force
`κ_θ·σ_W/(L_W·W_0)·(η_bg·w1 − 1)·Σ_i ∂_{ii}w1` (a `∇²w1`, grade σ¹); WL's admissibility operand is identically 0.
This is the standalone §2b/§3d background-order balance (NOT a term-provenance split). §3a is explicit: "a second
spatial derivative of `W_bg` is still first order in background amplitude, `O(η)`/`O(σ_W)`, and is not dropped."
⇒ PY spec-correct; WL's admissibility path lacks the `∇²W_bg` body force entirely. ADJUDICATED — the FINDING (genuine WL
under-retention) is SETTLED (WL operand 0 while §3a mandates the term); the MECHANISM below is a leading hypothesis
from a WL CODE READ, NOT yet run-confirmed. WL's θ body force = `applyBackgroundProfileWithGeneratedJets[eulerScalar[firstVariation, thetaField]]`
= `truncateBackground[eulerScalar[...] /. profileRules]` (WL audit L1339-1340, L360). `truncateBackground` is the
`Series[etaBg,0,1] + Series[sigmaW,0,1]` truncation (`truncateScalar`, L102-125) — it keeps only FIRST background
order in each bookkeeper. A `∇²W_bg` (second spatial derivative) is graded by WL at SECOND background order and is
truncated to 0, whereas §3a mandates it be retained at first shape order (PY encodes it at `σ_W¹` with a dedicated
second-jet symbol `w1_profile_d1d1`, per §3a "not dropped"). ⇒ GENUINE WL under-retention; the repair OBJECT: WL's
admissibility (and, if the same truncation bites elsewhere, the operator) must retain the background-amplitude-order
`∇²W_bg` body force per §3a — i.e. `truncateBackground`/its grading must NOT count a spatial derivative of the profile
as adding a background-bookkeeper order. ⚠ CONFIRM at repair time whether the exact mechanism is the `truncateBackground`
Series truncation or a missing energy term (run WL, inspect `eulerScalar[firstVariation, thetaField]` pre- vs
post-`truncateBackground`); either way the fix restores the §3a first-shape-order retention. This is the CLEAR repair
candidate, common to both scope options.

## 2. SLAB_OPERATOR_TERM_ORIGINS / ADVECTIVE ×4 — representational provenance re-bucketing (pending row check)
Computed operands (LAB_HELD/RHO4):
- PY advective = `(rho_br/W_0)·Σ_i u_{i,t}·∂_i w1_profile` = `(rho_br/W_0)·∇W_bg·u_t`  [grade (0,1)].
- WL advective = `rho_br·(1+eta_bg·w1)·Σ_i ∂_i u_{i,t}` = `(rho_br/W_0)·W_bg·∇·u_t` [grades (0,0),(1,0)] PLUS the
  same `∇W_bg·u_t` (0,1) term.
⇒ WL's ADVECTIVE bucket = the full divergence `∇·(W_bg u_t) = ∇W_bg·u_t + W_bg·∇·u_t`; PY's ADVECTIVE bucket = only
`∇W_bg·u_t`. The residual `{(0,0),(1,0)}` is the `W_bg·∇·u_t` (accumulation/compressibility) part PY assigns to a
DIFFERENT provenance bucket. The two engines use DIFFERENT bucket names: term-origin OBJECTs route
`ADVECTIVE→FLAG`, `KINETIC→FLAG`, `BULK_ENERGY→STRUCTURE_INCOMPLETE`, and `ACCUMULATION/FACE/FLUX/FACE_FLUX→COVERAGE`
(one-engine-only). ⇒ candidate REPRESENTATIONAL (same operator, different provenance bookkeeping) — CONFIRM by
checking the full `SLAB_OPERATOR` U_MOMENTUM row agrees modulo the ENERGY_BASIS quotient (§1d "modulo total in-plane
divergences"). ⚠ `PROTECTED_UNREDUCED` alone does NOT prove pure-representational (routing only tests presence of a
protected atom); the energy-basis quotient reconciliation (deferred/owed) is the decisive computation.

## 3. SLAB_OPERATOR_TERM_ORIGINS / KINETIC ×4 — likely genuine WL under-retention (pending row check)
THICKNESS_ROW: A(PY)=`{(0,0),(1,0),(2,0)}`, B(WL)=`{(0,0)}`; u-momentum rows AGREE (residual ∅).
PY thickness-kinetic = `μ·e_W_tt·W_bg²`; WL = `μ·e_W_tt·W_0²`. Residual `μ·e_W_tt·W_0²·(2η_bg·w1 + η_bg²·w1²)` =
`μ·e_W_tt·(W_bg²−W_0²)`. Direction PY-has-more (like admissibility, unlike advective). (1,0) is first-order (§2a
retain); (2,0)=η²w1² is the shape-order question (§3a: independent-factor count for `w1²`). ⇒ likely WL
under-retention (W_bg²→W_0² truncation) at least for (1,0); CONFIRM via the THICKNESS row modulo quotient + §2a/§3a
on (2,0).

## 4. COUPLING_KERNEL (SECTOR=TRANSVERSE_TO_THICKNESS, ±ADJOINTNESS_OPERAND_FORWARD) ×8 — certified non-IBP bulk; §2a which-engine verdict PENDING
Certified genuine non-IBP bulk by BOTH v2 build legs (from-scratch Euler operator + exact rational zero-test;
`REPRESENTATIONAL_DIVERGENCE=0` real). Grades: A(PY)=`{(0,1),(1,1),(2,1),(3,1),(4,1)}`+rem`(5,1)`;
B(WL)=`{(0,1),(0,2),(1,1),(1,2)}`. BIDIRECTIONAL — WL carries `b=2` cells `(0,2),(1,2)` PY lacks; PY carries `a≥2`
cells `(2,1),(3,1),(4,1),(5,1)…` WL lacks (the rational density `1/(1+η w1)` geometric tail); both differ at
`(0,1),(1,1)`. ⇒ the non-IBP BULK is certified (v2 legs) and the bidirectional grade pattern is MEASURED; the §2a
per-cell "which engine is spec-correct" verdict is PENDING. Adjudicate which retention is spec-correct per §3c
(weak variational restriction) + §2a (first shape order in EACH of η,σ): WL's `b=2` = second σ_W order (§2a would
truncate it → WL over-retains σ²?); PY's `a≥2` = higher η order from the density (§2a first η order → PY over-retains
η²?). Both may be over-retaining beyond §2a's "first shape order in each" — OR the density's full η-dependence is
mandated. OPEN — needs the §2a/§3a per-cell reading + possibly a WL/PY correctness check.

## Plan (full rigor)
- (A) ADMISSIBILITY — SETTLED (genuine WL under-retention; §3a/§3d mandate the ∇²W_bg body force; WL operand 0).
  Mechanism = leading-hypothesis `truncateBackground` truncation, run-confirm at repair. No further step-2 work — it
  is a step-3 repair target.
- (B) SLAB advective/kinetic: the ENERGY_BASIS quotient reconciliation — does each full operator ROW (U_MOMENTUM for
  advective, THICKNESS for kinetic) agree modulo the §1d total-in-plane-divergence quotient? ⛔⛔ DO NOT reuse
  `classify_total_divergence` on the SLAB rows — the v2 build legs established that classifier is valid ONLY for
  WEAK-PAIRING DENSITIES, NOT strong operators (`∂L/∂q − Div(...)` already carry a Div; testing them for total
  divergence conflates the operator's own Div with representational freedom). The strong-operator representational
  check is a VARIATIONAL-EQUIVALENCE / energy-basis-quotient test: two energy bases differing by a total in-plane
  divergence yield the SAME bulk operator (EL derivative of a total divergence = 0). Reconcile the 07/10-vs-08/11
  gamma-DivGrad quotient at the ENERGY level (the deferred/owed item) and verify each row residual is the operator
  image of an energy-quotient difference (⇒ representational) vs a genuine bulk operator difference. This is the hard
  deferred sub-project; reviewed instrument → legs. ⚠ direction already computed: advective = WL-has-more (PY
  re-buckets W_bg∇·u_t to ACCUMULATION/FLUX ⇒ likely representational); kinetic = PY-has-more W_bg²vsW_0² (⇒ likely
  genuine WL under-retention, NOT re-bucketing).
- (C) COUPLING: §2a/§3a per-cell adjudication of the bidirectional grades → which engine (or genuine disagreement).
- Then step 3 repairs ONLY the confirmed genuine engine gaps (⇒ admissibility ∇²w1 [SETTLED], kinetic (1,0) if
  confirmed genuine; NOT advective if representational), step 4 re-run + re-adjudicate, step 5 honest record.
- ⛔ a one-engine fix is a SPEC question FIRST — if §3d/§2a is ambiguous for a cell, that is a SPEC defect (fix spec),
  not an engine bug.

## 5. REFINED FRAMING (2026-08-29, code-grounded, post-compact) — the requested-truncation ROW-LEVEL residual instrument

Commands run this session (rule 2 — the claims carry the reads):
- `sed -n '95,135p' mathematica/S11c_b_brane_operator_mathematica_audit.wl` → `truncateScalar` = `Normal[Series[·,{etaBg,0,1}]]` then `Normal[Series[·,{sigmaW,0,1}]]`, PROTECTING `Inactive[_][…]` operator tokens across the Series (products of protected jets survive to higher combined order). WL truncates INTERNALLY at (η≤1 ∧ σ_W≤1).
- `sed -n '330,370p' …audit.wl` (profileRules): `Derivative[i,j,k][widthBase] -> sigmaW·profileJetSymbol["W_BG",orders]/LWidth^(i+j+k-1)`. ANY derivative order carries EXACTLY ONE `sigmaW` (extra spatial order absorbed into `LWidth^(order-1)`). ⇒ WL grades ∇²W_bg at σ_W¹ — §3a-consistent.
- `grep -n … scripts/S11c_b_brane_operator_sympy_audit.py`: PY `first_shape_series` (L713-726) = `sp.series(·,eta_bg,0,2).removeO()` + keep only `powers[eta_bg]≤1 and powers[sigma_W]≤1`. PY second jet L1865/L2117 = `sigma_W·w_profile_second/L_W` (σ_W¹). PY emits the COUPLING kernel at FULL η dependency (the raw `1/(1+η w1)` tail, η⁰…⁵) — first_shape_series is NOT applied to that operand.
- PY term-origin buckets (L1574-1583): `BULK_ENERGY`, `KINETIC`, `FACE_FLUX`, `ADVECTIVE`, `DIVERGENCE_FLUX`.

CONSEQUENCES (refine §§1-4 above):
- BOTH engines grade a 2nd spatial derivative of W_bg at σ_W¹ (first background order) — they AGREE on the jet-grading convention. ⇒ ADMISSIBILITY WL=∅ is NOT a `truncateBackground` truncation of a correctly-generated ∇²w1; it is a MISSING ENERGY TERM / missing variation in WL's §3a basis (the θ-variation `eulerScalar[firstVariation,thetaField]` produces no ∇²W_bg body force). The FINDING (genuine WL under-retention, §3a/§3d mandate it) is unchanged and still SETTLED; the MECHANISM is refined toward a missing §3a invariant (N15) rather than a Series truncation. ⚠ still CONFIRM at repair by running WL.
- BOTH engines implement a (η≤1 ∧ σ_W≤1) truncation primitive independently (PY `first_shape_series`, WL `truncateScalar`) ⇒ strong corroboration that the §2a requested truncation = keep (a,b) with a≤1 ∧ b≤1 (drop any η² or σ_W²), ε to first order. ⚠ two blind engines implementing the same reading is evidence, NOT proof (rule-7 shared blind spot possible) — the truncation reading is VERIFIED by the 2 decision-list legs on the instrument directive (they read §2a/§3a).
- The step-1 multigrade compares PER-BUCKET; it cannot see that PY's ADVECTIVE + ACCUMULATION/FLUX sums to WL's single ADVECTIVE bucket. The missing computation is the FULL-ROW STRONG-FORM residual.

⇒ INSTRUMENT (merges resume-prompt (A)+(B) into ONE rule-4 computation): per operator ROW × (anchoring, density), compute
  PY_row  := Σ_buckets PY term-origin operands of that row, with every Inactive in-plane divergence ACTIVATED/expanded to strong form (product rule ∇·(c F)=∇c·F + c∇·F),
  WL_row  := Σ_buckets WL term-origin operands of that row, likewise expanded,
  both TRUNCATED to the §2a requested truncation (ε¹, η≤1, σ_W≤1),
  RES_row := PY_row − WL_row     (emit all three; rule 2 — print, do not assert).
Verdict reads off RES_row:
  RES_row = 0  ⟺ engines AGREE within the deliverable scope (advective's re-bucketing / divergence-form arrangement is REPRESENTATIONAL).
  RES_row ≠ 0  ⟺ genuine in-scope difference (missing/frozen operator content) → adjudicate which engine is spec-correct against §3a/§3b/§3c; that engine's sibling is the repair target.
Predicted (COMPUTE, do not assume): ADVECTIVE RES_row=0 (product-rule identity ∇·(W_bg u_t)=∇W_bg·u_t+W_bg∇·u_t, pure re-bucket); KINETIC RES_row≠0 at (1,0) (WL froze thickness inertia W_bg²→W_0², §3b "do not freeze a coefficient at its constant binding before differentiation" ⇒ genuine WL gap); COUPLING RES_row=0 (agree at (0,1),(1,1); the a≥2 / b≥2 tails are OUT of the requested truncation); ADMISSIBILITY RES_row=PY_row≠0 (WL_row≡0, genuine gap). The kinetic (2,0) and coupling η²/σ² cells are OUT of scope under (η≤1∧σ_W≤1) and are NOT engine bugs unless a leg reads §2a as mandating full (unlinearized) coefficients — a SPEC question the instrument's decision-list legs adjudicate.
This is the reviewed Phase-A/B build: directive (name the OBJECT = RES_row; WITHHOLD the predicted residuals, rule 5) → 2 decision-list legs (Codex + Grok; they also verify the (η≤1∧σ_W≤1) reading of §2a) → Codex build → 2 build legs (fresh Claude agent + Grok; FORM-ablate the row-sum and the divergence expansion) → commit.

## 6. WL-side computability + kinetic mechanism (2026-08-29, `sed -n '815,875p' …_mathematica_audit.wl`)

- WL emits the COMPLETE assembled operator per DOF directly, in ACTIVE strong form (not only reconstructable):
  `operator["U_MOMENTUM_ROWS"] = kineticU + rows["U_INTERNAL"] + faceRows["U_FACE"]` (L832),
  `["THICKNESS_ROW"] = kineticEw + rows["EW_INTERNAL"] + faceRows["EW_FACE"]` (L834),
  `["MASS_EVOLUTION_ROW"] = evolutionTerms["ACCUMULATION"] + evolutionTerms["ADVECTIVE"] + projectedFaceFlux[…]` (L836).
  The `origins` (L849-862: KINETIC/BULK_ENERGY/FACE/FLUX/ADVECTIVE/ACCUMULATION) PARTITION these rows.
  ⇒ the row object is DIRECTLY EMITTED by both engines (WL `operator[row]`; PY `EXPANDED`+advective+faces),
  so the row-assembly is computable AND `BUCKET_PARTITION_CHECK` (origins-sum vs operator-row) is a genuine
  (non-tautological) check. WL's `Inactive[Div]` only appears in the SEPARATE `DIVERGENCE_FORM_SOURCE`
  presentation (L274-280), not in the active `operator`/`origins` used for the row-sum. Advective uses
  `divergence`=`Sum D` (active, L247), ADVECTIVE=`divergence[sigmaEuler·∂_t u]`=full ∇·(σ_E u_t) (L694).
- KINETIC mechanism CONFIRMED by direct read: WL `kineticEw = muW WZero^2 D[eWField,{time,2}]` (L841) —
  the thickness inertia coefficient FROZEN at W_0² — vs PY `mu_W·W_bg²·e_tt`. §2a: "Every explicit `W₀`
  factor in `U` and in the bindings … is the physical background thickness and becomes `W_bg(y)`"; §3b: "do
  not freeze a coefficient at its constant binding before differentiation." ⇒ PY-correct is the LEADING read
  (WL under-retains W_bg²→W_0²), BUT confirm the §2a "explicit W_0 factor in U" scope covers the KINETIC
  inertia coefficient (kinetic energy, not stored energy U) — a §2a/§3b reading the instrument's legs / a
  spec check settle (one-engine fix is a spec question first). The instrument computes RES_row =
  μ_W·e_W_tt·(W_bg²−W_0²); the WHICH-ENGINE verdict is the §2a/§3b reading.
- ADVECTIVE re-bucketing lives WITHIN `MASS_EVOLUTION_ROW`: PY ∇ρ_br·u_t (advective_constraint) vs WL full
  ∇·(σ_E u_t) split across ACCUMULATION+ADVECTIVE; the row-sum over {ACCUMULATION, ADVECTIVE, FLUX}
  reconciles it iff representational (RES_row=0). The fingerprint's `ADVECTIVE` family = this mass-row term.
