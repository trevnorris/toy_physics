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
