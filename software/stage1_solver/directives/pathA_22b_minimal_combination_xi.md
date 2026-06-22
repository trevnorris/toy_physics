# Directive pathA_22b — GATED derivation of the minimal combination ξ → the REAL GR-quadrupole verdict

**Status:** DRAFT v2 (REWRITTEN as a gated plan after Claude+Codex+GLM reconciliation) — for Codex design-review (read-only) →
fixes → confirm-pass → execute Gate 0. **Owner split:** Codex derives + codes + iterates; Claude reviews. Orchestrator owns this
directive + decisions.
**Gate context:** User-gated 2026-06-21 ("derive the minimal combination ξ holistically"). v1 (holistic single directive) was
design-reviewed `UNSOUND` (over-scoped + leaned on an unproven cancellation). GLM tertiary reframed it as a TWO-LAYER wall, one
layer breachable; a Codex SOURCE-VERIFICATION (`_scratch/codex_pathA_22b_glm_digest_verify.log`) confirmed the reframe only
PARTIALLY — so this v2 turns each optimistic claim into a FAIL-ABLE gate with a real proof obligation. Resume:
`decisions/13_emergent_constants_derivation.md` §0 (5g)–(5j) + §1–§5. Predecessor: `reports/pathA_22a_dimensional_skeleton.md`.

## 0. The target + the reconciled facts (what is VERIFIED vs the sources — do not re-assume)
Target (derive each factor TARGET-BLIND, assemble, compare LAST):
```
ξ·χ_Q·P0  =  P0 · χ_Q · g_mhat² · λγ⁵ / g_G     ?=?   54/5
```
`G=(a·c_s²/m_GNLS)·g_G`, `m̂0=(c_s/(a²√m_GNLS))·g_mhat`, `c=λγ·c_s`, `S_port ≡ χ_Q`. VERIFIED facts (Codex source-check; cite when
used):
- `P0=N0/D0` is a TARGET-BLIND INPUT from the EXISTING converged finite-core branch (do NOT re-solve a deep/empty throat); use the
  pathA_22a `(c_s/a)²`-normalized `P0`.
- ω⁵ outgoing coefficient depends ONLY on `P0` and `χ_Q` (`pde.tex:2034-2069`); canonical `σ_Q^can = 9/(8Ω_Q⁵) = 4a⁵/(27c_s⁵)`
  (`pde.tex:1985-1988`), odd fingerprint coeff `a⁵/(27c_s⁵)` (`pde.tex:1964,2053`). [GLM mis-cited `9a⁵/(8Ω_Q⁵)` — use the real value.]
- `Z(w)` = localized Maxwell weight (`pde.tex:357-416`) is DISTINCT from `W(w)` = brane projection kernel (`pde.tex:277,496-505`).
  Do NOT conflate. `Z_int=∫Z(w)dw` defined `pde.tex:289-295`; `μ₀^eff=μ₀/Z_int`, `q_eff=q_*/√Z_int` under zero-mode assumptions
  `pde.tex:541-563` (a COUPLING normalization, NOT a photon-speed law). `c_bulk²=C_B/C_E` `pde.tex:541-565`; the bulk→brane SPEED
  normalization is UNSPECIFIED (pathA_20b `BULK_METRIC_SPEED_NORMALIZATION_UNSPECIFIED`). `Z_w` exported
  `m1c_background_export.py:166-170`.
- `m̂` INCONSISTENCY (real): dimensionful `m̂²Γ5=2G/(5c⁵)` (`pde.tex:2075-2083`) vs dimensionless `m̂=1+O(a²/r²)`
  (`pde.tex:2095-2099`); code uses the dimensionful reading (`dimensional_check.py:4280-4288`). Reconciling MAY flip the pathA_22a
  `TUNABILITY_CHANNEL_PRESENT` reading (if `m̂` is dimensionless it is a FIELD factor, not the scale carrier).
- Actual outgoing-DtN data are branch-dependent + unresolved (`pde.tex:2845-2849`); current code `χ_Q` is a target-derived
  DIAGNOSTIC, not a dynamic DtN extraction (`m1c_physical_run.py:460-467`; `S_port` defaults 1.0 `patha_extraction.py:535-544`).

**Outcome classes (per gate + overall; ALL must be reachable — no forced clean result):** `DERIVED` / `CANCELS` vs
`DOES_NOT_CANCEL` / `STILL_TUNABLE_<knob>` / `CONDITIONAL_ON_<x>` / `BLOCKED_NEEDS_<x>` / (overall) `REAL_MATCH` / `REAL_MISS` /
`FAIL_ABLE_PENDING_<residuals>`.

## 1. The gates (HARD STOPS — a gate that yields an irreducible free knob / blocked derivation STOPS the chain; do NOT fake-close)

**GATE 0 — Foundations (SymPy, no numerics; the cheapest, fail-fast). Two deliverables:**
- **0a — reconcile `m̂`.** Determine which `m̂` reading is correct for the Path-A normalization law, and whether it changes the
  pathA_22a decomposition (dimensionful = scale carrier `m̂0²` carries `M⁻¹L⁻²T⁻²`; dimensionless = `m̂` is a field/source-map factor
  and the missing scale lives elsewhere). This GATES the meaning of every later factor — do it FIRST. Outcome:
  `MHAT_DIMENSIONFUL_CONFIRMED` / `MHAT_DIMENSIONLESS_REFRAME_<where the scale moves>`.
- **0b — prove-or-fail the `Z`-kernel cancellation in `g_mhat²/g_G`.** This is a REAL PROOF OBLIGATION, NOT an assumption: show
  whether `g_G` and `g_mhat²` BOTH reduce to the SAME factorizable scalar transverse functional `I_Z=∫√g_w Z(w)dw` multiplying
  SEPARATE field-content kernels (⇒ `I_Z` cancels in the ratio, `W_eff` off the critical path for the ratio) — OR are ratios of
  WEIGHTED AVERAGES `∫Z·K_stress / ∫Z·K_source` where `Z` does NOT cancel (⇒ full `W_eff` needed; wall real). Derive from the
  parent action; state the exact structural condition and whether it holds. Outcome: `CANCELS` / `DOES_NOT_CANCEL`. **Gate 0b
  gates Gate 4.**

**GATE 1 — `Z_int` evaluation (one integral).** Evaluate `I_Z = Z_int = ∫Z(w)dw` from the EXISTING exported `Z_w`
(`m1c_background_export.py:166-170`); note finite-domain quadrature vs the ideal integral. Report it as a COUPLING-normalization
artifact ONLY (`μ₀^eff=μ₀/Z_int`). DO NOT promote it to `λγ`.

**GATE 2 — `λγ = c_γ/c_s` (speed reduction — still a derivation).** Derive the bulk→brane PHOTON-CONE SPEED map (NOT just the
coupling `Z_int`). pathA_20b left `c_bulk²=C_B/C_E` with the speed normalization UNSPECIFIED — close it or carry `C_B/C_E`
explicitly. Do NOT force `c_γ=c_s` (pathA_20b's negative control rejects it). Outcome: `LAMBDAGAMMA_DERIVED` /
`STILL_TUNABLE_LAMBDAGAMMA` / `CONDITIONAL_ON_<C_B/C_E or speed-map>`.

**GATE 3 — actual-branch `χ_Q` (linear outgoing-DtN BVP).** Linearize the Euler–Lagrange equations around the FROZEN existing
finite-core branch, impose outgoing BCs, solve the low-ω linear BVP, extract the ω⁵ coefficient, set
`χ_Q := (actual ω⁵ coeff)/(a⁵/(27c_s⁵))`. This is a NEW LINEAR-RESPONSE solve (NOT a nonlinear profile re-solve; does NOT
re-introduce the off-path deep solve). Outcome: `CHI_Q=1_CANONICAL` / `DELTA_Q≠0_<value>` / `NEEDS_DYNAMIC_SOLVE`.

**GATE 4 — residual field-content ratio `g_mhat²/g_G` (ONLY IF Gate 0b = `CANCELS`).** Derive the O(1) ratio of field-content
kernels (`K_stress` = derivatives of the field, the Noether stress; vs `K_source` = field value / outgoing projection) from the
parent action + the existing profile. If Gate 0b = `DOES_NOT_CANCEL`, this gate becomes the full `W_eff` derivation — flag
`BLOCKED_NEEDS_W_EFF` and STOP (do not fake a ratio). Outcome: `RATIO_DERIVED_<value>` / `BLOCKED_NEEDS_W_EFF`.

**GATE 5 — assemble + compare (arithmetic).** Assemble `ξ·χ_Q·P0 = P0·χ_Q·g_mhat²·λγ⁵/g_G`; dimensional-check the assembled
combination (must be dimensionless and match pathA_22a); THEN — and only then — compare to `54/5`. Outcome:
`REAL_MATCH`/`REAL_MISS` (all factors derived) or `FAIL_ABLE_PENDING_<the still-open gates>` (with the named residuals + the
minimal next derivation to close each).

## 2. Method + discipline (binding)
- **TARGET-BLIND:** derive every factor with NO reference to `54/5` or `10.8/P0`; the comparison is Gate 5 only. Reverse-engineering
  any factor from `10.8/P0` is FORBIDDEN ([[feedback-calibrate-predict-methodology]]).
- **DERIVED-FORM GATE:** no `x==x` posing as a check; no hand-inserted field/coefficient/sign; no `g_*` set to 1 unless DERIVED; no
  restatement faking closure (the pathA_21 P1 trap). An HONEST negative/conditional is a VALID outcome — record it, do NOT paper over.
- **Dimensional-check** (units restored, SymPy `M,L,T`) on every factor + the assembly ([[feedback-dimensional-consistency-check]]).
- **Dual-engine** (`.wl` via `math -script` — NOT `wolframscript`, which exits 255 in this env) wherever Mathematica can
  independently verify ([[feedback-dual-engine-required]]).
- `timeout 600` on every script (SymPy/`.wl` DERIVATION; NOT a solver run). Gate 3's linear BVP is a LINEAR solve; if it threatens
  to exceed 600s, it is mis-scoped — reformulate, do not raise the cap ([[feedback-script-timeout-policy]]).
- Additive: new module(s)/group(s); do NOT modify faithful operators or the pathA_18/19/20/20b/21/22a groups. Sources under
  `research/pde/paper/pde.tex` (NOT `pde_ledger/paper`), `research/pde_ledger/notes/stages/...` (Stage 104/105),
  `src/stage1_solver/{patha_extraction,m1c_background_export,m1c_physical_run}.py`, the pathA_20/20b/21/21b/21c reports.
- Do NOT commit; do NOT touch git or `decisions/*`. Leave outputs for review.

## 3. Execution granularity (gates are SEPARATE executions with user/orchestrator gating between)
Run Gate 0 FIRST as its own execution (cheapest, and 0a/0b can reframe everything). STOP after Gate 0; bring the result to the
orchestrator before Gates 1–5. Subsequent gates each: derive → dimensional-check → dual-engine → STOP at the gate's outcome.

## 4. Deliverables (per gate)
- A report section in `reports/pathA_22b_minimal_combination_xi.md` (append per gate) — the derivation, provenance line refs, the
  target-blind attestation, the dimensional check, and the gate outcome (§0 classes) + residual ledger.
- New SymPy module/group(s) + `.wl` cross-checks (`math -script`); tests incl. a target-blindness guard (the comparison value is
  not used upstream) and, for Gate 0b, a negative control (a non-factorizing kernel pair must yield `DOES_NOT_CANCEL`).
- Run logs; do NOT commit.

## 5. Review plan (after each gate's Codex exit 0)
Per-factor dimensional-fidelity audit ([[feedback-transliteration-fidelity-audit]]: derived factor vs source, term by term) +
adversarial review (restatement? reverse-engineered? `g_*`-set-to-1? verdict forced? is Gate 0b a real proof, not an assumed
cancellation?) — separate clean agents. For the FINAL assembled combination (Gate 5), route through the GLM tertiary consult
([[feedback-decisive-test-not-tautological]]). Claude synthesizes; bring each gate outcome + (at Gate 5) the GATE-A
freeze-amendment implication (installing any derived value un-pins `m̂0²·S_port=1`, hash `ed3585…`, decision-11 §2/§3 — a user
methodology call) to the user.
