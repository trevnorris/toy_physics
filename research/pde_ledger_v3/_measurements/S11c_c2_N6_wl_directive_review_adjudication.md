# S11c-c2 blind Wolfram N6 build directive — decision-review adjudication (2026-09-07)

Artifact: `directives/S11c_c2_N6_wl_build_directive.md` (orchestrator-written, physics-bearing build directive for the
NEW blind Mathematica N6 engine). Reviewed until clear (physics-bearing → G2/G4, ⛔ not one-pass). 2 decision legs
(orchestrator-written → **Codex `gpt-5.6-sol` xhigh + Grok `grok-4.6` high**), identical prompt
`directives/_legs/S11c_c2_N6_wl_directive_review_prompt.md`. Logs (ephemeral /tmp):
`scratchpad/{codex,grok}_wl_n6_decision.log`.

## Round 1 — both legs EXIT=0, both FOLD-REQUIRED (9 findings each; ~12 distinct, all verified)
The gate paid off: the two legs were **complementary** (one leg alone would have shipped a directive missing half the
defects). I verified every finding against the sources (directive, cleared reconcile/covariance directives, route-2
spec, §5c, the SymPy N6 instruments, the c1 WL directive) and grounded three fold-decisions in the code
(`ANCHORING_L_MINUS_M` home; the slot-guard names; the PIT contract). All verified TRUE. Folds applied one pass →
re-review (round 2).

| # | Issue (verified) | Legs | Fold applied |
|---|---|---|---|
| 1 | **Answer leakage** — `RESOLVED.md` handed as a physics authority (it states the computed residuals + verdict); inline leaks (`C_E=C_M`/"representation-independent"/"must reduce to zero"; `V_E≡V_M`; `ΔC=0`; `R_N6=0` case; `R_cov` no-nonzero; SymPy `δ≈2.6e-22`) | Grok F1,F6 · Codex 2,3 | Dropped `RESOLVED.md` from the authority list; blindness clause now states result-bearing measurements are **absent from the builder's context**; scrubbed every inline outcome (carrier→"testing whether…"; compute `V_E`,`V_M` both; junk-knife located by `a_ρ=h_α=0`, not `R_N6=0`; caveats reworded; `δ` now WL-derived) |
| 2 | **WL tag namespace can't join** — used flat `WL_S11CC2_N6_*`; SymPy objects are `S11CC2_N6RC_*` / `S11CC2_N6COV_*` | Codex 1 | Split into two namespaces `WL_S11CC2_N6RC_*` / `WL_S11CC2_N6COV_*` matching the SymPy object names (interface identifiers, not representation copying); updated the emit scaffolding |
| 3 | **"signatures 6/9/12" untranslated** (SymPy diagnostic key a blind WL builder can't act on) + **missing §3c weak extract / six blocks** in `I`/`B` | Grok F2,F3 · Codex 4 | Replaced integers with the named formal-integral kernel families (non-integral remainder / one-momentum / two-momentum / second-order-middle, by phase/bound-vars/measure); stated `I`,`B` are the S11c-b §3c weak restriction (six blocks) of the pressure-slot increment; emit key = (weak-block, kernel-family, grade, face) |
| 4 | **`μ_M` pullback not translated** — engine-neutral section defined only Φ (prediction); conflating `μ_M` with `μ_E∘Φ` → circular `R_cov` | Grok F4 | Added the route-2 energy-pullback+EL construction of `μ_M` (Jacobian only in the quadratic scalar, degree-2 projection); kept `μ_M` and `μ_E∘Φ` separately parameterizable |
| 5 | **`S11c_b_SHARED_PHYSICS` missing** from the sibling-spec list (`μ_E`/`μ_M`/`C_E` come from S11c-b §3a/§3b/§3c, not §1c) | Grok F5 | Added it as a named sibling spec (energy §3a, operator §3b, weak restriction §3c, live jets); noted §1c is the already-folded θ-row, ⛔ not the energy |
| 6 | **Source amplitude recoded** — presented `b_{r,s}` as the object; it is the slot-IDENTIFICATION of the c1 source (recoding double-counts ε / drops slot factors) | Grok F7 | `es`/`ms` = the μ/V-coefficients of the re-derived c1 source (`Z`+resolvent stripped); the formula identifies the slots, ⛔ not a second expression |
| 7 | **PIT too weak** — "finite-field / numeric" licensed floats; `k_out=k_in` collapses off-diagonal; transferred the SymPy `δ` | Grok F9 · Codex 7 | Exact finite fields only; several primes; `k_out≠k_in` independently sampled; shared coords/inverses across routes; joint singular rejection; real-branch-cell coverage; **WL-derived** `D`/`E`/`δ`; bad-prime handling |
| 8 | **Slot-linearity assumed, not guarded** — asserted `S_P=ΣC·p` with no guard | Codex 5 | Emit `SLOT_GUARD` (`G_lin`), `SLOT_GUARD_CROSS` (`G_cross`), `CLOSURE_EQUIVALENCE` (`P→χ(P)`), denominator pressure-dependency check — measurements, not assertions |
| 9 | **Independence/provenance not emitted** — relegated to the prose report | Codex 6 | Added emitted `FROZEN_RELATIONS`, `PROVENANCE` (fingerprints + independent-builder census), `SOURCE/R_COV baseline`+`control-delta`, `ACTUAL_CONTROL_PARAMETERS`, `SPLIT_SUM` |
| 10 | **Dimensions labelled, not able-to-fail** | Codex 8 | `DIMENSIONS` now requires independently computed units + a deliberate extra-`W_0` incompatible summand (able-to-fail); unknown units reported, not coerced |
| 11 | **Cross-anchoring contract** `S11CC2_ANCHORING_L_MINUS_M` acknowledged but not emitted | Codex 9 | Grounded: it is emitted by the FULL c2 self-energy engine (`selfenergy_fold_sympy_audit.py:1090`), ⛔ not the N6 instruments → scoped OUT (full-c2 WL engine's object); this build not claimed complete for all §5c |
| 12 | **Control header over-forbade coefficient knives** — the Φ-coefficient knife (`2·a_ρ`) IS a deliberate one-sided coefficient change | Grok (sound-note) | Reworded the header: carrier/junk/source knives are FORM; the Φ-coefficient knife is the deliberate coefficient exception (what the truth table cannot see) |

**Sound (both legs, preserved):** the central algebraic translation (coefficient-level `C_E−C_M`; `es=(μ_E,V_E)`/
`ms=(μ_M,V_M)`; CARRIER uses `ms`; SOURCE/CROSS closed-response-only; rank-2 Φ prolongation; `R_cov` sign + non-circular
construction; `R_COV_INCREMENT` excludes the bare term); N6-only scope sufficient (route-2 permits the face/pressure-slot
fold); the four controls one-sided on the correct routes; no rule-17 freeze; the self-discrepancy-as-covariance
prohibition present; the three caveats accurately left open.

## Round 2 — both legs EXIT=0, both FOLD-REQUIRED; convergence tightening (12 → 4). Prompt `_r2.md`, logs `scratchpad/{codex,grok}_wl_n6_decision_r2.log`
Folds 3,4,5,6,7,10,11,12 confirmed LANDED by both legs. Four remaining (all verified against sources; 2 are regressions
my round-1 folds bred — precisely what review-until-clear exists to catch):

| # | Issue (verified) | Legs | Fold applied |
|---|---|---|---|
| R2-1 | **Fold-1 incomplete** — the two SymPy N6 **build directives** (reconcile, covariance) stayed builder-facing "winning authorities" while STATING the outcomes (`R_N6` ~18 cols, `C_E=C_M`, `V_E≡V_M`); a file cannot be both an authority and "absent from context" | Grok F1 · Codex 1 | Reclassified both (+`RESOLVED.md`) as **reviewer-side sources only, ⛔ not handed to the build**; builder-facing authorities = `SHARED_PHYSICS` §5c + route-2 spec + sibling specs + this translation; removed their bullets; scrubbed the outcome-naming from the rationale line itself |
| R2-2 | **Regression (from fold-1)** — removing `V_E≡V_M` left "the square commutes **iff** Φ is θ-independent…", now both answer-bearing (implies `R_cov=0`) and wrong (ignores the velocity channel the source square carries) | Codex 2 | Deleted the `iff`; `R_cov` now stated to measure the **complete combined-source square** incl. independently-computed `V_E`,`V_M`; μ-isolation is downstream after the velocity bridge is observed |
| R2-3 | **Regression (from folds 2+8)** — slot guards named `SLOT_GUARD`/`_CROSS`/`CLOSURE_EQUIVALENCE` under `N6RC_` don't match the SymPy **diagnostic** objects `S11CC2_N6_SLOT_GUARD_{NATIVE,CARRIER,RESIDUAL}`/`CLOSURE_GUARD_*` → can't T7-join | Codex 3 | Renamed to the diagnostic namespace `WL_S11CC2_N6_SLOT_GUARD_{NATIVE,CARRIER,RESIDUAL}` + `CLOSURE_GUARD_*`; `G_cross`/denominator data → slot-guard residual payload |
| R2-4 | **Fold-9 incomplete** — `FROZEN_RELATIONS` omitted the four **pressure-symbol identity/assumption census** (route-2:97; reconcile_sympy:60-73); distinct `δp` symbols could be minted → carrier silently zeroes | Codex 4 | Added the four pressure-symbol identity/assumption census (same inherited symbols incl. assumptions, ⛔ never minted by spelling) to `FROZEN_RELATIONS` + provenance |

Everything else round-2-confirmed sound (the central algebraic translation, both checks, the six-block restriction, Φ
prolongation, non-circularity, PIT, dimensions, caveats, scope). No leg contradicted another (Codex's F2–F4 are
complementary refinements Grok did not reach; no conflicting "sound" claim). ⚠ Rule-15 watch: 2 of the 4 were regressions
from my folds — if round 3 breeds further regressions from the round-2 folds, change the fold author (astra/Codex-sol).

## Round 3 — CLEAR-TO-BUILD. Prompt `_r3.md`, logs `scratchpad/{codex,grok}_wl_n6_decision_r3.log`
Both legs did a fresh full pass and confirmed all four round-2 folds landed with no recurrence of the two regressions.
- **Grok r3: CLEAR-TO-BUILD** — 0 findings; every fold + the central translation confirmed sound.
- **Codex-sol r3: FOLD-REQUIRED — 1 item** (verified, adjudicated): fold R2-3 landed in the object contract but the
  emit-scaffolding description (`standardEmissionName`) + the section heading still said "two namespaces", omitting the
  diagnostic-guard family `WL_S11CC2_N6_SLOT_GUARD_*`/`CLOSURE_GUARD_*` — a builder wiring the scaffold could route the
  guards into `N6RC_` and break the T7 join. ⚠ Grok saw the same spot and filed it non-blocking ("the emit contract
  wins"); the legs' substance did not conflict. G4: Codex right — a physics-neutral consistency fix (Codex: "No payload
  physics needs changing"). **Folded:** heading → "matched namespaces"; intro + `standardEmissionName` now wire all
  three prefixes (`N6RC_`/`N6COV_`/diagnostic `N6_`), routing the six guard objects through the diagnostic-family case.
  Mechanically verified (three-prefix wiring consistent, no "two namespaces" residue, final leaked-value lint empty).

**Stopping rule met (substantive, ⛔ not "both legs green"):** after the round-3 fold nothing outstanding changes what is
computed or may be claimed — the object contract was cleared by both legs, and the only round-3 change was a
physics-neutral scaffolding-consistency fix both legs' substance agreed on. Convergence 12 → 4 → 1(trivial). Rule-15
watch clear (round 3 bred no new regression). ⇒ **DIRECTIVE CLEAR-TO-BUILD.**

## Next
Commit the reviewed directive baseline → astra WL build (`gpt-6-astra` high; Mathematica; detached; 2-seat;
`--sandbox danger-full-access`) → verify deliverable → 2 build legs SERIALIZED (fresh Claude agent + Grok; kernel
ablation) → repairs/clearance → commit both `.out`.
