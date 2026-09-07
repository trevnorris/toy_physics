# Decision review (ROUND 2) — the blind Wolfram N6 build directive, after fold

## Artifact
`research/pde_ledger_v3/directives/S11c_c2_N6_wl_build_directive.md` (orchestrator-written; physics-bearing build
directive for a NEW blind Mathematica engine). This is **round 2**: round 1 (two decision legs) returned FOLD-REQUIRED
with 9 findings each; all were verified and folded. Review the CURRENT directive **until clear**. Two obligations:
(a) do your own fresh pass against the sources (a fold can breed a new defect — review with full rigor, ⛔ do not assume
the fold is correct); (b) confirm each round-1 fold actually landed and is faithful. This is a DECISION review of a
directive — the `.wl` engine does not exist yet, so ⛔ no fictional-script ablation; executable tests are deferred to the
build. ⛔ You are reviewing the DIRECTIVE, not re-adjudicating the settled covariance verdict (Reading B).

## The round-1 folds to confirm (recorded in `_measurements/S11c_c2_N6_wl_directive_review_adjudication.md`)
Confirm each landed correctly and introduced no regression:
1. **Answer leakage removed** — `RESOLVED.md` is no longer a handed authority (now named only as "absent from the
   builder's context"); no inline computed outcome remains (`C_E=C_M`, `V_E≡V_M`, `R_N6=0`-case, `R_cov` value, SymPy
   `δ`). The junk knife is located by `a_ρ=h_α=0`, not by a residual value.
2. **Two WL namespaces** `WL_S11CC2_N6RC_*` / `WL_S11CC2_N6COV_*` matching the SymPy object names for the T7 join.
3. **Kernel families named by phase/bound-vars/measure** (not "6/9/12"); `I`/`B` = the S11c-b §3c weak restriction
   (six blocks) of the pressure-slot increment; emit key = (weak-block, kernel-family, grade, face).
4. **`μ_M` stated** as the route-2 energy pullback + EL, distinct and separately parameterizable from `μ_E∘Φ`.
5. **`S11c_b_SHARED_PHYSICS`** added as a sibling spec (energy §3a / operator §3b / weak restriction §3c).
6. **Source amplitude** = the μ/V-slot coefficients of the re-derived c1 source (not a recoded parallel formula).
7. **PIT** = exact finite fields only, `k_out≠k_in`, branch-cell coverage, WL-derived `δ`.
8. **Slot guards** `G_lin`/`G_cross`/`CLOSURE_EQUIVALENCE` emitted.
9. **Provenance/independence emitted** (`FROZEN_RELATIONS`, `PROVENANCE`, baseline/control-delta, `SPLIT_SUM`).
10. **`DIMENSIONS` able-to-fail** (extra-`W_0` incompatible summand).
11. **Cross-anchoring** `S11CC2_ANCHORING_L_MINUS_M` scoped OUT (full c2 engine's object; build not claimed complete for
    all §5c).
12. **Control header** no longer forbids the deliberate Φ-coefficient knife.

## Source-of-truth (read to judge; quote both sides with `file:line`)
- `research/pde_ledger_v3/directives/S11c_c2_SHARED_PHYSICS.md` §5c, §§1–2, §6, §7.
- `research/pde_ledger_v3/directives/S11c_c2_N6_reconcile_directive.md`, `.../S11c_c2_N6_covariance_directive.md` (cleared).
- `research/pde_ledger_v3/scripts/S11c_c2_N6_reconcile_sympy.py`, `.../S11c_c2_N6_covariance_sympy.py` (cleared instruments — the objects the WL engine reproduces blind).
- `research/pde_ledger_v3/_measurements/S11c_c2_N6_route2_spec_astra.md` (material carrier + PIT + slot guards).
- `research/pde_ledger_v3/directives/S11c_b_SHARED_PHYSICS.md` (energy/operator/§3c restriction), `S11c_a_SHARED_PHYSICS.md`, `S11c_c1_SHARED_PHYSICS.md`, `S11b_SHARED_PHYSICS.md`.
- `research/pde_ledger_v3/directives/S11c_c1_wl_build_directive.md` (blind-WL precedent).

## Required method
For each finding, quote the directive AND the contradicting/omitted source with `file:line`; a physics claim without a
source `file:line` is discarded. Report a finding only if it catches a way the WL engine, built to this directive, could
compute the wrong thing or let a wrong claim be made. Cover: faithful engine-neutral translation of both checks; anything
still missing for the carrier bridge + `R_cov`; blindness by absence (no result-bearing measurement in the builder's
context; "re-derived" not "imported"; no designed-to-agree); rule-17 freezes; controls able-to-fail and one-sided; no
leaked value / residual-zero exit / VERDICT; the three caveats left open. ⛔ Do not propose making the WL engine agree
with SymPy.

## Output
Findings (directive `file:line`, source `file:line`, why it changes what is computed/claimed, minimal fix); note any
round-1 fold that did NOT land or that bred a new defect. Sound sections briefly. End with: CLEAR-TO-BUILD, or
FOLD-REQUIRED with the blocking items.
