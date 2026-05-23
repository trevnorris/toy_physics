---
unit_id: 057
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 057

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:107-112` — the previous `expect_zero("y_req identity", y_req_sq - ((Omega_Pe**2/zeta_req)*(kappa+pi^2/4) - kappa))` self-subtraction was replaced with a round-trip check `expect_zero("y_req defining equation", zeta_req - sp.simplify(Omega_Pe**2*(kappa + sp.pi**2/4)/(kappa + y_req_sq)))`. The definition of `y_req_sq` at line 93 and the surrounding `print` at line 101 are unchanged.
- `mathematica/moving_throat_pde_stage057_physical_parameter_map_mathematica_audit.wl:96-102` — the previous `expectZero["y_req identity", yReqSq - ((omegaPe^2/zetaReq)(kappa + Pi^2/4) - kappa)]` self-subtraction was replaced with `expectZero["y_req defining equation", zetaReq - FullSimplify[omegaPe^2 (kappa + Pi^2/4)/(kappa + yReqSq), Assumptions -> $Assumptions]]`. The definition of `yReqSq` at line 82 and the `Print` at line 86 are unchanged.

**Assessment:**
The edit matches the directive's required patch character-for-character in both engines. The check is now non-tautological: substituting `y² → y_req_sq` into `Ω²(κ+π²/4)/(κ+y²)` yields `(κ+y_req_sq) = (Ω²/ζ_req)(κ+π²/4)`, so the residual collapses to `ζ_req − ζ_req = 0` only because the closed form `y_req_sq = (Ω²/ζ_req)(κ+π²/4) − κ` is correctly derived — any sign error or factor slip in `y_req_sq` would leave a nonzero residual. No collateral edits to either script beyond the named block; the orchestrator's preemptive `ConditionalExpression -> e` strip in the Mathematica `expectZero` helper (lines 22-26) is a generic idiom fix and does not alter substantive logic. Saved outputs (sympy line 22, mathematica lines 37-38) confirm `y_req defining equation = 0` / `PASS: y_req defining equation` printing in place of the prior `y_req identity` label.

## Exec log assessment

**SymPy:** exit=0 (per `fix_batch_III.2.log` line 97; canonical transcript timestamp 2026-05-22 17:43:15, newer than script mtime 17:42:17). Notable lines:
- `A_K - (kappa+pi^2/4)/(kappa+y^2) = 0`
- `kappa_req identity = 0`
- `y_req defining equation = 0`
- `Stage 40 audit passed.`

**Mathematica:** exit=0 (per `fix_batch_III.2.log` line 98; canonical transcript timestamp 2026-05-22 17:43:21, newer than script mtime 17:42:18). Notable lines:
- `PASS: kappa_req defining equation`
- `PASS: kappa_req identity`
- `y_req defining equation = 0`
- `PASS: y_req defining equation`
- `Stage 057 Mathematica audit passed.`

The Mathematica `Limit::alimv` warnings on the `zeta_max - limit(Pe->inf, y->0)` line are nominal (auditor already accepted them — Mathematica drops the assumption on the limit variable; the residual still reduces to zero).

**Output freshness:** Confirmed. Script mtimes are 17:42:17 (sympy) and 17:42:18 (mathematica); transcript mtimes are 17:43:15 and 17:43:21, both strictly newer. Both transcripts contain the renamed assertion `y_req defining equation`, confirming they were regenerated against the patched scripts.

## Material-change assessment

`material_change`: false.

The patch swaps one assertion form for another that exercises the same closed form `y_req² = (Ω²/ζ_req)(κ+π²/4) − κ`. The printed values of `Omega_req^2`, `y_req^2`, and `kappa_req` are unchanged from the pre-fix transcript (definition lines untouched). No downstream unit can observe a derived-quantity change from this edit, only a stronger algebraic guarantee on the existing formula.

## Side observations (non-blocking)

- The sympy banner still reads `STAGE 40 — PHYSICAL (Pe, kappa, eta) PLACEMENT MAP` and the closing print is `Stage 40 audit passed.`, while the file is named `stage057_*`. The mathematica counterpart correctly says `STAGE 040` and `Stage 057 Mathematica audit passed.`. This is a pre-existing cosmetic mismatch unrelated to F1, flagged here for the auditor's attention on a future pass.
- The `expectZero` helper in the Mathematica script now contains a `ConditionalExpression[e_, _] :> e` strip after the first FullSimplify. This was applied by the orchestrator as a batch-wide idiom fix; for unit 057 specifically, none of the residuals exercise it (all reach pure `0`), so the strip is a no-op here. Not a finding.

## Verdict justification

The single finding (F1, tautological_check) was applied verbatim per directive in both engines, with no collateral edits to the surrounding definition lines. Both exec transcripts are post-fix (mtimes newer than scripts), exit 0, and contain the renamed `y_req defining equation = 0` / `PASS: y_req defining equation` assertion as the directive required. The new check is genuinely non-tautological — it round-trips `y_req_sq` through `Ω²(κ+π²/4)/(κ+y²)` and only vanishes because the closed form is correctly derived. No regressions visible in the transcripts; all prior passing assertions (A1–A7, B1–B9) continue to print `0`/`PASS`. Verified.
