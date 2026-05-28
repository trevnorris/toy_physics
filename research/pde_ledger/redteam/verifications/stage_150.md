---
unit_id: 150
batch: IV.5
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 150

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage150_full_profile_residual_sympy_audit.py:37` — replaced
  `Sq = sp.simplify(sp.diff(Tq, x).subs(x, 0))` with two comment lines plus
  `Sq = Aq*k - Cq*Pi` (the hand-derived closed form for `T_q'(0)`). Now appears at line 39 with comments at lines 37–38.
- `mathematica/moving_throat_pde_stage150_full_profile_residual_mathematica_audit.wl:37` — replaced
  `sQ = FullSimplify[D[tq, x] /. x -> 0, Assumptions -> $Assumptions];` with a two-line
  comment plus `sQ = aq*k - cq*p;` (the hand-derived analogue). Now appears at line 39 with comments at lines 37–38.

Both edits exactly match the directive's "After" blocks; no collateral edits anywhere else in the diff.

**Assessment:**
The edit is correct and addresses the finding. The previous form `X - simplify(X) == 0` was a definitional tautology; the new form compares the symbolic derivative
`diff(Tq,x).subs(x,0)` against an independently-written closed-form combination of the constants `Aq`, `Cq`, `Pi`. Any transcription error in `Sq` (e.g. swapping `Aq*Pi - Cq*k`) would now fail the assertion. The exec logs confirm both engines still simplify `T_q'(0) - S_q` to 0 (sympy line 19, mathematica lines 12–13 `PASS: T_q'(0)-S_q`), as expected since the derivative of the closed-form `Tq` at `x=0` does equal `Aq*k - Cq*Pi`. The load-bearing `R''(0) - target = 0` check is unaffected and passes in both engines (sympy line 28, mathematica lines 19–20 `PASS: R''(0) - target`).

Note: Codex did not append the `## Applied: F1` block to the directive file as requested by the directive header, but the actual file edits are present and correct. Flagged as a procedural side observation, not a verification failure.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `S_q(Pi) =` block (lines 5–15) now displays the compact closed form `-Π^2/((1-e^(-Π))(-Π^2 + π^2/4)) + Π·(Π·e^(-Π) + π·sinh(π/2)/2)/((1-e^(-Π))(-Π^2 + π^2/4)·cosh(π/2))` rather than the prior expanded derivative.
- `T_q'(0)-S_q = 0` (line 19) — non-tautological check passes.
- `R''(0) - target = 0` (line 28) — load-bearing assertion passes.

**Mathematica:** exit=0. Notable lines:
- `S_q(Pi) = -(p^2/((1 - E^(-p))*(-p^2 + Pi^2/4))) + (p*Sech[Pi/2]*(p/E^p + (Pi*Sinh[Pi/2])/2))/((1 - E^(-p))*(-p^2 + Pi^2/4))` (line 5) — compact closed form mirrors the SymPy display.
- `T_q'(0)-S_q = 0` / `PASS: T_q'(0)-S_q` (lines 12–13).
- `R''(0) - target = 0` / `PASS: R''(0) - target` (lines 19–20).
- `Stage 150 Mathematica audit passed.` (line 26).

**Output freshness:** confirmed. SymPy `.py` mtime 2026-05-27 19:50:15; `.txt` mtime 2026-05-27 19:51:36 (re-generated after edit). Mathematica `.wl` mtime 2026-05-27 19:50:20; `.txt` mtime 2026-05-27 19:54:50 (re-generated after edit).

## Material-change assessment

`material_change`: false.

The fix is to a check-assertion construction only. No derived numerical/symbolic results that downstream stages would consume have changed — the residual `R(x)`, the constants `Aq`, `Cq`, and the curvature `R''(0)` are all algebraically identical before and after the edit. The change only affects whether the `T_q'(0) - S_q` check can fail on a transcription error; it does not alter any output that other stages depend on.

## Side observations (non-blocking)

- Codex did not append the `## Applied: F1` block to `redteam/directives/stage_150.md` as the directive's apply-protocol header requested. The file edits themselves are exactly as specified, so this is a paperwork omission rather than a functional defect. Orchestrator may want to flag this for codex iteration discipline but verification of the underlying finding stands.
- The auditor's documented Mathematica-transliteration observation (lines 119–140 of the report) was intentionally not raised as a directive finding by the auditor and remains an open structural note for future work, not a defect of this iteration.

## Verdict justification

The single finding (F1) is fully resolved: both engines now define `S_q` from the hand-derived closed form `Aq*k - Cq*Pi` (sympy) / `aq*k - cq*p` (mathematica), matching the directive verbatim. Both exec logs exit 0 with all assertions passing, including the (now non-tautological) `T_q'(0) - S_q = 0` check and the load-bearing `R''(0) - target = 0` check. Output files are fresh. No regressions in the diff. Verdict: **verified**.
