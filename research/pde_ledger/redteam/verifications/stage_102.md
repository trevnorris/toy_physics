---
unit_id: 102
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 102

## Per-finding outcomes

### F1 — missing_verification_script (script_doesnt_cover_claim, sympy side)

**Classification:** resolved

**What changed:**
The SymPy script `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage102_higher_odd_irrelevance_sympy_audit.py` was rewritten (orchestrator-direct application) to introduce a local `expect_zero(name, expr, tol=None)` helper (lines 20-32) that raises `AssertionError` if `sp.simplify(expr) != 0`, and three call sites (lines 51-54) that exercise the three deliverables:

- `expect_zero("D1: tauQ irrelevance at omega^5", tau5)` where `tau5 = sp.simplify(sp.diff(sp.im(Yser5.coeff(omega, 5)), tauQ))` (line 46) — matches M1.
- `expect_zero("D2: tauQ coefficient at omega^7 - 1/4", tau7 - sp.Rational(1, 4))` where `tau7 = sp.simplify(sp.diff(sp.im(Yser8.coeff(omega, 7)), tauQ))` (line 47) — matches M2.
- `expect_zero("D3: canonical odd coefficient at omega^5", sp.im(Yser5.coeff(omega, 5)) - chiQ * sp.Rational(9, 32) / Omega**5)` — matches M3.

The Mathematica file was also touched for the orchestrator's banner sweep: line 26 `banner["STAGE 102 — HIGHER ODD IRRELEVANCE"]` (was `"STAGE 085 — HIGHER ODD IRRELEVANCE"` per the auditor's cosmetic note). No math content in the `.wl` changed.

**Assessment:**
The three new asserts gate exactly the three claim-manifest items M1, M2, M3 from the directive. The expressions are non-tautological: each derives from the closed-form `Y = 3/4 + (1/4)/(1 - omega^2/Omega^2 - I*chiQ*sigma_can*omega^5 - I*tauQ*omega^7)` via `sp.series` and coefficient extraction, with comparison against the externally-stated rationals (`0`, `1/4`, `9/(32 Omega^5)`). The directive's prescribed structure used `assert` directly; the orchestrator instead routed through an `expect_zero` helper that performs the same `sp.simplify(...) == 0` check and raises `AssertionError` on failure — semantically equivalent gating, with the bonus of printing a `PASS: ...` line for transcript clarity. Sole deviation from the literal directive: helper-routing and the final banner reads `STAGE 102 AUDIT PASSED` instead of `Stage 102 SymPy audit passed.` — cosmetic only. Mathematica banner sweep on line 26 is a non-substantive label fix the auditor flagged as cosmetic.

## Exec log assessment

**SymPy:** exit=0 (transcript ends in `STAGE 102 AUDIT PASSED` and shell exit 0 is implied by clean termination of the script — no `AssertionError` traceback present). Notable lines:

- `tauQ coefficient in omega^5 term = 0` (line 3)
- `tauQ coefficient in omega^7 term = 1/4` (line 4)
- `PASS: D1: tauQ irrelevance at omega^5 (residual = 0)` (line 5)
- `PASS: D2: tauQ coefficient at omega^7 - 1/4 (residual = 0)` (line 6)
- `PASS: D3: canonical odd coefficient at omega^5 (residual = 0)` (line 7)

PASS line count: 3 (matches expected 3 for D1, D2, D3).

**Mathematica:** exit=0 (`Exit[0]` reached at line 53 of `.wl`; no `Exit[1]` from `fail[]`). Notable lines:

- `tauQ irrelevance at omega^5 = 0` then `PASS: tauQ irrelevance at omega^5` (lines 9-10)
- `tauQ coefficient at omega^7 - 1/4 = 0` then `PASS: tauQ coefficient at omega^7 - 1/4` (lines 11-12)
- `check canonical odd coefficient = 0` then `PASS: check canonical odd coefficient` (lines 13-14)

PASS line count: 3 (matches expected 3). Engines agree term-by-term modulo variable-name (`Omega` vs `omegaQ`) and ordering — same agreement as in the original audit.

**Output freshness:** confirmed.
- SymPy: script mtime 1779901959, output mtime 1779913725 (output newer).
- Mathematica: script mtime 1779901947, output mtime 1779913861 (output newer).

Both outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The edit only adds gating around already-computed (and previously printed) residuals that were known to equal `0`, `1/4`, `0`. No derived numeric result, symbol, or downstream quantity changes. The Mathematica banner-label sweep is cosmetic. No downstream unit can depend on the print/assert distinction.

## Side observations (non-blocking)

- The orchestrator-direct rewrite used a helper function (`expect_zero`) rather than raw `assert` lines. Functionally equivalent — both raise `AssertionError` on failure, which is what gating requires — but the directive prescribed literal `assert ... == 0` statements. Not a defect; flagging only because future verifier readers comparing directive text to file diff may notice the divergence.
- The auditor's other cosmetic notes (paper card title text reading "Stage~119" while `\label{stage:102}` is unit 102) are outside the verifier's scripts-only scope and were not addressed in this directive; that is consistent with the auditor not raising them as findings.

## Verdict justification

The single finding F1 is fully resolved: three `expect_zero` calls now gate D1, D2, D3 in the SymPy script with non-tautological residual expressions; the post-fix transcript shows three `PASS:` lines plus the new closing banner; both engine exit codes are 0; outputs are fresh; the Mathematica side already gated these and continues to pass. No regression in either script or log. Verdict: `verified`.
