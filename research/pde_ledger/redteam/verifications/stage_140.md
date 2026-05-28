---
unit_id: 140
batch: IV.5
verifier_model: claude-opus-4-7
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 140

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:26` now reads `banner["STAGE 140 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE"];` (the orchestrator note confirms this was carried by the Cluster A mass-renumber pass).

**Assessment:**
The Mathematica transcript at line 3 of the output confirms `STAGE 140 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE` is printed. The `expectZero["Sigma_0_hat - 20 That^2/9", ...]` assertion still holds (transcript line 7-8: residual = 0, PASS). The directive's optional structural reordering of the algebra was not taken, which is acceptable (the directive marked it discretionary).

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- sympy `.py:25-28`: three new `assert sp.Abs(... - sp.Float(notes_value, 30)) < sp.Float('1e-12', 30)` checks for `That_nat`, `That_comp`, and the fractional enhancement, plus a `print('numeric fixed-point checks PASS')`.
- mathematica `.wl:55-63`: a new `Module[{diff1, diff2, diff3, tol}, ...]` block that issues three `If[Abs[diffN] < tol, pass[...], fail[...]]` checks against `SetPrecision[0.866512630228382, 30]`, `SetPrecision[0.901484054174206, 30]`, and `SetPrecision[0.0403588161624, 30]`.

**Assessment:**
The edits match the directive verbatim (file:line ranges, before/after blocks, tolerance `1e-12`). The new assertions are non-tautological: they compare floating-point quantities derived from `sqrt(9*Ms/20)` with `Ms_nat/Ms_comp` literals against the independently-recorded notes values; perturbing either the `Ms_nat`/`Ms_comp` input literals or the `sqrt(9*Ms/20)` formula would break the new asserts. The sympy transcript ends with `numeric fixed-point checks PASS` (line 6), and the Mathematica transcript prints `PASS: That_nat matches notes`, `PASS: That_comp matches notes`, `PASS: fractional enhancement matches notes` (lines 12-14). No collateral edits outside the named ranges.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Sigma_0 in terms of That = 20*That**2/9`
- `That_nat = 0.866512630228381532619923771658`
- `numeric fixed-point checks PASS`

**Mathematica:** exit=0. Notable lines:
- `STAGE 140 — SELFMATCHED MOUTH SUSCEPTIBILITY CLOSURE` (banner correct)
- `PASS: Sigma_0_hat - 20 That^2/9`
- `PASS: That_nat matches notes`, `PASS: That_comp matches notes`, `PASS: fractional enhancement matches notes`
- `Stage 140 Mathematica audit passed.`

**Output freshness:** confirmed. Scripts mtime 19:49; sympy output mtime 19:51; mathematica output mtime 19:52. Both outputs are newer than their corresponding scripts.

## Material-change assessment

`material_change`: false.

The edits add cosmetic banner text and three new numeric checks against pre-recorded notes values; no derived result, closed-form identity, or numerical constant changed. Downstream units consuming stage 140 outputs (`Sigma_0 = (20/9) That^2`, the `That_nat`/`That_comp` numerics) see identical values to pre-fix.

## Side observations (non-blocking)

None. The `.py` script remains 28 lines and the `.wl` script remains 68 lines; both retain their original derivation choreography (the optional structural reordering in F1 was declined, consistent with the directive marking it discretionary).

## Verdict justification

Both findings are `resolved`. The Mathematica banner now correctly reads STAGE 140; the three new numeric assertions (sympy + mathematica) guard the Section-3 numerics with a 1e-12 tolerance and are non-tautological. Both engines exit 0 with their fresh outputs containing the expected PASS lines. No regressions or collateral edits.
