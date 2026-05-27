---
unit_id: 092
batch: IV.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 092

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl:33-77` was restructured to work in the dimensionless `(eps2, eps4)` variables from the outset. The script now (lines 40–48) defines `k0Full/k2Full/k4Full` from the conservative-carrier expansion, substitutes `kg2 -> eps2*kPole/omegaQ^2` and `kg4 -> eps4*kPole/omegaQ^4` directly, then derives `K_0` via `4*k2Eps^2/k4Eps` (line 52) — no `Solve[branch == 0, kg0]` anywhere. It verifies the closed form against `4 kPole (1+eps2)^2/(1+eps4)` (line 57), defines `cPole = kPole/k0FromBranch` (line 60), and asserts the static-geometry limit on both `c_pole = 1/4` (line 66) and `c_geom = 3/4` (line 70).

**Assessment:**
The new algebraic path is genuinely independent from the SymPy script — the SymPy file still uses `Solve[branch == 0, Kg0]` (sympy:50) whereas the rewritten Mathematica file derives `K_0 = 4 K_2^2/K_4` directly from the eps-substituted forms. The directive's five verification criteria all hold: (a) no `Solve[branch` substring (confirmed by absence in the file); (b) `PASS: K_0 closed form matches 4 K_pole (1+eps2)^2/(1+eps4)` (output line 9); (c) `PASS: c_pole - (1+eps4)/(4(1+eps2)^2)` (output line 12); (d) `PASS: static-geometry limit c_pole = 1/4` and `PASS: static-geometry limit c_geom = 3/4` (output lines 14, 17); (e) exit 0. The new `K_0` assertion is non-tautological: `k0FromBranch` is built from `4 k2Eps^2/k4Eps` (with `k2Eps`, `k4Eps` arising from the omega-series of the conservative carrier with eps-substitutions) and compared against an independently-formed candidate `4*kPole*(1+eps2)^2/(1+eps4)`.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
SymPy: three `expect_zero` calls appended at `scripts/moving_throat_pde_stage092_dynamic_geometry_obstruction_sympy_audit.py:80-82` after the `Dropped higher-order tail` print, before the FINAL LEDGER banner — exactly the snippet specified in the directive (eps^0 coefficient = 1/4, eps2 coefficient = -1/2, eps4 coefficient = 1/4). Mathematica: three matching `expectZero` calls at `moving_throat_pde_stage092_dynamic_geometry_obstruction_mathematica_audit.wl:79-81` using `Coefficient[Coefficient[smallSeries, eps2, i], eps4, j]`.

**Assessment:**
Assertions match the directive's prescribed text byte-for-byte. Non-tautological: `small_series` is constructed by series-expanding `cpole_expected = (1+eps4)/(4*(1+eps2)^2)` (sympy:72-73; smallSeries in .wl:73) and the asserted coefficients (1/4, -1/2, 1/4) are independently-stated rational targets, not derived from `linear_part`. A regression that altered the predicted series would produce a non-zero coefficient and fail the script. Both transcripts show all three asserts evaluating to zero (sympy output lines 20-22, Mathematica output lines 21-26 with three PASS lines).

## Exec log assessment

**SymPy:** exit=0 (clean run; no AssertionError; FINAL LEDGER banner reached). Notable lines:

- `c_pole - (1+eps4)/(4(1+eps2)^2) = 0` (line 15)
- `first-order eps^0 coefficient = 0` (line 20)
- `first-order eps2 coefficient = 0` (line 21)
- `first-order eps4 coefficient = 0` (line 22)
- `Dropped higher-order tail  = -eps_2*eps_4/2` (line 19) — second-order as expected.

**Mathematica:** exit=0 (script reached `Stage 092 Mathematica audit passed.` at output line 28). Seven PASS lines counted, matching the seven `expectZero` calls in the file (4 from F1: K_0 closed form, c_pole formula, static c_pole=1/4, static c_geom=3/4; 3 from F2: three first-order coefficients). Notable lines:

- `PASS: K_0 closed form matches 4 K_pole (1+eps2)^2/(1+eps4)` (line 9)
- `PASS: c_pole - (1+eps4)/(4(1+eps2)^2)` (line 12)
- `PASS: static-geometry limit c_pole = 1/4` (line 14)
- `PASS: static-geometry limit c_geom = 3/4` (line 17)
- `PASS: first-order eps^0 coefficient`, `PASS: first-order eps2 coefficient`, `PASS: first-order eps4 coefficient` (lines 22, 24, 26)

**Output freshness:** Confirmed. SymPy script mtime 1779902001 vs output mtime 1779913705 (output ~3.25 h newer). Mathematica script mtime 1779902021 vs output mtime 1779913748 (output ~3.26 h newer). Both transcripts regenerated post-edit.

## Material-change assessment

`material_change`: false.

The edits hardened verification surfaces (new assertions, restructured derivation path) but did not modify any derived result. The load-bearing identity `c_pole = (1+eps_4)/[4(1+eps_2)^2]` was already proven in the prior versions and remains the conclusion; F1's restructure derives the same `K_0` closed form via a different algebraic route and gets the same answer; F2's coefficient checks merely codify what was already printed. Downstream units that consume the dynamic-geometry obstruction split are unaffected.

## Side observations (non-blocking)

- The SymPy script retains the `Solve[sp.Eq(branch, 0), Kg0]` path (sympy:50); only the Mathematica engine was restructured per F1. This is consistent with the directive (F1 targeted only the `.wl`) — the two engines now reach the same `c_pole` formula by two genuinely different algebraic routes, satisfying the second-engine independence policy.
- Both engines independently report `Dropped higher-order tail = -eps_2*eps_4/2` (sympy:19; .wl output:20), confirming the second-order content agrees across engines.

## Verdict justification

Both findings are fully resolved by orchestrator-direct edits matching the directive's prescribed text. The Mathematica derivation no longer mirrors the SymPy pipeline — it works in `(eps2, eps4)` from the start and reaches `K_0 = 4 K_pole (1+eps_2)^2/(1+eps_4)` via `4 K_2^2/K_4`, an algebraically independent route. Three new first-order coefficient assertions on each engine close the F2 gap; the asserted targets (1/4, -1/2, 1/4) are independent of the linear-part construction, so they are non-tautological. Both exec logs exit 0; the expected 7 Mathematica PASS lines and 3 new SymPy zero-assert lines are present; outputs are freshly regenerated. No regressions, no new findings.
