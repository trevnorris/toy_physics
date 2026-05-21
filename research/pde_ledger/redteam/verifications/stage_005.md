---
unit_id: 005
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 005

## Per-finding outcomes

### F1 — missing_verification_script (subtype `missing_mathematica`)

**Classification:** resolved

**What changed:**
A new Mathematica audit script was created at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.wl` (lines 1-160). The diff (`exec_logs/stage_005_diff.patch`) shows this is the only file added; no other files were touched. The script defines `Wg`, `Wgp`, `projW`, `projWPrime`, `boundaryValue` from first principles using `Integrate[..., {w, -Infinity, Infinity}, Assumptions -> ...]`, `D[...]`, and `Limit[...]`, then asserts each of M1–M5 with `If[Simplify[lhs - rhs] =!= 0, Print["FAIL: ..."]; Exit[1], Print["PASS: ..."]]`, and emits `STATUS: PASS` as the final line.

**Assessment:**
The edit is substantive, not a transliteration of the SymPy choreography. Concrete checks:

- **Mathematica-native primitives.** Projections call `Integrate[..., {w, -Infinity, Infinity}, Assumptions -> projectionAssumptions]` directly (lines 20-28), wrapped in `FullSimplify`. Boundary values use `Limit[expr, w -> +/- Infinity]` (lines 30-33) — these are Mathematica idioms; SymPy does not have a direct analog and the SymPy script computes boundary by polynomial substitution at `w = +/- oo`. No Python helper names are mirrored.

- **Test profiles differ from SymPy** (per directive § F1 step 3):
  - M1: `Phi = (Sin[t] + x^2)(w^2 + 3)` vs. SymPy `(t^2 + x)(w^2 + 1)` — different t- and x-dependence and different polynomial weight in w.
  - M2: `Q = w Exp[-w^2/4]` vs. SymPy `w^3 + w` — Mathematica's profile has nontrivial Gaussian decay (so `Integrate` must evaluate a real Gaussian moment) rather than a polynomial that vanishes only because of the kernel.
  - M3: `G0 = Cos[t](w^2+2)`, `Gx = x^2 w^2`, `Gy = y(w^4+1)`, `Gz = z w^2`, `Gw = w + w^3/3`, `Gamma = (Sin[t] + x)(w^2 + 1)` — substantially different from the SymPy polynomial test fields.
  - M4: `Gw = w` — intentionally identical, the documented cross-check point with SymPy line 278.
  - M5: same `J^a` set as M3 minus `Gamma` — also differs from SymPy.

- **Non-tautological assertions.** Each residual is `Simplify[lhs - rhs]` where lhs and rhs are computed by independent `Integrate` / `D` calls, not by re-substituting the same intermediate. The mathematica log shows real, nonzero mutant residuals on M1–M5 (e.g. M2 mutant = `16/(5 Sqrt[5])`, M3 mutant = `3`, M4 mutant = `-2`, M5 mutant = `3`), confirming the algebra exercises live integrations rather than collapsing to `0 == 0`.

- **Boundary discharge verified, not assumed.** Lines 55-60 (M2), 83-88 (M3), 136-141 (M5) compute `boundaryValue[Wg[w] * profile]` and assert it equals zero before invoking the boundary-decay IBP form. Mathematica log confirms all three boundary values are `0`.

- **Mutant guards present for each item.** M2/M3/M4/M5 each have a sign-flip mutant that produces a nonzero residual (logged). M4 additionally has a "nonzero leakage" guard (line 124). M1's mutant is `lhsM1 + rhsM1` — weaker than a sign-flip on the commutator because both sides are equal and nonzero, so the sum is twice the value; this catches the trivial-zero degeneracy and is the pattern the directive prescribed. Not a regression.

- **Cross-check point.** M4 reports `leakage = 1` (mathematica log line 18), matching SymPy line 278's `-Proj_{W'}[w] = 1`. Sign and value agree across engines.

- **Output freshness.** Script mtime 2026-05-21 01:13; output mtime 2026-05-21 11:50; mathematica log timestamped 2026-05-21T11:12 with `exit_code: 0`. Output is newer than the script.

No collateral edits in the diff. No other files were modified. SymPy script was untouched and no new SymPy exec was needed.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy exec log was generated for unit 005 in this iteration (only the new `.wl` required execution). The SymPy script was not modified, so this is expected.

**Mathematica:** exit=0. Notable lines from `exec_logs/stage_005_mathematica.log`:
- `M4 leakage cross-check = 1` (matches SymPy line 278, confirms engine agreement on the universal Gaussian moment).
- `M2 mutant residual = 16/(5*Sqrt[5])` and `M3 mutant residual = 3` (mutant guards report real, nonzero residuals — non-tautological).
- `M3 residual = 0` and `M5 residual = 0` (projected inhomogeneous Maxwell law and projected charge continuity both hold under Mathematica's independent integration).
- Final line `STATUS: PASS`.

**Output freshness:** Confirmed. `mathematica/output/moving_throat_pde_stage005_projected_maxwell_covariant_mathematica_audit.txt` mtime is 2026-05-21 11:50, newer than the script's 2026-05-21 01:13.

## Material-change assessment

`material_change`: false.

No derived results changed. The single edit adds an independent verification engine; it does not modify any identity, sign, or numerical value that downstream units consume. Downstream units that depend on the projected inhomogeneous Maxwell law or projected charge continuity already saw the SymPy form; the new Mathematica script confirms the same identities (residuals and leakage value) on a different test set. No re-audit of downstream units is required on substantive grounds.

## Side observations (non-blocking)

- The M1 mutant `lhsM1 + rhsM1` is a degeneracy guard (catches the trivial `0 == 0` case) rather than a sign-flip guard on the commutator. The other M2–M5 mutants are genuine sign-flip guards. This matches the directive's wording ("at least one mutant guard per assertion") and is not a finding, just an observation that M1's guard is weaker than the others.
- The script uses `FullSimplify` inside `projW`/`projWPrime`. For the polynomial-in-w test profiles in M3 and M5, this is unnecessary (the integrals close in elementary form) but harmless. Not a concern.
- M3 and M5 share the same `J^a = G^a` profile shape (`Cos[t](w^2+2), x^2 w^2, y(w^4+1), z w^2, w + w^3/3`). This is a stylistic overlap inside the same script and is acceptable — the directive listed them separately and both pass on their own assertion.

## Verdict justification

The single finding F1 (`missing_mathematica`) is resolved. The new `.wl` file exists at the prescribed path, exits 0, emits `STATUS: PASS` plus `PASS:` lines for each of M1–M5 and their mutant guards, uses Mathematica-native `Integrate`/`D`/`Limit` primitives rather than transliterated SymPy helpers, employs test profiles that genuinely differ from the SymPy script for M1/M2/M3/M5 (M4 intentionally overlaps as the cross-check point), reports the same Gaussian leakage value `1` as SymPy line 278, and produces real nonzero mutant residuals demonstrating the assertions are non-tautological. No regressions are visible in the diff; the SymPy script was not touched.
