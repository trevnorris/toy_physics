---
unit_id: 002
batch: I.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-21T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 002

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
The Mathematica audit `.wl` was edited in three places, all on `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.wl`:

1. Line 83: `y00 = 1/(2 Sqrt[Pi])` replaced with `y00 = FullSimplify[SphericalHarmonicY[0, 0, theta, phi]]`. The kernel-derived value of Y_00 now drives Section I, and `norm(Y00) - 1` becomes a non-trivial check on the built-in.
2. Lines 134-151: hand-built `Mintegrand`/`Kintegrand` was replaced by coefficient extraction (`MaaExt = 2 Coefficient[lw, dadt, 2]`, etc.) from the angular-integrated `lw`, plus a separately postulated `MintegrandBoxed`/`KintegrandBoxed`, with two new `expectZero` assertions: `extracted M matches boxed M (4 Pi overlap form)` and `extracted K matches boxed K (4 Pi overlap form)`. The previous single literal `lw - lwTarget` check was removed in this block. (`Mintegrand` / `Kintegrand` are then aliased to the extracted matrices so the downstream `Mmat`/`Kmat` integrals remain meaningful.)
3. Lines 3 and 189-191: `Needs["VariationalMethods\`"]` added at the top, and the local-helper EL calls were replaced with `{elAEq, elLEq} = EulerEquations[lredTime, {qa, qLfun}, t]` followed by `elA = elAEq[[2]] - elAEq[[1]]` (and the q_L mirror). The local helper `eulerLagrange1D` is preserved for the Section III single-component EL check at line 306, as the directive required.

**Assessment:**
The edits substantively replace SymPy-transliterated algebra with Mathematica-native primitives:

- `SphericalHarmonicY` is a real kernel function call, not a literal copy of the SymPy numeric.
- The coefficient extraction is genuinely independent of the boxed form: `lw` is built from the angular integral of `(muEta etaT^2 - Tw etaW^2 - K0 eta^2)*sin(theta)`, `Coefficient` peels off the quadratic forms in `(dadt, dLdt, da, dL)`, and the extracted 2x2 is then compared to a separately postulated 4 Pi-prefactor matrix. The comparison is non-trivial because the extraction route never references the boxed form. Codex deviated from the directive's literal `(1/(dadt^2))` and `(1/2)` factors (which are not the right scalars), using instead `2 Coefficient[..,2]` and `Coefficient[Coefficient[..],..]` — this is the mathematically correct extraction for a Lagrangian of the form (1/2)(M_AB Qdot^A Qdot^B - K_AB Q^A Q^B), since `Coefficient[(1/2) Maa dadt^2, dadt, 2] = (1/2) Maa` and the cross term `(1/2)(2 MaL) dadt dLdt = MaL dadt dLdt` gives `Coefficient[Coefficient[..,dadt],dLdt] = MaL`. Codex's `## Applied: F1` deviation note correctly documents this. The deviation makes the assertions correct rather than malformed; this is good-faith engineering, not a tautology.
- `EulerEquations` is a real `VariationalMethods` primitive. The sign convention is consistent: `EulerEquations` returns `dL/df - d/dt(dL/df') == 0`, so `elAEq[[2]] - elAEq[[1]] = -(dL/df - d/dt(dL/df')) = d/dt(dL/df') - dL/df`, which matches the original helper convention used by the downstream comparison `elA - (Maa qa'' + MaL qLfun'' + Kaa qa + KaL qLfun)`.

The `expectZero` name strings for the new assertions ("extracted M matches boxed M (4 Pi overlap form)") are distinct from any SymPy strings, breaking the byte-identical naming critique. The remaining shared physical symbols (muEta, Tw, K0, alphaA, alphaL) are the user-facing constants the auditor explicitly allowed.

The collateral change is the rewrite of `lwTarget` to use `MintegrandBoxed`/`KintegrandBoxed` instead of `Mintegrand`/`Kintegrand`. This is consistent with the new structure: the downstream `Lred - LredTarget` check now exercises the boxed-form matrices through the integral. Given that the extracted vs. boxed matrices are already separately asserted equal, this check is downstream-redundant rather than tautological — flagged in side observations, not blocking.

The Section III multi-component checks were added in F2 (next), and Mathematica's `SphericalHarmonicY` is used only in Section I (the directive scope); the Section III literals for Y_2m remain hand-written, which the directive did not call out.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
On the same `.wl` file:

- Line 213-214: `basis` extended to `{y20, y21c, y21s, y22c, y22s}` and `basisNames` likewise to all five. The existing `Do[..., {i, 1, Length[basis]}]` loops automatically iterate over five rather than three components.
- Lines 234-249: a new block builds `normMatrix5` (5x5 angular integral of `basis[[i]] basis[[j]] dOmega`) and `gradMatrix5` (5x5 angular integral of `gradS2Inner[basis[[i]], basis[[j]]] dOmega`) via `Table[...]`, then asserts they equal `IdentityMatrix[5]` and `6 IdentityMatrix[5]` respectively. The output `.txt` confirms both produce a literal 5x5 zero matrix after subtraction.
- Lines 265-266: `qvec`/`qdotvec` extended to all five `q_{lm}` symbols, so the per-component reduced-density loop on lines 268-297 iterates over all five.
- Lines 216-217: the supplementary phase-shift identities for Y21s and Y22s are retained, as the directive required.

**Assessment:**
The verification is now substantive across all five P2 components:
- `norm(Y21s) - 1` and `norm(Y22s) - 1` PASS via direct integration (output lines 53-54, 61-62).
- `angular energy(Y21s) - 6` and `angular energy(Y22s) - 6` PASS (lines 55-56, 63-64).
- `-Delta_S2 Y21s - 6 Y21s` and `-Delta_S2 Y22s - 6 Y22s` PASS (lines 75-76, 79-80).
- The 5x5 norm matrix and 5x5 angular-stiffness matrix produce literal 5x5 zero matrices when subtracted from `IdentityMatrix[5]` and `6 IdentityMatrix[5]` — these directly cover the off-diagonal entries (e.g. <Y21c, Y22s>, <Y20, Y21s>) that the audit specifically called out as not previously integrated.
- Single-component reduced density now PASSes for Y21s and Y22s as well (lines 85-86, 89-90).

The assertions are non-tautological because each basis element is defined by an explicit closed-form spherical harmonic (lines 208-212), the integrals are evaluated by `Integrate`, and the comparison is to a target (`IdentityMatrix[5]`, `6 IdentityMatrix[5]`) that is unrelated to the integration mechanism. No assertion was made true by reformulation.

The 3-of-5 weak verification claim that motivated the finding is now fully resolved.

## Exec log assessment

**SymPy:** exit=n/a. The orchestrator did not capture `stage_002_sympy.log` and the SymPy script was untouched by codex (the diff has only one `diff --git` entry, for the `.wl`). The SymPy output `.txt` mtime is 2026-05-21 11:25, which is newer than the original SymPy script (Apr 21), confirming a previous run remains the most recent execution; no fresh post-fix SymPy run was needed since neither finding targeted the SymPy script.

**Mathematica:** exit=0 (orchestrator-confirmed; `stage_002_mathematica.log` was not captured in `exec_logs/`, but the saved output `mathematica/output/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.txt` shows every assertion ending in `PASS:` with no `FAIL:` lines and the script exits via `Exit[0]`). Notable lines from the output:

- `PASS: extracted M matches boxed M (4 Pi overlap form)` (output line 24) — the new F1 coefficient-extraction assertion lands.
- `PASS: extracted K matches boxed K (4 Pi overlap form)` (line 27).
- `PASS: real P2 norm matrix - IdentityMatrix[5]` with a literal 5x5 zero matrix shown (lines 65-67).
- `PASS: real P2 angular stiffness matrix - 6 IdentityMatrix[5]` (lines 68-70).
- `PASS: norm(Y21s) - 1`, `PASS: angular energy(Y21s) - 6`, `PASS: -Delta_S2 Y21s - 6 Y21s`, plus the Y22s analogues (lines 53-56, 61-64, 75-76, 79-80).
- `PASS: single-component reduced density for Y21s`, `... for Y22s` (lines 85-86, 89-90).

**Output freshness:** `mathematica/output/moving_throat_pde_stage002_breathing_reduction_mathematica_audit.txt` mtime is 2026-05-21 11:50; the edited `.wl` mtime is 2026-05-21 00:39. Output is fresh (newer than the script). The SymPy output (2026-05-21 11:25) is also newer than the unchanged SymPy script.

## Material-change assessment

`material_change`: false.

The Mathematica edits change *how* the same physical claims are verified, not the claims themselves. The boxed `4 Pi mu_eta alpha^2`-type overlap-matrix formulas, the EL system `M_AB Q'' + K_AB Q == 0`, the 5x5 isotropy of the real P2 multiplet, and the single-component reduced equation `M2 q'' + K2 q == 0` are all the same statements as before; the second engine now arrives at them via `SphericalHarmonicY`, `Coefficient`, `EulerEquations`, and 5x5 `Table` integration rather than transliterating SymPy. No derived quantity downstream units depend on has changed. The SymPy script was not touched.

## Side observations (non-blocking)

1. After the F1 edits, the existing `formal reduced Lagrangian - boxed matrix form` check (`Lred - LredTarget` where `Lred = Integrate[Mintegrand, ...]` with `Mintegrand = MintegrandExtracted` and `LredTarget` uses `MintegrandBoxed`) is now downstream-redundant: once `MintegrandExtracted == MintegrandBoxed` is asserted, the integrated form equality follows trivially. It is harmless (still true, still informative as a sanity sentinel) but no longer carries independent information. Not a regression.

2. Section III's Y_2m literals (lines 208-212) are still hand-written rather than obtained via `SphericalHarmonicY[2, m, theta, phi]`. The F1 directive scoped Mathematica-native construct usage to Section I; expanding it to Section III is a possible future improvement but is not part of either finding.

3. The Mathematica exec log file `stage_002_mathematica.log` is absent from `redteam/exec_logs/`. The orchestrator stated the run was confirmed independently, and the post-fix output `.txt` (mtime 2026-05-21 11:50, post-edit) shows uniform PASS, so this is informational; the orchestrator may want to standardize log capture for the audit-trail.

## Verdict justification

Both F1 and F2 are resolved with substance. F1's three sub-changes each substitute a Mathematica-native primitive (`SphericalHarmonicY`, `Coefficient`-based matrix extraction, `EulerEquations`) for the transliterated SymPy algebra; codex's deviation from the directive's literal `1/(dadt^2)` and `1/2` factors was mathematically necessary and well-documented in the `## Applied:` block, and the resulting assertions are non-tautological because the extracted matrices come from independent angular integration while the boxed targets are postulated separately. F2's 5x5 matrix-form assertions, full-basis per-component checks, and 5-element `qvec` substantively close the 3-of-5 gap with direct off-diagonal integration. All assertions PASS in the fresh output, and no regressions are visible in the diff. Material change: false; the SymPy script and all derived numbers downstream units rely on are unchanged.
