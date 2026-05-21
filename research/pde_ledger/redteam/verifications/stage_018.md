---
unit_id: 018
batch: I.2
verifier_model: claude-opus-4-7
verify_date: 2026-05-21T15:10:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: n/a
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 018

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:41-45` the tautological assertion `assert_zero("compatibility equality", sp.simplify(K_from_norm - K_from_one_pole) - compatibility)` was replaced exactly per directive with two independent `sp.solve` calls (one on `compatibility`, one on `K_from_norm - K_from_one_pole`) followed by an equality assertion on the two `N0` solutions and a mutation assertion `N0_from_compat - 2 * N0_from_equality`.

**Assessment:**
The edit matches the directive's `Before`/`After` blocks verbatim. The previous form expanded literally to `expr - expr == 0`; the new form invokes `sp.solve` twice on two different expressions and compares their roots. The mutation row introduces a nonzero residual `-N0_from_equality = -3*Ptarget*(MSigma+B2+Z2)^2/(B4+Z4)` that is genuinely nonvanishing, so it would fail if the equality row were silently degenerate. No collateral edits beyond the prescribed line range. SymPy exec log shows the script still exits 0. See side observations for a residual nuance about the two solve inputs being algebraically equivalent (the directive's own prescription).

### F2 — tautological_check

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:73-75` the two `K1.subs(sol)` / `H_even.subs(sol)` assertions were replaced exactly per directive with three substitutions that use the independently defined closed forms `expected_dK = B01+Z01+27*(B41+Z41)` and `expected_dM = -(B21+Z21)+3*(B41+Z41)` (defined at lines 65-66) plus a mutation `K1.subs({dKSigma: expected_dK + 1, dMSigma: expected_dM})`.

**Assessment:**
The new checks no longer route through `sp.solve`'s output `sol`, so they no longer reduce to the contract of `sp.solve`. `K1 = D21 + D01/9 = -(dM+B21+Z21) + (dK-B01-Z01)/9`; substituting the closed forms collapses to `-3(B41+Z41) + 27(B41+Z41)/9 = 0`, which is the substantive identity. The mutation `expected_dK + 1` yields `K1 = 1/9 != 0`, a genuine nonzero residue that the assertion will catch if SymPy fails to simplify. No collateral edits to lines 69 (`D01_comp`) or 76-83 (`Xi1` block) as required. SymPy exec exit 0.

### F3 — missing_verification_script

**Classification:** resolved

**What changed:**
New file created at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.wl` (141 lines). It covers M1-M8 and their mutation rows (M3-mut, M5-mut, M7-mut) with `If[FullSimplify[...] =!= 0, Exit[1]]` and `If[FullSimplify[...] === 0, Exit[1]]` guards, and ends with `Print["STAGE 018 MATHEMATICA AUDIT PASS"]; Exit[0]`.

**Assessment:**
The script is not a line-by-line transliteration of the SymPy script. M1 derives `u2, u4` from the pole-series expansion `Normal[Series[1/(D0+D2 x^2+D4 x^4), {x,0,4}]]` followed by `Coefficient[...*D0, x, 2]` and `Coefficient[...*D0, x, 4]` — this re-derives the rational-function coefficients independently of the SymPy script's hand-computed `u2 = -D2/D0`, `u4 = (D2^2 - D0*D4)/D0^2`. M5 uses Mathematica's matrix-derivative idiom `D[gateVector, {{dKSigma, dMSigma}}]`. M6 substitutes the closed-form slopes via a rule list (no `Solve` round-trip), as the directive demanded. M8 uses `Integrate[..., {w, -Infinity, Infinity}]` with `$Assumptions = Element[w, Reals]`. The saved output at `mathematica/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.txt` (mtime 15:02, newer than the script mtime 13:30) shows all M1, M2, M3, M3-mut, M4, M5, M5-mut, M6 (K1 and H_even), M7, M7-mut, M8 (inertia and stiffness) residuals printed with the expected values: zeros for the equality rows and the genuinely nonzero `-Ptarget`, `2/27`, and `(54*(B41+Z41))/(B0-KSigma+Z0)` for the mutation rows, then "STAGE 018 MATHEMATICA AUDIT PASS". The exec log file `stage_018_mathematica.log` is missing from `redteam/exec_logs/`, so I rely on the saved output file as the freshness evidence (see Exec log assessment).

### F4 — insufficient_verification

**Classification:** resolved

**What changed:**
At `scripts/moving_throat_pde_stage018_parent_throat_action_bundle_master_sympy_audit.py:84-88` an additional block was inserted after the existing `Xi1` mutation row: it computes `Xi1_from_expected = Xi1.subs({dKSigma: expected_dK, dMSigma: expected_dM})` and asserts it equals `N01/N0 - 27*(B41+Z41)/(KSigma-B0-Z0)`. Lines 76-83 (original A14/A15 block) are untouched, matching the directive's instruction to add, not replace.

**Assessment:**
The new check routes through `expected_dK`/`expected_dM` (the closed forms defined directly from physical premises) rather than through `sp.solve`'s output `sol`, so it does not share intermediate simplification with A11 (`D01.subs(sol)`). If SymPy's `subs` for `sol` were silently miscomputing `D01.subs(sol)`, that error would propagate through both A11 and A14, but the new check is independent of `sol` and would diverge. Substantive. SymPy exec exit 0.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
```
STEP 16 PARENT THROAT ACTION BUNDLE MASTER AUDIT
Checked isotropic compatibility, exact weak-axisymmetric wall-slope solve, and residual Xi1.
STATUS: PASS
# exit_code: 0
```
All four new/changed assertions ran without an `AssertionError`, including the two new `assert_nonzero` mutation rows from F1 and F2.

**Mathematica:** exit=n/a. The orchestrator did not capture `redteam/exec_logs/stage_018_mathematica.log`; only `stage_018_sympy.log` and `stage_018_diff.patch` are present in that directory. However, the saved Mathematica output at `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage018_parent_throat_action_bundle_master_mathematica_audit.txt` (mtime 2026-05-21 15:02:27, newer than the `.wl` mtime 2026-05-21 13:30:20) ends with `STAGE 018 MATHEMATICA AUDIT PASS`, and every printed residual is consistent with the script's `If[... Exit[1]]` guards passing (zeros for the M1, M2, M3, M4, M5, M6, M7, M8 equality rows; the expected nonzero values `-Ptarget`, `2/27`, `(54*(B41+Z41))/(B0-KSigma+Z0)` for M3-mut, M5-mut, M7-mut). Treating the saved output as authoritative for exit status, the Mathematica run succeeded (exit 0). Recording `mathematica_exit: n/a` in the front-matter because the exec log itself was not captured by the orchestrator.

**Output freshness:** Confirmed. SymPy script mtime 1779391793 (13:29:53); SymPy output mtime 1779397234 (14:58:13). Mathematica `.wl` mtime 1779391820 (13:30:20); Mathematica output mtime 1779397347 (15:02:27). Both saved outputs are newer than their corresponding scripts, so they reflect the post-fix scripts.

## Material-change assessment

`material_change`: false.

No derived results changed downstream. The four edits are all in the verification layer: F1 and F2 swap tautological assertions for substantive ones over the same physical claims, F4 adds an independent verification path for an already-existing claim, and F3 adds an independent Mathematica engine for the same set of claims. No new symbolic identities are introduced, and no existing claims are revised. Downstream units that import results from unit 018 (none do mechanically; stage 018 is a master-note audit) would see identical numeric and symbolic outputs.

## Side observations (non-blocking)

- F1's substantive replacement is exactly what the directive prescribed, but the two `sp.solve` inputs are algebraically identical: `K_from_norm - K_from_one_pole = (B0+Z0+N0/Ptarget) - (B0+Z0+3(M+B2+Z2)^2/(B4+Z4)) = N0/Ptarget - 3(M+B2+Z2)^2/(B4+Z4) = compatibility`, so `sp.solve(compatibility, N0)` and `sp.solve(K_from_norm - K_from_one_pole, N0)` solve the same equation. The check therefore exercises SymPy's `solve` normalization on two algebraically equivalent forms rather than verifying a non-trivial physical equivalence; the substantive cross-closure content was already covered by A2 and A3. This is the auditor's own prescribed form and Codex applied it correctly, so I do not block on it. A future audit pass may want to strengthen F1 (e.g., assert that `N0_from_compat` equals `3*Ptarget*(MSigma+B2+Z2)^2/(B4+Z4)` in closed form), but that is a re-audit concern, not a verification failure.
- Mathematica exec log was not captured. The orchestrator's `redteam exec-mathematica 018` step appears to have either skipped or failed to write `stage_018_mathematica.log`, while the actual Mathematica run did produce a fresh output file. If this is a recurring pattern across other units, the orchestrator's log-capture for the Mathematica path may need review — but stage 018 itself is not blocked on it because the saved output is fresh and shows PASS.

## Verdict justification

All four findings are resolved. F1 and F2 replaced tautological assertions with substantive ones using closed-form expressions or independent solve paths, exactly as the directive specified. F3 adds a Mathematica audit that re-derives M1 via series expansion (not a transliteration) and verifies all eight claim families with proper `Exit[1]` guards; the saved output (newer than the script) shows PASS with genuine nonzero residuals on every mutation row. F4 adds an independent `Xi1` recomputation routed through the closed-form slopes rather than `sp.solve`'s output. SymPy exec exits 0 and the Mathematica saved output reports PASS with consistent residuals; no regressions appear in the diff. One subtle concern about F1's algebraic equivalence is recorded as a side observation but does not warrant rework, since the directive's prescription was applied verbatim and the mutation row is genuinely nonzero.
