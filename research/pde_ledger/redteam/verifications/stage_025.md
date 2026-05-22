---
unit_id: 025
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 8
findings_total: 8
material_change: false
---

# Verification — unit 025

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:62-64` defines `SAMPLE_POINT` at module scope with the directive's exact values.
- `scripts/...sympy_audit.py:110-117` keeps the legacy `expect_zero("P0 - P0_compact", ...)` at line 108 and adds the new `P0 raw at sample = 1/3` / `P0 compact at sample = 1/3` checks with explicit `raise AssertionError` if either rational is not `1/3`.
- `mathematica/...mathematica_audit.wl:36-37` defines `samplePoint`.
- `mathematica/...mathematica_audit.wl:70-76` adds the same sample-point evaluation and `fail[]` guard, ending with `pass["P0 numerical at sample point"]`.

**Assessment:**
The new checks are non-tautological: they substitute concrete integers and demand the rational evaluate to exactly `1/3`. A sign flip in `N0`, `D0`, or `P` would now break the test (the result would no longer be `1/3`). Saved sympy output line 23-24 shows `P0 raw at sample = 1/3` and `P0 compact at sample = 1/3`. Saved mathematica output shows `PASS: P0 numerical at sample point`. No collateral edits beyond the directive's scope.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:153-164` retains the two legacy `expect_zero` identities (151-152) and adds `Delta on sample`, `D0 on sample`, `P0 on sample` numerical reads with `raise AssertionError` if any is `<= 0`.
- `mathematica/...mathematica_audit.wl:109-118` mirrors with `FullSimplify[... /. samplePoint]` and `If[!TrueQ[... > 0], fail[...]]` guards, ending with `pass["Stability positivity on sample point"]`.

**Assessment:**
The positivity check is real: it converts the prose claim "Delta > 0 and D0 > 0 imply P0 > 0" into three executed `<=0` traps at a concrete point. Saved sympy output shows `Delta on sample = 15`, `D0 on sample = 1/3`, `P0 on sample = 1/3` (all match the directive's hand-computed expected values). Saved mathematica output shows `PASS: Stability positivity on sample point`. The hand-computed expectations (`Delta = 15`, `D0 = 1/3`, `P0 = 1/3`) match the output exactly — strong evidence the substitutions evaluated as intended.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:190-201` adds `sample_iv` (extending `SAMPLE_POINT` with `X = 1`) and asserts `dP0/dK < 0`, `dP0/dX > 0`, `dP0/dK + dP0/dX == 0` numerically.
- `mathematica/...mathematica_audit.wl:140-148` mirrors with `Join[samplePoint, {x -> 1}]` and `If[!TrueQ[... < 0], fail[...]]` / `If[!TrueQ[... > 0], fail[...]]` / sum-to-zero guard, ending with `pass["Monotonic sign checks on sample point"]`.

**Assessment:**
The new sign checks are non-tautological: they convert the prose "monotonic" claim into executed numerical-sign assertions. Saved sympy output shows `dP0/dK on sample = -1` and `dP0/dX on sample = 1` (matching the directive's hand-computed expected values). Saved mathematica output shows `PASS: Monotonic sign checks on sample point`. A sign flip in `N0` would now flip both derivative signs and trip the `<= 0` / `>= 0` traps.

### F4 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:119-122` adds a comment block above `target = ...` citing `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:321-342` and naming the upstream derivation chain (`Gamma5_port = a^5/(27*c_s^5)`, `gamma_GR = 2*G/(5*c^5)`, ratio = `54*G*c_s^5/(5*a^5*c^5)`).
- `mathematica/...mathematica_audit.wl:78-83` adds the analogous `(* ... *)` comment block with the same citation.

**Assessment:**
I read `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:315-343` directly to verify the citation. Lines 333 (`expect_zero("Stage-5 Gamma5_port anchor", Gamma5_port - a**5 / (27 * c_s**5))`), 334 (`gamma_GR = sp.simplify(2 * G / (5 * c**5))`), and 340-342 (`expect_zero("ratio target at mhat=1", ratio_target.subs(mhat, 1) - 54 * G * c_s**5 / (5 * a**5 * c**5))`) confirm the citation is accurate and the cited upstream stage genuinely derives the 54/5 coefficient. The citation is non-tautological (a reader can trace and check), and an asserted upstream derivation (not just stated) exists.

### F5 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/...sympy_audit.py:129-138` adds a solvability assertion: computes `mhat_sq = target / P0_compact`, substitutes `SAMPLE_POINT` and `{G, c_s, a, c} -> 1`, and asserts `mhat_sq_at_sample > 0`.
- `mathematica/...mathematica_audit.wl:89-100` adds both the sample-point positivity assertion AND the F8-required `Reduce` solvability check.

**Assessment:**
The directive's expected value `mhat^2 = 162/5` matches the saved output exactly in both engines (sympy output line 43, mathematica output `mhat^2 on sample = 162/5`). The check is non-tautological: it forms the ratio `target/P0_compact` (which would be negative if a sign error existed in either side of the target equation) and asserts positivity at a definite rational. A sign flip in `54/5` or in `P0_compact` would now break this check.

### F6 — insufficient_verification

**Classification:** resolved

**What changed:**
No separate edits — closed by F2's positivity check, as the directive specified. The `## Applied: F6` block in the directive correctly references the F2 patch with `summary: closed by F2 patch; no separate edits`.

**Assessment:**
F2's sample-point positivity assertions (`Delta > 0`, `D0 > 0`, `P0 > 0` at the sample point) directly test the prose claim from F6 ("If Delta > 0 and D0 > 0, then P0 > 0 whenever P != 0"). The structural overlap is exactly what the auditor anticipated.

### F7 — insufficient_verification

**Classification:** resolved

**What changed:**
No separate edits — closed by F3's sign checks, as the directive specified. The `## Applied: F7` block in the directive correctly references the F3 patch with `summary: closed by F3 patch; no separate edits`.

**Assessment:**
F3's sample-point sign assertions (`dP0/dK < 0`, `dP0/dX > 0`) directly test the prose "monotonic" claim from F7. Coverage overlap is correct.

### F8 — mathematica_transliteration

**Classification:** resolved

**What changed (mathematica only):**
- `mathematica/...mathematica_audit.wl:41`: `delta = Factor[omegaU^2*omegaW^2 - r^2]` (was `FullSimplify[omegaU^2*omegaW^2 - r^2, ...]`). Output now shows `Delta = (omegaU*omegaW - r)*(omegaU*omegaW + r)` — structurally different from SymPy's `OmegaU**2*OmegaW**2 - R**2`.
- `mathematica/...mathematica_audit.wl:62-63`: `p0Combined = Together[n0/d0]; p0Compact = Apart[p0Combined, k]` (was a hand-written compact form mirroring the SymPy literal). Output shows the Apart-canonical form using the factored Delta, not the SymPy-style nested form.
- `mathematica/...mathematica_audit.wl:126-129`: `dP0dK = Limit[(p0 /. k -> k + h) - p0)/h, h -> 0]` (was `D[p0, k]`). Lines 128-129 keep a `D[p0, k]` computation and lines 135-136 add `expectZero["Limit dP0/dK - D[p0,k]", dP0dK - dP0dKDirect]` as a second-source cross-check. Both `PASS`.
- `mathematica/...mathematica_audit.wl:40, 59, 65, 87, 104, 123, 131`: banner strings rewritten from "SECTION X — TITLE" em-dash to "Section X: title" colon form. No verbatim match to the SymPy banner strings remains.
- `mathematica/...mathematica_audit.wl:94-100`: adds `Reduce[mhat^2 == target/p0Compact && mhat > 0 && delta > 0 && d0 > 0, mhat, Reals]` and asserts the result is not `False`. This is a Mathematica-native primitive not used on the SymPy side.

**Assessment:**
All four required Mathematica primitives are present and exercised (`Factor`, `Apart`, `Limit`, `Reduce`). The cross-form check `P0 - P0_compact == 0` (line 68) is now a real second-source check rather than a tautology: SymPy's `simplify(N0/D0)` and Mathematica's `Apart[Together[n0/d0], k]` reach the same canonical content through different machinery, and the engines now agree on the canonical form, not by transliteration but by independent canonicalization. The `Limit dP0/dK - D[p0,k] = 0` and `Limit dP0/dX - D[p0,x] = 0` cross-checks (PASS in the output) are real second-source derivative confirmations absent in the SymPy side.

The `Reduce` solvability check returned a non-False conditional (a long disjunction over sign cases), so `solvability === False` is not triggered and `pass["II.2 target equation solvable for mhat > 0"]` fires. Note: `$Assumptions` is not consulted by `Reduce` (only by `FullSimplify`), so Reduce explores the full sign lattice — this is acceptable because the inline constraints `mhat > 0 && delta > 0 && d0 > 0` are what the test needs, and Reduce's task is to confirm that *some* parameter regime makes the equation solvable for `mhat > 0`. If the target had a wrong sign such that `mhat^2 = target/p0Compact` were forced negative on the stability branch, Reduce would return `False` and the assertion would fire. The check is non-vacuous.

## Exec log assessment

**SymPy:** exit log not captured by the orchestrator at `/var/projects/toy_physics/research/pde_ledger/redteam/exec_logs/stage_025_sympy.log` (file absent). However, the saved output `scripts/output/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.txt` is fresher than the script (mtime 1779404313 > 1779404129) and shows all expected new lines:
- line 22: `P0 - P0_compact = 0`
- lines 23-24: `P0 raw at sample = 1/3` / `P0 compact at sample = 1/3`
- line 43: `mhat^2 on sample = 162/5`
- lines 61-63: `Delta on sample = 15` / `D0 on sample = 1/3` / `P0 on sample = 1/3`
- lines 89-90: `dP0/dK on sample = -1` / `dP0/dX on sample = 1`
None of the `raise AssertionError` traps fired (the script reached the end and emitted Section IV output). Treating saved-output success as `exit=0`.

**Mathematica:** exit log not captured at `redteam/exec_logs/stage_025_mathematica.log` (file absent). Saved output `mathematica/output/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.txt` is fresher than the script (mtime 1779404322 > 1779404200) and contains every required `PASS:` line:
- `PASS: P0 - P0_compact`
- `PASS: P0 numerical at sample point`
- `PASS: II.2 target equation solvability at sample point`
- `PASS: II.2 target equation solvable for mhat > 0`
- `PASS: Delta*D0 - (K*Delta - Delta*C^2/varpi^2 - Q)`
- `PASS: N0 - P^2/Delta^2`
- `PASS: Stability positivity on sample point`
- `PASS: Limit dP0/dK - D[p0,k]`
- `PASS: Limit dP0/dX - D[p0,x]`
- `PASS: dP0/dK + N0/(K - X - Q/Delta)^2`
- `PASS: dP0/dX - N0/(K - X - Q/Delta)^2`
- `PASS: dP0/dX + dP0/dK`
- `PASS: Monotonic sign checks on sample point`
The trailing line `Stage 8 Mathematica audit passed.` confirms `Exit[0]` was reached. Treating saved-output success as `exit=0`.

**Output freshness:** confirmed. Both saved `.txt` outputs have mtimes strictly newer than their corresponding script mtimes (sympy: 1779404313 > 1779404129; mathematica: 1779404322 > 1779404200).

## Material-change assessment

`material_change`: false.

The edits add positivity, sign, and solvability assertions at a concrete sample point; they do not alter any symbolic form of `P0`, `Delta`, `D0`, `N0`, the target coefficient, or the derivative identities that downstream units might consume. The Mathematica-side rewrite of `delta = Factor[...]`, `p0Compact = Apart[...]`, and limit-form derivatives produces canonical forms that are algebraically equal to (not different from) the prior forms — verified by the `expectZero` cross-checks all passing. No downstream unit re-audit triggered by this change.

## Side observations (non-blocking)

- The `Reduce` solvability check's `$Assumptions` are not honored (Reduce uses inline constraints only), so the printed result is a comprehensive case analysis. The check itself works as designed (asserts `=!= False`), but the printed output is ~100KB. This is not a finding — it is a known Mathematica idiom — but if future audits want a more focused solvability witness, passing `delta > 0 && d0 > 0` as `Element` constraints or pre-substituting positive-parameter assumptions into the equation before `Reduce` would shorten the output.
- The exec logs at `redteam/exec_logs/stage_025_sympy.log` and `redteam/exec_logs/stage_025_mathematica.log` are absent. The orchestrator should populate these — current verification leans on saved `.txt` outputs which are fresh. Not a blocking issue for this finding set.

## Verdict justification

Every finding in the original report is `resolved`. Codex applied each directive verbatim, added the directive-specified numerical sample-point checks in both engines, restructured the Mathematica witness with independent primitives (`Factor`, `Apart`, `Limit`, `Reduce`) and distinct banners, and cited the upstream stage where `54/5` is derived (verified by direct read of stage 023). The saved outputs show every required new line at the directive's hand-computed expected values (1/3, 15, 162/5, -1, +1) and every required `PASS:` marker; no `AssertionError` or `FAIL` appears. Output mtimes confirm the saved files are post-fix. No regressions visible in the diff. No collateral edits beyond the directive's scope.
