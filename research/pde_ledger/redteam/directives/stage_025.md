---
unit_id: 025
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 8
stop_cold: null
applied: true
applied_at: 2026-05-21T22:56:53Z
findings_applied: 8
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 025

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## Shared sample point (used by F1, F2, F3, F5)

Both engines must use the same isotropic sample point:
```
K = 2,  varpi = 1,  C = 1,  OmegaU = 2,  OmegaW = 2,  R = 1,  GU = 1,  GW = 1
```
At this point: `Delta = 15`, `Q = 10`, `P = 5`, `B0 = 1`, `Z0 = 2/3`, `D0 = 1/3`, `N0 = 1/9`, `P0 = 1/3`, and with `X = 1`: `dP0/dK = -1`, `dP0/dX = +1`. Codex must store this dict as a module-level constant in each script (named `SAMPLE_POINT` in sympy, `samplePoint` association in Mathematica) and reuse it across F1/F2/F3/F5.

If any reused assertion below substitutes the sample point and the resulting rational does NOT match the expected value listed above, the patch is wrong — recompute by hand before applying.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:97-104`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:58-65`

**Issue:**
Section II.1 checks `P0 - P0_compact == 0` where `P0 = simplify(N0/D0)` and `P0_compact = simplify(P^2/(Delta*(K*Delta - Delta*C^2/varpi^2 - Q)))`. Since `D0 * Delta = K*Delta - Delta*C^2/varpi^2 - Q` by the definition of D0 at line 77, the two forms are algebraically identical by construction; the check is tautological.

**Required change:**

SymPy (script `..._sympy_audit.py`):
1. Immediately after line 99 (`P0_compact = sp.simplify(...)`), insert a constant `SAMPLE_POINT` dictionary (or use the module-level constant defined per the Shared Sample Point section above):
   ```python
   SAMPLE_POINT = {K: sp.Integer(2), varpi: sp.Integer(1), C: sp.Integer(1),
                   OmegaU: sp.Integer(2), OmegaW: sp.Integer(2), R: sp.Integer(1),
                   GU: sp.Integer(1), GW: sp.Integer(1)}
   ```
   (place `SAMPLE_POINT` at module scope just below the symbol declarations near line 60 instead of inside the function, so F2/F3/F5 can reuse it).
2. Keep the existing `expect_zero("P0 - P0_compact", P0 - P0_compact)` at line 104 as a redundant identity check, but immediately *after* it add:
   ```python
   p0_value_raw = sp.nsimplify(P0.subs(SAMPLE_POINT))
   p0_value_compact = sp.nsimplify(P0_compact.subs(SAMPLE_POINT))
   print(f"P0 raw at sample     = {p0_value_raw}")
   print(f"P0 compact at sample = {p0_value_compact}")
   if p0_value_raw != sp.Rational(1, 3):
       raise AssertionError(f"P0 raw at sample point != 1/3: got {p0_value_raw}")
   if p0_value_compact != sp.Rational(1, 3):
       raise AssertionError(f"P0 compact at sample point != 1/3: got {p0_value_compact}")
   ```

Mathematica (script `..._mathematica_audit.wl`):
1. Define `samplePoint` at top-of-file scope (just below the `$Assumptions` definition near line 35):
   ```mathematica
   samplePoint = {k -> 2, varpi -> 1, cCoupling -> 1, omegaU -> 2, omegaW -> 2,
                  r -> 1, gU -> 1, gW -> 1};
   ```
2. After the `expectZero["P0 - P0_compact", ...]` call at wl:65, add:
   ```mathematica
   p0Raw = (p0 /. samplePoint);
   p0Compactv = (p0Compact /. samplePoint);
   Print["P0 raw at sample     = ", fmt[p0Raw]];
   Print["P0 compact at sample = ", fmt[p0Compactv]];
   If[p0Raw =!= 1/3, fail["P0 raw at sample point != 1/3", p0Raw]];
   If[p0Compactv =!= 1/3, fail["P0 compact at sample point != 1/3", p0Compactv]];
   pass["P0 numerical at sample point"];
   ```

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 025` and `redteam exec-mathematica 025`. New lines `P0 raw at sample     = 1/3` and `P0 compact at sample = 1/3` must appear in both transcripts, followed by a PASS marker (sympy: no AssertionError; mathematica: the explicit `pass[]` line). Both scripts must exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Added the shared sample point constants and asserted raw and compact P0 both evaluate to 1/3 at the sample.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:118-126`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:74-81`

**Issue:**
Section III's two `expect_zero` calls are algebraic identities by the definitions at sympy lines 76-77 (and wl lines 44-45). The "positivity structure" advertised in the section banner and the line-126 print is never asserted.

**Required change:**

SymPy:
1. Keep the existing two `expect_zero` lines (123-124) for legacy.
2. After line 124, before line 126, add:
   ```python
   delta_value = sp.nsimplify(Delta.subs(SAMPLE_POINT))
   d0_value = sp.nsimplify(D0.subs(SAMPLE_POINT))
   p0_pos = sp.nsimplify((N0 / D0).subs(SAMPLE_POINT))
   print(f"Delta on sample = {delta_value}")
   print(f"D0    on sample = {d0_value}")
   print(f"P0    on sample = {p0_pos}")
   if delta_value <= 0:
       raise AssertionError(f"Delta on sample is not positive: {delta_value}")
   if d0_value <= 0:
       raise AssertionError(f"D0 on sample is not positive: {d0_value}")
   if p0_pos <= 0:
       raise AssertionError(f"P0 on sample is not positive: {p0_pos}")
   ```
   Expected values: `Delta = 15`, `D0 = 1/3`, `P0 = 1/3` — all positive.

Mathematica:
1. Keep the existing two `expectZero` calls (wl:78-79) for legacy.
2. After wl:79, before wl:80, add:
   ```mathematica
   deltaVal = (delta /. samplePoint);
   d0Val = (d0 /. samplePoint);
   p0PosVal = ((n0/d0) /. samplePoint);
   Print["Delta on sample = ", fmt[deltaVal]];
   Print["D0    on sample = ", fmt[d0Val]];
   Print["P0    on sample = ", fmt[p0PosVal]];
   If[!TrueQ[deltaVal > 0], fail["Delta on sample is not positive", deltaVal]];
   If[!TrueQ[d0Val > 0], fail["D0 on sample is not positive", d0Val]];
   If[!TrueQ[p0PosVal > 0], fail["P0 on sample is not positive", p0PosVal]];
   pass["Stability positivity on sample point"];
   ```

**Verification command:**
After Codex applies, transcripts must show three new numeric lines (`Delta on sample = 15`, `D0 on sample = 1/3`, `P0 on sample = 1/3`) and a positivity PASS in each engine. Scripts exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Added sample-point positivity checks for Delta, D0, and P0 in both engines.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:133-148`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:83-96`

**Issue:**
The three derivative identities are forced by the explicit form `P0 = N0/(K - X - Q/Delta)` at sympy:137 (wl:86). The "monotonic" sign claim in the section banner is not tested.

**Required change:**

SymPy:
1. Keep the existing three `expect_zero` calls (146-148) for legacy.
2. Extend `SAMPLE_POINT` for the IV-section symbol `X`: just before the new sign checks add
   ```python
   sample_iv = dict(SAMPLE_POINT)
   sample_iv[X] = sp.Integer(1)  # X = C^2/varpi^2 = 1 at sample
   ```
3. After line 148, add:
   ```python
   dP0_dK_val = sp.nsimplify(dP0_dK.subs(sample_iv))
   dP0_dX_val = sp.nsimplify(dP0_dX.subs(sample_iv))
   print(f"dP0/dK on sample = {dP0_dK_val}")
   print(f"dP0/dX on sample = {dP0_dX_val}")
   if dP0_dK_val >= 0:
       raise AssertionError(f"dP0/dK on sample is not negative: {dP0_dK_val}")
   if dP0_dX_val <= 0:
       raise AssertionError(f"dP0/dX on sample is not positive: {dP0_dX_val}")
   if dP0_dK_val + dP0_dX_val != 0:
       raise AssertionError(f"dP0/dK + dP0/dX on sample is not zero: {dP0_dK_val + dP0_dX_val}")
   ```
   Expected: `dP0/dK = -1`, `dP0/dX = +1`, sum = 0.

Mathematica:
1. Keep the existing three `expectZero` calls (wl:94-96) for legacy.
2. After wl:96, add:
   ```mathematica
   sampleIV = Join[samplePoint, {x -> 1}];
   dP0dKVal = (dP0dK /. sampleIV);
   dP0dXVal = (dP0dX /. sampleIV);
   Print["dP0/dK on sample = ", fmt[dP0dKVal]];
   Print["dP0/dX on sample = ", fmt[dP0dXVal]];
   If[!TrueQ[dP0dKVal < 0], fail["dP0/dK on sample is not negative", dP0dKVal]];
   If[!TrueQ[dP0dXVal > 0], fail["dP0/dX on sample is not positive", dP0dXVal]];
   If[dP0dKVal + dP0dXVal =!= 0, fail["dP0/dK + dP0/dX on sample is not zero", dP0dKVal + dP0dXVal]];
   pass["Monotonic sign checks on sample point"];
   ```

**Verification command:**
After Codex applies, transcripts must show `dP0/dK on sample = -1` and `dP0/dX on sample = 1` and a monotonic-sign PASS. Scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Added sample-point derivative sign checks for dP0/dK and dP0/dX in both engines.
- deviation: none

## F4 — hardcoded_result

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:106`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:67`

**Issue:**
The constant `54/5` in `target = 54*G*c_s^5/(5*a^5*c^5)` is pre-baked. Neither engine derives it; neither cites the upstream stage where it is derived.

**Required change:**

SymPy: insert a comment block immediately above line 106:
```python
# target coefficient 54/5 = ?
# This coefficient is carried forward from upstream and not derived in this script.
# Codex / red-team note: add a citation here pointing to the upstream stage where
# 54/5 is computed (file path and line). If no such derivation exists upstream,
# raise a separate finding and downgrade this Section II.2 to a definition rather
# than a "verification".
```
If Codex can locate the upstream derivation by reading other unit scripts (e.g. via grep for `54` in `scripts/moving_throat_pde_stage*.py`), replace the placeholder with the actual `file:line` citation. If not located, leave the placeholder text and emit a `## Blocked: F4` directive entry explaining "no upstream derivation of 54/5 found in `scripts/`."

Mathematica: insert the analogous comment block above wl:67. Use `(*  ...  *)` Mathematica comment syntax.

**Verification command:**
After Codex applies, both scripts must contain a comment block above the `target = ...` line. The verifier inspects the diff for the new comment. If Codex blocked, the verifier reads the `## Blocked: F4` reasoning before marking F4 closed.

## Applied: F4

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Added comments citing `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:321-342` as the upstream derivation of the 54/5 target coefficient.
- deviation: none

## F5 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py:107-111`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl:67-71`

**Issue:**
Section II.2 prints `equation_residual` but never asserts. The single physically meaningful claim in the script (the target equation) is decoratively printed only.

**Required change:**

SymPy: after line 111 (the `sp.pprint(equation_residual)`), add:
```python
# Solvability check: mhat^2 = target / P0_compact must be positive on the
# stability branch (Delta > 0, D0 > 0, so P0_compact > 0; and target > 0
# since G, c_s, a, c > 0). On the sample point this is a definite rational.
mhat_sq = sp.simplify(target / P0_compact)
mhat_sq_at_sample = sp.nsimplify(mhat_sq.subs(SAMPLE_POINT))
print(f"mhat^2 on sample = {mhat_sq_at_sample}")
if mhat_sq_at_sample <= 0:
    raise AssertionError(f"mhat^2 on sample is not positive: {mhat_sq_at_sample}")
```
(`target` and `P0_compact` are in scope inside `normalization_formula`; if `SAMPLE_POINT` was defined module-scope per F1, it is accessible here.)

Expected: at the F1 sample point with `mhat`, `G`, `c_s`, `a`, `c` still symbolic, `mhat_sq.subs(SAMPLE_POINT)` is a symbolic rational in those variables. Codex must therefore also substitute `G -> 1, c_s -> 1, a -> 1, c -> 1` for the numeric assertion (extend the dict or do an additional `.subs({G: 1, c_s: 1, a: 1, c: 1})` before `nsimplify`). With those substitutions: `target = 54/5`, `P0_compact = 1/3`, so `mhat^2 = (54/5) / (1/3) = 162/5 = 32.4`. Expected output line: `mhat^2 on sample = 162/5`.

Mathematica: after wl:71 (the `Print["Target residual = ", ...]`), add:
```mathematica
mhatSq = FullSimplify[target / p0Compact, Assumptions -> $Assumptions];
mhatSqSample = mhatSq /. Join[samplePoint, {gConst -> 1, cs -> 1, a -> 1, cSpeed -> 1}];
Print["mhat^2 on sample = ", fmt[mhatSqSample]];
If[!TrueQ[mhatSqSample > 0], fail["mhat^2 on sample is not positive", mhatSqSample]];
pass["II.2 target equation solvability at sample point"];
```
Expected output line: `mhat^2 on sample = 162/5`.

**Verification command:**
After Codex applies, both transcripts must contain `mhat^2 on sample = 162/5` and a PASS. Scripts exit 0.

## Applied: F5

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Added positive solvability assertions for mhat^2 at the shared sample point in both engines.
- deviation: none

## F6 — insufficient_verification

**Target:** absorbed by F2. No additional Codex edits required.

**Issue:** Same gap as F2 (positivity claim prose-only).

**Required change:** none beyond F2. Codex must add an `## Applied: F6` block referencing the F2 patch with `summary: closed by F2 patch; no separate edits` and `deviation: none`.

**Verification command:** same as F2.

## Applied: F6

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: closed by F2 patch; no separate edits
- deviation: none

## F7 — insufficient_verification

**Target:** absorbed by F3. No additional Codex edits required.

**Issue:** Same gap as F3 (monotonicity sign prose-only).

**Required change:** none beyond F3. Codex must add an `## Applied: F7` block referencing the F3 patch with `summary: closed by F3 patch; no separate edits` and `deviation: none`.

**Verification command:** same as F3.

## Applied: F7

- files_changed:
  - `scripts/moving_throat_pde_stage025_minimal_isotropic_normalization_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: closed by F3 patch; no separate edits
- deviation: none

## F8 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl` (whole file)

**Issue:**
The `.wl` file's banner strings, term orderings, module decomposition, and assertion list are line-by-line ports of the `.py` file. No Mathematica-native primitive (`Factor`, `Apart`, `Limit`-based derivative definition, `Reduce`) is used. Both engines walk the same path and therefore cannot catch each other's translation errors.

**Required change:**
Rewrite the Mathematica `.wl` so it independently computes the four sectional results from the same stated symbolic inputs but does NOT mirror the SymPy term ordering or compact-form choice. Specifically:

1. At wl:39, replace `delta = FullSimplify[omegaU^2*omegaW^2 - r^2, ...]` with
   ```mathematica
   delta = Factor[omegaU^2*omegaW^2 - r^2];
   ```
   Result will be `(omegaU*omegaW - r)*(omegaU*omegaW + r)`, structurally different from the SymPy unfactored form. Use this factored form throughout the wl file.

2. At wl:60, replace the hand-written compact form with a derivation:
   ```mathematica
   p0Combined = Together[n0/d0];
   p0Compact = Apart[p0Combined, k];
   ```
   Then verify `FullSimplify[p0 - p0Compact] === 0` (this becomes a real cross-form check rather than a tautology, because `Apart` returns a Mathematica-chosen canonical form, not the SymPy hand-written one).

3. At wl:87-88, replace `D[p0, k]` and `D[p0, x]` with limit-of-quotient definitions:
   ```mathematica
   dP0dK = FullSimplify[Limit[((p0 /. k -> k + h) - p0)/h, h -> 0]];
   dP0dX = FullSimplify[Limit[((p0 /. x -> x + h) - p0)/h, h -> 0]];
   ```
   Then add `expectZero` calls comparing these against `D[p0, k]` and `D[p0, x]` to confirm the two derivative definitions agree (a real second-source check, not present in the SymPy script).

4. In the banners, change "SECTION I — MINIMAL ISOTROPIC ZERO-FREQUENCY COEFFICIENTS" etc. to either use a different separator (e.g., colon instead of em-dash) OR restate the section title in different wording (e.g., "Zero-Frequency Coefficients (Section I)"). This is not stylistic noise: it documents that the script was authored, not transliterated. Apply consistently to all four section banners (wl:38, wl:57, wl:62, wl:70, wl:75, wl:84, wl:90).

5. After F5's Section II.2 solvability check, add an independent re-derivation via `Reduce`:
   ```mathematica
   solvability = Reduce[mhat^2 == target/p0Compact && mhat > 0 && delta > 0 && d0 > 0,
                        mhat, Reals];
   Print["II.2 solvability (Reduce form) = ", fmt[solvability]];
   If[solvability === False, fail["II.2 target equation unsolvable for mhat > 0"]];
   pass["II.2 target equation solvable for mhat > 0"];
   ```
   This is a Mathematica-native test absent from the Python script — it makes engine agreement on II.2 a real signal rather than parallel silence.

**Verification command:**
After Codex applies, the verifier (a) greps the `.wl` file to confirm no line matches the SymPy banner strings verbatim (compare to the `.py` file's literal banner strings); (b) confirms `Factor`, `Apart`, and `Limit` each appear at least once in the `.wl`; (c) runs `redteam exec-mathematica 025` and confirms the new `II.2 solvability` PASS line appears. The script must exit 0.

If after the rewrite the engines disagree (Mathematica's `Apart`-canonical form does not match the SymPy compact form after `FullSimplify`), that is a real `engine_disagreement` finding to surface in the next audit pass — do NOT mask it by re-mirroring SymPy's form.

## Applied: F8

- files_changed:
  - `mathematica/moving_throat_pde_stage025_minimal_isotropic_normalization_mathematica_audit.wl`
- summary: Reworked the Mathematica witness to use factored Delta, Apart-derived P0, quotient-limit derivatives, distinct banners, and Reduce solvability.
- deviation: none
