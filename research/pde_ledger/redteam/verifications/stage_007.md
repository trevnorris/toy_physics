---
unit_id: 007
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

# Verification — unit 007

## Per-finding outcomes

### F1 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
Codex created a new file `mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl` (148 lines; the directive nominated `scripts/...wl` but the orchestrator relocated it into the repository's canonical `mathematica/` directory, consistent with the sibling layout for other stages). The diff (`exec_logs/stage_007_diff.patch`) shows exactly one file added, no collateral edits to any other script or the sympy file.

The new script:

- Declares native Mathematica `$Assumptions` with positivity for `lambda, mu0, sigma, tau, epsilon` and reality + nonzero for `etaSym` (lines 9-11).
- Defines `Z[w_] := Exp[-w^2/lambda^2]` and uses `Integrate[..., {w, -Infinity, Infinity}]` directly to obtain `Zint`, `Z2int` (lines 13-16).
- Each of M1-M11 is asserted by computing `FullSimplify[residual]`, printing it, and gating an `Exit[1]` on a separate `FullSimplify[... ] === 0` (or `=== 1` for M11) check.

Specifically:
- M1 (`Zint - Sqrt[Pi]*lambda`) at lines 19-24
- M2 (`Z2int - Sqrt[Pi/2]*lambda`) at lines 27-32
- M3 (`IWZsmooth - lambda/Sqrt[lambda^2 + sigma^2]`) at lines 41-46
- M4 (`IWSsmooth - 1/(Sqrt[Pi]*Sqrt[sigma^2 + tau^2])`) at lines 49-54
- M5 (`fieldMutationDelta - etaSym*lambda^3*sigma^2/(2*(lambda^2 + sigma^2)^(3/2))`) at lines 56-71
- M6 (`sourceMutationDelta - etaSym*mu0*x*sigma^2*tau^2/(2*Sqrt[Pi]*(sigma^2 + tau^2)^(3/2))`) at lines 73-89
- M7 (`IWZmatch - 1/Sqrt[2]`, computed by re-integrating `Integrate[Wmatch[w] Z[w], ...]` rather than restating `Z2int/Zint`) at lines 91-100
- M8 (`mu0ProjMatch/mu0Red - Sqrt[2]` with `IWSmatch = Wmatch[0]` for the delta-sampling step) at lines 102-112
- M9 (`IWZeps - lambda/Sqrt[epsilon^2 + lambda^2]`) at lines 114-124
- M10 (`IWSeps - Sqrt[2]/(2*Sqrt[Pi]*epsilon)`) at lines 127-132
- M11 (`Limit[IWZeps, epsilon -> 0, Direction -> "FromAbove", Assumptions -> lambda > 0] === 1`) at lines 134-145

**Assessment:**
The script is substantively independent, not a SymPy transliteration:

- Variable naming is camelCase Mathematica-native (`Wsmooth`, `IWZsmooth`, `Fmut`, `Wmatch`, `IWZmatch`, `mu0ProjMatch`, `Weps`) versus the SymPy script's snake_case (`W_smooth`, `I_WZ_smooth`, `F_w_dependent`, `W_match`, `I_WZ_match`, `mu0_proj_match`, `W_eps`). The directive explicitly prescribed the camelCase names, and the script obeys.
- All Gaussian integrals are computed by Mathematica's native `Integrate[..., {w, -Infinity, Infinity}]`, not by manual moment formulas pasted in. The script never sets `Zint := Sqrt[Pi]*lambda` and then "checks" it — it integrates first, then asserts the closed form.
- M7 is computed by re-integrating `Integrate[Wmatch[w] Z[w], ...]` and asserted against `1/Sqrt[2]` directly, not by symbolically restating `Z2int/Zint` (which would be a tautology after M1, M2). This satisfies the directive's explicit non-tautology demand on this claim.
- M5 and M6 require Mathematica to evaluate the second Gaussian moment `Integrate[Wsmooth*Z*w^2, ...]` and `Integrate[Wsmooth*Ssmooth*w^2, ...]` natively; the assertion targets are the closed-form moment results, so a wrong evaluation would be caught.
- M11 uses `Limit[..., Direction -> "FromAbove"]`, an independent Mathematica operation, rather than substituting `epsilon -> 0` after simplification.
- M8 uses the same hand-installed `Wmatch[0]` delta-sampling convention as the SymPy script (line 124 there, line 102 here). This is a modeling choice for the delta-source case shared by both engines and was explicitly prescribed by the directive (item 6), so it does not violate independence.

All eleven M-claims produce residual `0` in the exec log; the script terminates with `STATUS: PASS` and exit code 0.

The directive's nominal output-path requirement (a `.txt` in `scripts/output/`) is satisfied at the parallel path `mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt` (mtime 2026-05-21 11:51, newer than the script's 11:01).

No side effects: the SymPy `.py` is untouched (diff confirms), and no other files in the repo were modified.

## Exec log assessment

**SymPy:** exit=n/a. No `stage_007_sympy.log` was generated this iteration because the directive did not modify the SymPy script and the orchestrator only re-ran the engine whose script changed. The pre-existing SymPy output (`scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt`, mtime 2026-05-21 11:26) is unchanged and still ends in `STATUS: PASS`.

**Mathematica:** exit=0. Notable lines from `exec_logs/stage_007_mathematica.log`:

```
M1 residual = 0
PASS: M1 Gaussian profile area
...
M7 residual = 0
PASS: M7 matched-observer overlap
M8 residual = 0
PASS: M8 matched projection/reduction ratio
...
M11 residual = 0
PASS: M11 sharp-sampling limit
STATUS: PASS
# exit_code: 0
```

Every M1-M11 residual prints as `0` and the corresponding `PASS:` line follows. No warnings, no `FAIL:` lines, no `Exit[1]` triggers.

**Output freshness:**
- `mathematica/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.wl` mtime: 2026-05-21 11:01
- `mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt` mtime: 2026-05-21 11:51 (newer than script — fresh)
- `redteam/exec_logs/stage_007_mathematica.log` mtime: 2026-05-21 11:16 (newer than script; slightly older than the saved `.txt`, indicating the engine was re-invoked one extra time after the log was captured — content-identical, both PASS)
- `scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt` mtime: 2026-05-21 11:26 (newer than the unchanged sympy script at 2026-05-04)

All outputs are fresh relative to their scripts.

## Material-change assessment

`material_change`: false.

The fix only adds a new second-engine verification. The SymPy-derived results (the `sqrt(2)` projection-vs-reduction ratio, the smooth-overlap closed forms, the regulator divergence) are now independently confirmed by Mathematica, not changed. Downstream units that depend on stage-007 conclusions see exactly the same numerical/symbolic content; this change strengthens confidence rather than altering any derived quantity.

## Side observations (non-blocking)

- The new file lives in `mathematica/` rather than the `scripts/` path nominated by the directive. The orchestrator's exec engine clearly handles the canonical `mathematica/` layout (the log records `math -script /.../mathematica/.../*.wl` and the run succeeds). The directive's path string is a stale stylistic artifact — verification does not block on it. If desired, future directives could specify `mathematica/` directly to match the repository layout.
- M11's `Limit[..., Assumptions -> lambda > 0]` redeclares the positivity assumption locally even though `$Assumptions` already carries it. Harmless redundancy; possibly defensive against a Mathematica idiosyncrasy where `Limit` does not always honor `$Assumptions`.
- The script's `M7` check is implemented exactly as the directive asked (`FullSimplify[IWZmatch - 1/Sqrt[2]] === 0`) and crucially does not short-circuit through `Z2int/Zint`. This was the auditor's stated worry-point about tautology, and the implementation cleanly avoids it.
- The SymPy script's A10 (`I_WZ_match - Z2_int/Z_int == 0`, a partially-tautological restatement) has no analog in the new `.wl` — the Mathematica script went straight to the substantive `1/Sqrt[2]` form, which is arguably stronger than the original SymPy on that one assertion. Not a finding, just a note.

## Verdict justification

The single finding (missing Mathematica engine) is resolved: a new, substantively independent `.wl` script exists; it uses Mathematica-native `Integrate`/`FullSimplify`/`Limit` to derive every one of the eleven claim-manifest items M1-M11; each assertion is non-tautological (in particular M7 re-integrates rather than restating `Z2int/Zint`); the exec log shows exit code 0 with all eleven `PASS:` lines and `STATUS: PASS`; outputs are fresh; the diff shows no collateral edits; no regressions appear. The path-location deviation (`mathematica/` vs `scripts/`) is a stylistic relocation the orchestrator handled and does not affect substance.
