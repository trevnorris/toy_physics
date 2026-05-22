---
unit_id: 024
batch: II.1
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 024

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Codex restructured the angular-algebra setup in
`/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
to remove the line-by-line port of the SymPy script. Concretely (post-edit line
numbers from the current file):

1. Lines 31-43: The `pairings[list_List]` recursive helper and the hand-decomposed
   Wick-pair `i4`/`i6` definitions were replaced with native Cartesian-component
   surface integrals over the 2-sphere. `n[1] = Sin[theta] Cos[phi];
   n[2] = Sin[theta] Sin[phi]; n[3] = Cos[theta];` and `i4`/`i6` are now
   `Integrate[n[i] n[j] ... Sin[theta], {theta, 0, Pi}, {phi, 0, 2 Pi}]`.
2. Lines 50-53 (`tripleOverlap`): the inner index variable previously named `n`
   was renamed to `nn` to avoid shadowing the new `n[_]` Cartesian function.
3. Lines 117-119: Only the three physical group amplitudes remain
   (`gU = lambdaU*iEtaU; gW = lambdaW*iEtaW; rPair = lambdaR*iUW;`). The five
   SymPy shorthand symbols `deltaPair, sPair, qPair, hPair, pPair` were deleted.
4. Lines 121-133 (`zResp`, `nResp`) and lines 142-147 (the six `expectZero`
   targets for Z0/Z2/Z4/N0/N2/N4): all rewritten inline in
   `omegaU, omegaW, rPair, gU, gW` only. No shorthand identifiers appear.
5. Lines 160-165: the three pre-substituted lane forms were replaced with a
   single `xLane[lam_] := x0 + eps*lam*x1;` constructor called with
   `xLane[1]`, `xLane[1/2]`, `xLane[-1]`.

A `grep` for any of `deltaPair|sPair|qPair|hPair|pPair|pairings` in the file
returns zero hits, confirming complete removal of the SymPy-shorthand layer
and the Wick-pair helper.

**Assessment:**
The edit fully and non-tautologically addresses the finding. The
F1 verification checklist asked for: (a) `Integrate[..., {theta, 0, Pi},
{phi, 0, 2 Pi}]` calls in i4/i6 — present at lines 35-43; (b) no
`deltaPair, sPair, qPair, hPair, pPair` bindings — confirmed absent;
(c) the recursive `pairings` helper removed — confirmed absent;
(d) a single `xLane[lam_]` parameterizer called with three λ values —
present at lines 162-165.

The Mathematica engine now constructs the angular moments from a structurally
independent route: native symbolic spherical integration rather than a Python-
style Wick-pair recursion. The `Z*`/`N*` residuals reference no engine-shared
shorthand — they are written directly in `omegaU, omegaW, rPair, gU, gW` and
canonicalized by `FullSimplify[Together[Expand[...]]]` (the `expectZero`
helper). Because the assertions check residuals = 0 of fully expanded rational
functions, passing them requires the angular-integration layer and the
coupled-pair series to land on the same rational form via two structurally
different intermediate algebras. This is the substantive independence the
auditor asked for.

The 28-PASS-line directive language was an undercount; the saved output
contains 34 `PASS:` lines, matching the 34 `expectZero` invocations in the
edited `.wl` (35 total occurrences including the function definition). The
assertion names and residual content match the pre-edit list exactly — no
collateral renames or reorderings.

No side effects: the banner strings, section ordering (I → II → III → IV → V),
`$Assumptions` blocks, `Clear[...]` calls, final ledger Print block, and
`Exit[0]` line are all preserved. No new features were introduced.

## Exec log assessment

**SymPy:** exit=n/a. No `stage_024_sympy.log` is present in
`redteam/exec_logs/`; the finding did not require any SymPy script change, and
the saved SymPy output
(`scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt`,
mtime May 21 16:03) still shows all assertion residuals printing `= 0`
unchanged.

**Mathematica:** exit=n/a. No `stage_024_mathematica.log` was captured by
the orchestrator. However, the saved output
(`mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt`,
mtime May 21 16:53) was regenerated *after* the `.wl` edit (mtime May 21 16:01)
and contains 34 `PASS:` lines with all residuals reported as `0` (matrix-valued
ones as all-zero arrays). Notable lines:

- `Gram - I5 = {{0, 0, 0, 0, 0}, ...} ... PASS: Gram - I5`
- `Z0 formula = 0 ... PASS: Z0 formula` through `N4 formula = 0 ... PASS: N4 formula`
  (the six residuals that exercise the inline-expanded coupled-pair series)
- `M - M_target = {{0, 0, 0, 0, 0}, ...} ... PASS: M - M_target`
  (exercises the new sphere-integral `i6` via `tripleOverlap`)
- `b_x - 3 a_x = 0 ... PASS: b_x - 3 a_x` (exercises the new `xLane`
  constructor)
- Final ledger Print block emitted; no `FAIL:` lines.

Because the `expectZero` helper calls `Exit[1]` on any non-zero residual and
the saved output ends with the final ledger Print, the run exited 0.

**Output freshness:** confirmed. `.wl` mtime is May 21 16:01; output `.txt`
mtime is May 21 16:53 — output is post-fix. The SymPy script (unchanged) at
May 11 11:58 is older than its output at May 21 16:03, also consistent.

## Material-change assessment

`material_change`: false.

The edit is purely a setup-layer restructuring of the Mathematica audit. The
assertion list, residual mathematical content, and saved PASS lines are
identical to the pre-edit state. No derived numerical constant, closed-form
expression, or symbolic result is altered downstream. Downstream units are
unaffected.

## Side observations (non-blocking)

- The orchestrator did not capture `stage_024_sympy.log` or
  `stage_024_mathematica.log` in `redteam/exec_logs/`. The directive
  explicitly forbade Codex from running scripts, so the absence is expected,
  but if the harness intends those logs to be part of the verifier inputs
  per `verify_prompt_024.md`, the runner that produced the saved `.txt`
  outputs is not currently tee'ing into the `exec_logs/` directory. Freshness
  was confirmed via output mtimes instead.
- The `expectZero` helper applies `FullSimplify[Together[Expand[expr]],
  Assumptions -> $Assumptions]` before the `=== 0` test, which is the
  appropriate canonicalizer for the inline-expanded `Z*`/`N*` residuals (the
  directive flagged that `Together` may be needed; the helper already
  composes `Together` with `Expand`/`FullSimplify` so no per-call wrapper
  was required).
- The directive's assertion-count of "28 PASS lines" is an undercount; the
  pre-edit and post-edit output both contain 34. This is not a defect —
  Codex preserved the full assertion list — but the directive's stated
  count is inaccurate and the orchestrator should not key off it.

## Verdict justification

Codex applied the full five-part edit non-tautologically: the Wick-pair
helper and SymPy-named shorthand symbols are gone, replaced by native
sphere integrals and inline rational expressions; the lane forms are
factored through a single `xLane[lam_]` constructor; and the
assertion list, names, and residual content are preserved verbatim. The
saved Mathematica output is fresh (post-edit mtime) and all 34 assertions
PASS. The two engines now reach the same residuals = 0 conclusion via
structurally independent intermediate algebras, which is the substance
the auditor requested.
