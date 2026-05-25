---
unit_id: 080
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 080

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:62-83` —
Codex inserted exactly the directive's block: a fresh `_omega_explicit(p)` helper,
a `A_F1_recomputed` constructed from `Rational(12321,5)`, `pi`, and `y_F1`, and a
loop over four `(name, pe_val, expected)` tuples computing `zeta_val =
sp.N(A_F1_recomputed * _omega_explicit(pe_val)**2, 25)` and asserting
`|zeta_val - expected| < 1e-14`. Insertion is between the existing `vals` print
loop and the large-lambda saturation block, matching the directive verbatim.

**Assessment:**
The check is non-tautological: each `pe_val` is a concrete `sp.Float(c, 30)`
with `c` equal to the Stage-61-derived numerical Pe constant, evaluated through
a freshly named helper (`_omega_explicit`) and an `A_F1_recomputed` that is
algebraically reconstructed rather than reused from the existing `A_F1`. The
expected target is the printed `zeta_*(1)` value (string form, 25 digits), so
mutating any of the four Pe constants (e.g., changing `96.5285247264386` to
`0.001`) would change `_zeta_val` away from `_expected` by far more than
`1e-14` and the assertion would fail. The refreshed sympy output shows the
four diffs are `3.7e-16`, `1.6e-16`, `4.9e-16`, `3.0e-16` — all well under
`1e-14` — and the script ran to the `FINAL LEDGER` section without raising.
No collateral edits beyond the inserted block.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:55-60`
inserts the new independent-path block: `omegaIndep[p_]`, `aF1Indep` (analytic
form evaluated with `N[..., 50]`), and four `zetaTarget*` values computed from
`N[aF1Indep * omegaIndep[ToExpression["c`30"]]^2, 40]`. Lines 78-81 replace the
four literal `ToExpression["2.4xxx`25"]` targets with `zetaTargetSuffChi`,
`zetaTargetFailChi`, `zetaTargetSuffJ`, `zetaTargetFailJ` respectively.

**Assessment:**
Grep of the .wl file for the literal SymPy-output target strings
(`2.46622291347846389`, `2.46752913273870287`, `2.44257571477179128`,
`2.46752736855057761`) returns zero matches — the directive's verification
condition is satisfied. The refreshed Mathematica output reports diff=0 for
all four numeric checks and PASS lines for each. The independent path is
genuinely a different code route (`omegaIndep` defined separately, `aF1Indep`
constructed analytically rather than reusing `aF1`), so engine-agreement
within Mathematica is now substantive rather than literal-string transcription
from SymPy. No collateral edits beyond the inserted block and the four-line
replacement.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:97-100`
preserves the original `zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max`
check and adds three new `expectTrue` lines for the J-pair monotonicity,
J<=chi suff comparison, and J<=chi fail comparison, exactly matching the
directive's required-change block.

**Assessment:**
The refreshed Mathematica output prints four `PASS:` lines for these
inequalities (in addition to the four numeric checks and four limit checks).
Each `expectTrue` evaluates a concrete `<` or `<=` between independently
computed numerical zeta values, so it would catch a sign mistake in any of
the four `Pe_*` constants that swapped the ordering. The directive explicitly
scoped this fix to Mathematica only (since F1 added the substantive numerical
checks on the SymPy side), and Codex respected that scoping. No collateral
edits.

## Exec log assessment

**SymPy:** exit log not captured (file missing). Output `.txt` was regenerated
post-fix (mtime 23:38 vs script mtime 23:36) and contains the new
`zeta numeric check ...` lines plus the entire FINAL LEDGER block, which is
only printed after the `for ... raise AssertionError` block exits without
raising. So the post-fix sympy script completed successfully:

```
zeta numeric check zeta_suff^(chi)(1): diff = 3.71602824196198313509429994768E-16
zeta numeric check zeta_fail^(chi)(1): diff = 1.56302152908550079327091729582E-16
zeta numeric check zeta_suff^(J)(1):   diff = 4.90133341060895771213370823195E-16
zeta numeric check zeta_fail^(J)(1):   diff = 3.01805693496816131515201589094E-16
```

All diffs `< 1e-14`, as required.

**Mathematica:** exit log not captured. Output `.txt` was regenerated post-fix
(mtime 23:38 vs script mtime 23:36) and ends with
`Stage 080 Mathematica audit passed.`, which is only printed before
`Exit[0]` after every `expectApprox`/`expectTrue` PASSes. Notable lines:

```
PASS: zeta_suff^(chi)(1) numeric check  (diff = 0)
PASS: zeta_fail^(chi)(1) numeric check  (diff = 0)
PASS: zeta_suff^(J)(1) numeric check    (diff = 0)
PASS: zeta_fail^(J)(1) numeric check    (diff = 0)
PASS: limit zeta_*(... -> zeta_max)     (4 lines)
PASS: zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max
PASS: zeta_suff^(J)(1) < zeta_fail^(J)(1) < zeta_max
PASS: zeta_suff^(J)(1) <= zeta_suff^(chi)(1)
PASS: zeta_fail^(J)(1) <= zeta_fail^(chi)(1)
```

The only warnings are benign `Limit::alimv` notices (the limit variable
`lambdaMu` is in `$Assumptions` and Mathematica ignores those during `Limit`)
— pre-existing and unrelated to Codex's edits.

**Output freshness:** confirmed. Both `.txt` outputs at mtime 23:38 are newer
than their corresponding script mtimes at 23:36 — so the refresh happened
after Codex applied the edits.

## Material-change assessment

`material_change`: false.

The five printed `zeta_*(1)` and `zeta_max^(F1)` numerical values are
identical to those in the pre-fix output (compare report's "Engine
cross-check" block to the current `.txt` outputs). All Codex edits added
assertions (and one independent recomputation path in Mathematica); no
derivation route changed any printed symbolic content or numeric value.
Downstream units consuming the Family-1 zeta thresholds receive the same
numbers as before.

## Side observations (non-blocking)

- The `Limit::alimv` warnings in the Mathematica output are pre-existing and
  not caused by F1/F2/F3 edits; the `$Assumptions = ... lambdaMu > 0` line is
  the source. Not blocking — `Limit` still returns the correct values to ~30
  digits and all four `limit -> zeta_max` checks PASS at diff=0.
- The directive's verification command refers to running `redteam
  exec-sympy 080` / `redteam exec-mathematica 080` but the exec_logs/ for
  unit 080 contains only `stage_080_diff.patch` (no `*_sympy.log` /
  `*_mathematica.log`). The refreshed `.txt` outputs nonetheless make the
  outcome unambiguous (script-end banners are printed only on full success),
  so verdict stands.

## Verdict justification

All three findings are mechanically resolved: Codex's edits match the
directive's required changes line-for-line, the post-fix script outputs
include the new assertions printing PASS / sub-1e-14 diffs, and both engine
output banners ("FINAL LEDGER" for SymPy, "Stage 080 Mathematica audit
passed." for Mathematica) appear only on successful completion. The fixes
add assertion strength without altering any derived numerical value, so no
downstream impact.

stage 080: verified
