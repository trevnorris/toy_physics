---
unit_id: 055
batch: III.2
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T18:05:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 055

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:56` now reads
  `expect_zero("KX/KW equivalence", (1 / AK).subs(y, 0).subs(x, x_floor) - pi**2 / (4 * zeta_req))`,
  matching the directive's "After" string exactly. Surrounding lines 54-55 (the comment and the
  now-vestigial `KX_over_KW = sp.symbols(...)` placeholder) are left intact as the directive
  explicitly required.
- `mathematica/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.wl:60` now reads
  `expectZero["KX/KW equivalence", ((1/aK) /. y -> 0 /. x -> xFloor) - Pi^2/(4 zetaReq)];`.
  The directive notes the literal line number in the working tree had shifted from 54 to 60 due
  to pre-existing batch edits; the substantive content matches the required "After" string. No
  other Mathematica lines were modified.

**Assessment:**
The edit correctly re-anchors the equivalence check to the `AK` definition. The new LHS,
`1/AK` evaluated at `y = 0`, depends on the actual algebraic form of `A_K` (numerator/denominator
in `1 - x/4 + x*y^2/pi^2`); a regression in the `A_K` definition (e.g. an off-by-one in the
denominator's polynomial in `x`) would now propagate into a nonzero residual instead of being
masked by the hand-typed `1 - x/4`. The check still vanishes when `x = x_floor` because, by
construction at `y = 0`, `1/AK = 1 - x/4`, and substituting `x_floor = 4 - pi^2/zeta_req`
yields exactly `pi^2/(4*zeta_req)`. Thus the assertion is non-tautological in the auditor's
sense — it now exercises two upstream objects (A_K's form and the x_floor closed form)
simultaneously — yet it still passes on the current consistent definitions, as confirmed by
both fresh outputs reporting `KX/KW equivalence = 0` (SymPy) and `PASS: KX/KW equivalence`
(Mathematica). No collateral edits beyond the two cited lines.

## Exec log assessment

**SymPy:** exit=0 (per fix_batch_III.2.log line 85, refresh sympy succeeded at 17:41:25; no
HALT line for stage 055). Notable lines from the refreshed
`scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt`:
- `twin value = 0`
- `closure maximum = 0`
- `x floor = 4 - pi^2/zeta_req = 0`
- `KX/KW equivalence = 0`

**Mathematica:** exit=0 (per fix_batch_III.2.log line 86, refresh mathematica succeeded at
17:41:26; no HALT line for stage 055). Notable lines from the refreshed
`mathematica/output/moving_throat_pde_stage055_explicit_reachability_mathematica_audit.txt`:
- `PASS: twin value`
- `PASS: closure maximum`
- `PASS: x floor = 4 - Pi^2/zeta_req`
- `PASS: zeta_max(x_floor) - zeta_req`
- `KX/KW equivalence = 0` followed by `PASS: KX/KW equivalence`
- `Stage 055 Mathematica audit passed.`
The `Limit::alimv` warnings about assumptions involving the limit variable are pre-existing,
the same as those flagged in the original auditor report, and do not affect the residual.

**Output freshness:** confirmed. The SymPy script mtime is 2026-05-22 17:40 and its output
mtime is 2026-05-22 17:41 (output ~1 min newer than script). The Mathematica script mtime is
2026-05-22 17:40 and its output mtime is 2026-05-22 17:41. Both transcripts are post-fix.

## Material-change assessment

`material_change`: false.

The edit only strengthens an assertion; it does not change any derived numerical or symbolic
result that a downstream stage could consume. `Omega_exp`, `A_K`, `zeta_twin`, `zeta_max`,
and `x_floor` are unchanged, and the printed regime-split narrative is identical. No
downstream unit pulls from this file's outputs in a way that would be perturbed.

## Side observations (non-blocking)

- The `KX_over_KW = sp.symbols("KX_over_KW", positive=True, real=True)` line at script line 55
  is now an unused placeholder (the new assertion does not reference it). The directive
  explicitly forbade modifying that line, so this is correct verifier behavior, but the user
  may wish to clean it up in a future pass.
- The orchestrator's preemptive `ConditionalExpression -> e` strip in `expectZero` (lines
  22-27 of the `.wl`) is unrelated to F1 and does not affect the new assertion, since
  `FullSimplify` already collapses the residual to a bare `0` here. Noted only because the
  orchestrator flagged it in the prompt.
- The Mathematica banner still reads `STAGE 038`; the auditor did not flag this and it is
  outside F1's scope, so it does not affect verification.

## Verdict justification

The single low-severity finding (F1, tautological_check) was applied exactly as specified in
both engines. The SymPy line matches the directive's "After" string verbatim at line 56;
the Mathematica change lands at line 60 instead of line 54, but the directive's `## Applied:`
block flagged this drift up front and the substantive content is identical to the prescribed
"After" string. The new assertion is genuinely non-tautological (it depends on the algebraic
form of `A_K`, not just on `x_floor`'s closed form), and both freshly refreshed transcripts
show the assertion still evaluates to zero with passing exit. No collateral edits, no
regressions, no downstream impact.
