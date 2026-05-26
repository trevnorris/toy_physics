---
unit_id: 048
batch: III.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-26T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 048

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
In `scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py:59-66`, Codex inserted, immediately after the existing `print("limit xi->1^- of F_tr =", ...)` at line 58, exactly the prescribed block: a `soft_coeff = sp.simplify(sp.limit((1 - xi) * F_tr, xi, 1, dir="-"))` computation, a `print("softening coefficient for F_tr =", soft_coeff)` line, and an `expect_zero("(1-xi) F_tr softening coefficient", soft_coeff - (9*delta + 9 + 2*R**2)**2*(9*delta + 9 + 2*R)**2 / (81*(9*delta**2 + 18*delta + 9 + 2*R**2)**2))` assertion. The diff at `redteam/exec_logs/stage_048_diff.patch` shows exactly an 8-line insertion at line 58 of the original; no other lines in the file were touched (the existing `expect_zero("F_tr(xi=0)-1", ...)` at line 57, the limit print at line 58, the `M_crit - G_tr` block formerly at lines 59-66, the `S` pole-direction assertion, and the inverse-map / margin / `dxi/dzeta` blocks are all unchanged).

**Assessment:**
The edit matches the directive verbatim — same RHS coefficients, same Python idioms, same insertion point. No collateral edits. The new assertion is non-tautological: the LHS is a `sp.simplify(sp.limit(...))` of a non-trivial one-sided limit of `(1 - xi) * F_tr`, and the RHS is a hand-written canonical form taken from the Mathematica `softCoeffExpected`. Substituting any other algebraic value on either side would yield a nonzero residual.

Cross-engine textual agreement is confirmed: the SymPy exec log line 20 reports `softening coefficient for F_tr = (2*R + 9*delta + 9)**2*(2*R**2 + 9*delta + 9)**2/(81*(2*R**2 + 9*delta**2 + 18*delta + 9)**2)`; the Mathematica exec log line 17 reports `softening coefficient for F_tr = ((9 + 9*delta + 2*r)^2*(9 + 9*delta + 2*r^2)^2)/(81*(9*(1 + delta)^2 + 2*r^2)^2)`. After expanding `9*(1+delta)^2 = 9 + 18*delta + 9*delta^2`, the two engines' closed forms are algebraically identical, matching the prescribed RHS. The SymPy script still passes (`# exit_code: 0`, "All Stage-048 symbolic checks passed.").

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- Line 19: `limit xi->1^- of F_tr = oo` (unchanged from pre-fix)
- Line 20: `softening coefficient for F_tr = (2*R + 9*delta + 9)**2*(2*R**2 + 9*delta + 9)**2/(81*(2*R**2 + 9*delta**2 + 18*delta + 9)**2)` (new)
- Line 21: `(1-xi) F_tr softening coefficient = 0` (new — the assertion itself)
- Line 52: `All Stage-048 symbolic checks passed.`

**Mathematica:** exit=0. Notable lines:
- Line 17: `softening coefficient for F_tr = ((9 + 9*delta + 2*r)^2*(9 + 9*delta + 2*r^2)^2)/(81*(9*(1 + delta)^2 + 2*r^2)^2)`
- Line 18: `(1-xi) F_tr softening coefficient = 0`
- Line 19: `PASS: (1-xi) F_tr softening coefficient`
- Line 47: `Stage 048 Mathematica audit passed.`

Both engines now assert the strengthened softening-coefficient identity, and the residuals printed on both sides are `0`. The pre-existing checks (`dG_tr/dxi`, `F_tr(xi=0)-1`, `M_crit - G_tr`, `S(zeta=0)-1`, `dS/dzeta`, inverse maps, pole margin, branch margin, `dxi_phys/dzeta`) all continue to pass — no regression.

**Output freshness:** The orchestrator-captured exec log `redteam/exec_logs/stage_048_sympy.log` is dated `2026-05-26T01:52:47-06:00`, four minutes after the script mtime (`2026-05-26 01:48:24`), so the post-fix run is what was captured. The saved canonical `.txt` output at `scripts/output/moving_throat_pde_stage048_support_compensation_sympy_audit.txt` is still dated May 11 (i.e. older than the script) — see Side observations.

## Material-change assessment

`material_change`: false.

The added assertion is a stronger statement of an already-claimed limit (`F_tr` diverges as `xi -> 1^-`). It does not modify any derived result downstream stages consume: `zeta_req`, `zeta_crit`, `dxi_phys/dzeta`, `M_crit`, and the pole/branch margins are byte-identical to the pre-fix exec log. Stages 049-057, which consume `zeta_req` and the `zeta_phys >= zeta_req` success condition (per `\stagefield{Downstream use}`), see no change in any algebraic input. The fix only closes a verification gap on the IVT existence argument for `xi_req in (0,1)`; the existence itself is reused but the closed form is not.

## Side observations (non-blocking)

The canonical `scripts/output/moving_throat_pde_stage048_support_compensation_sympy_audit.txt` (the file under `scripts/output/`, not the orchestrator exec log) has mtime `May 11 12:48`, which predates the script mtime `May 26 01:48`. The directive correctly instructed Codex not to run python, and the orchestrator captured a fresh exec log to `redteam/exec_logs/stage_048_sympy.log` (May 26 01:52), so verification is on solid ground. If downstream tooling reads the `scripts/output/` `.txt` (rather than the exec log), it will be stale. This is non-blocking for unit 048 but may want a post-batch refresh of saved outputs.

The script and notes both still carry an internal banner label "STAGE 31" / "STAGE 031" (sympy line 8 in log; mathematica line 8 in log); the auditor already flagged this as historical numbering drift, not a math issue, and the file/print correctly identify the unit as 048. No action needed for unit 048.

## Verdict justification

The single finding F1 (low-severity `insufficient_verification`) was applied exactly as specified — the directive's required insertion appears verbatim at the prescribed location, with no deviations, no collateral edits, and no unsafe shortcuts. The new SymPy assertion is non-tautological (a non-trivial `sp.simplify(sp.limit(...))` matched against an independently-derived Mathematica canonical form), and both engines now print `0` for the residual and exit 0. The two engines' closed forms for the softening coefficient agree after trivial expansion of `9*(1+delta)^2`. No regressions: every previously-passing assertion still passes. `material_change: false` because no downstream-consumed quantity moved. Verdict: `verified`.
