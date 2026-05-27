---
unit_id: 075
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 075 (v2)

This supersedes the v1 verifier report for unit 075. The v1 report classified the v1
patches as "resolved", but the v2 audit re-opened the findings on the grounds that the v1
"identity check" and "round-trip check" were both structurally tautological (i.e., they
both reduce to algebraic-fraction cancellations regardless of physics). The v2 directive
addresses those tautologies via genuine asymptotic-limit checks; this report verifies the
v2 fixes.

## Per-finding outcomes

### F1 — tautological_check (re-opened by v2 audit, fixed by v2 directive)

**Classification:** resolved

**What changed:**
- SymPy: `scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:68-82`
  inserted a new asymptotic-limit block immediately after the (pre-existing, now
  informational) free-symbol identity asserts. The new code computes:
  - `large_alpha_check = sp.limit(alpha_sym * Deltainf_sym, alpha_sym, sp.oo)` followed by
    `assert large_alpha_check == 1`
  - `small_alpha_check_delta0 = sp.limit(Delta0_sym, alpha_sym, 0)` followed by
    `assert small_alpha_check_delta0 == sp.Rational(1, 2)`
  Both limits are taken on the free-symbol closed forms `Delta0_sym` / `Deltainf_sym`
  (declared `positive=True, real=True`), not on the numerically-substituted expressions.
- Mathematica: `mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:96-114`
  inserted the parallel block as a `Module[{aSym, eSym, ...}, ...]` that re-declares the
  free-symbol forms locally, calls `Limit[aSym*deltaInfSym, aSym -> Infinity, Assumptions ->
  eSym > 0]` and `Limit[delta0Sym, aSym -> 0, Assumptions -> eSym > 0]`, and gates each on
  `TrueQ[... === 1]` / `TrueQ[... === 1/2]` with `pass[...]` / `fail[...]`.

**Assessment:**
The asymptotic limits are non-trivial consequences of the closed forms — a wrong constant
in the numerator or denominator of `Delta_inf` would change `lim alpha*Delta_inf` away
from 1, and a wrong factor in `Delta_0` would change the small-alpha limit away from 1/2.
SymPy's `sp.limit` and Mathematica's `Limit` are independent algorithms (different
internal series/asymptotic-handling code paths), so the two engines now exercise a
genuinely independent quantity rather than two copies of a CAS-level fraction
cancellation. Both saved transcripts confirm the limits print the expected values and
both engines PASS:
- SymPy output lines 8-10:
  ```
  alpha * Delta_inf large-alpha limit = 1
  Delta_0 small-alpha limit = 1/2
  PASS: alpha * Delta_inf -> 1 (large-alpha) and Delta_0 -> 1/2 (small-alpha)
  ```
- Mathematica output lines 25-28:
  ```
  alpha * Delta_inf large-alpha limit = 1
  PASS: alpha * Delta_inf -> 1 (large alpha)
  Delta_0 small-alpha limit = 1/2
  PASS: Delta_0 -> 1/2 (small alpha)
  ```
The pre-existing v1 tautological identity asserts (A1/A2, M1/M2) and the v1 round-trip
asserts (A3/A4, M3/M4) were retained as informational per the directive's "may remain in
place as informational consistency checks" allowance; they are now correctly labeled as
tautological in the Mathematica comment (see F3 below). No collateral edits beyond what
the directive specified.

### F2 — paper_misalignment

**Classification:** resolved (user direction (a), 2026-05-27)

**What changed:**
Per the directive's `## Approved by user` block, the orchestrator applied paper-side and
notes-side edits replacing the stale `117 Theta_w` (paper Inputs and body) and `168
Theta_w` (notes section 3 and section 4 arithmetic) text with `100 Theta_w`, matching the
script's actual `alpha_r = 10 -> alpha_r^2 = 100` and the paper's own boxed `Theta_fail` /
`Theta_suff` numerics. No script change was needed for the value itself; the script-side
companion is F4. (The verifier does not read paper.tex / notes/, but the directive records
the destination edits, and the script-side `alpha_r^2 = 100` is independently confirmed by
the F4 lock below.)

**Assessment:**
Direction (a) is internally consistent: the script's `alpha_r = 10` ⇒ `alpha_r^2 = 100`
reproduces the paper's boxed `Theta_fail / Pe_req = 3.62605617972939e-4` and `Theta_suff /
Pe_req = 4.2149534156997728721e-2` exactly (saved SymPy output lines 32-33; saved
Mathematica numeric checks lines 45-48 report `diff = 0`). The two stale text values (117
and 168) did not match those boxed numbers, so removing them was the only direction that
preserved every load-bearing number. F2 is closed by user resolution; the script-side
drift-prevention work is F4.

### F3 — mathematica_transliteration (subsumed by F1 per directive)

**Classification:** resolved

**What changed:**
- Mathematica `..._audit.wl:79-83`: the misleading v1 comment block ("This identity check
  is the independent-derivation leg required by the second-engine policy: Mathematica's
  FullSimplify must prove the identity for symbolic alpha and eta...") was replaced with
  an accurate comment:
  ```
  (* Note: the algebraic identity below is structurally tautological (it follows
     from the definition of Delta_0 / Delta_inf by canceling a common factor).
     The genuine independent check is the asymptotic-limit block further below,
     which exercises a non-trivial property of the closed forms via Mathematica's
     Limit operator (computed independently from SymPy's sp.limit). *)
  ```
- The new asymptotic-limit block from F1 sits immediately below the tautological-identity
  block and provides the actual independent check; the SymPy script obtains the same
  property via `sp.limit`.

**Assessment:**
F3 was structurally subsumed by F1 per the directive. With F1 in place, the Mathematica
script now performs at least one substantive check (the two asymptotic limits) computed by
an algorithm independent of SymPy's `sp.limit`. The rewritten comment correctly describes
the situation and no longer claims the tautological identity is the independent leg.

### F4 — script-side `alpha_r^2 == 100` lock (new finding per user direction (a))

**Classification:** resolved

**What changed:**
- SymPy `..._sympy_audit.py:24-29`: immediately after `alpha_r = sp.Integer(10)`, the
  directive's lock was inserted:
  ```python
  assert alpha_r**2 == 100, "paper Inputs line lock: alpha_r^2 must equal 100"
  print("alpha_r^2 (paper Inputs line lock) =", alpha_r**2, "  PASS")
  ```
- Mathematica `..._audit.wl:37-41`: immediately after `alphaR = 10;`, the directive's lock
  was inserted:
  ```
  expectZero["alpha_r^2 - 100 (paper Inputs line lock)", alphaR^2 - 100];
  ```
- Both engines also got their banners relabeled `STAGE 58 / 058` -> `STAGE 075` per the
  rename portion of the orchestrator's work (clears the v1 side-observation about the
  stale banner).

**Assessment:**
The lock surfaces any future drift between the paper Inputs line (`Upsilon_w = 100
Theta_w`) and the script's `alpha_r` integer as an explicit assertion failure. The
assertion is non-tautological in the relevant sense: it constrains the literal integer the
script uses, so a future edit that changed `alpha_r` to e.g. `sp.sqrt(117)` or `11` would
trip the assert. SymPy output line 5:
```
alpha_r^2 (paper Inputs line lock) = 100   PASS
```
Mathematica output lines 5-6:
```
alpha_r^2 - 100 (paper Inputs line lock) = 0
PASS: alpha_r^2 - 100 (paper Inputs line lock)
```
Both lock checks fire and PASS. No deviation from directive.

## Exec log assessment

**SymPy:** exit=0 (per orchestrator-confirmed run; the per-finding `redteam/exec_logs/stage_075_sympy.log`
file is absent on disk because the orchestrator applied F1/F3/F4 directly rather than via
a Codex iteration). The refreshed transcript at
`scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt`
(mtime 2026-05-27 02:15, fresher than the script mtime 02:12) shows the post-fix run.
Notable lines:
- Line 5: `alpha_r^2 (paper Inputs line lock) = 100   PASS` — F4 lock fires.
- Lines 6-7: `Delta_0 algebraic identity (free alpha, eta) = 0` / `Delta_inf algebraic
  identity (free alpha, eta) = 0` — v1 informational asserts still pass.
- Lines 8-10: `alpha * Delta_inf large-alpha limit = 1` / `Delta_0 small-alpha limit =
  1/2` / `PASS: alpha * Delta_inf -> 1 (large-alpha) and Delta_0 -> 1/2 (small-alpha)` —
  F1 substantive asymptotic checks pass.
- Lines 35-36: `Upsilon_fail - alpha_r^2 * Theta_fail = 0` / `Upsilon_suff - alpha_r^2 *
  Theta_suff = 0` — v1 round-trip asserts (informational) still pass.

The Python interpreter reached the final ledger print without raising `AssertionError`,
which means all five `assert` statements (alpha_r^2 == 100, delta0_identity == 0,
deltainf_identity == 0, large_alpha_check == 1, small_alpha_check_delta0 ==
sp.Rational(1,2), plus the two round-trip asserts) passed; the script's exit was 0.

**Mathematica:** exit=0 (`Exit[0];` at script line 131; transcript ends with `Stage 075
Mathematica audit passed.`). The `redteam/exec_logs/stage_075_mathematica.log` file is
also absent on disk; the refreshed transcript at
`mathematica/output/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.txt`
(mtime 2026-05-27 02:19, fresher than the script mtime 02:12) shows 15 PASS lines:
- Lines 5-6: F4 lock PASS.
- Lines 21-22 / 23-24: v1 algebraic identity PASS lines (informational).
- Lines 25-26 / 27-28: F1 asymptotic-limit PASS lines (substantive — large-alpha → 1,
  small-alpha → 1/2).
- Lines 29-30 / 31-32: v1 round-trip PASS lines (informational, tautological).
- Lines 33-48: eight `expectApprox` numeric checks all report `diff = 0` and PASS.
- Line 50: `Stage 075 Mathematica audit passed.`

**Output freshness:** confirmed.
- SymPy: script mtime 2026-05-27 02:12, output mtime 02:15 (output 3 min fresher).
- Mathematica: script mtime 2026-05-27 02:12, output mtime 02:19 (output 7 min fresher).
Both `.txt` outputs were re-generated post-fix.

## Material-change assessment

`material_change`: false.

The four changes are: (i) two non-tautological asymptotic-limit assertions on the
free-symbol closed forms, (ii) a comment-only rewrite, (iii) a defensive script-side
constant lock, and (iv) a banner relabel. None of these touch the numeric values
`Delta_0`, `Delta_inf`, `Upsilon_fail`, `Upsilon_suff`, `Xi_fail`, `Xi_suff`, `Theta_fail`,
`Theta_suff` that downstream Stages 076-078 consume. Comparing the v1-era values quoted in
the v2 audit report (`Delta_0 = 0.00017330207902152514906`, `Delta_inf =
0.020144756554052159427`, `Upsilon_fail/Pe_req = 0.036260561797293886969`,
`Theta_fail/Pe_req = 0.00036260561797293886969`, `Theta_suff/Pe_req =
0.042149534156997728721`) against the v2 post-fix transcript (SymPy lines 17-18, 28-33;
Mathematica lines 13-14, 19-20): byte-for-byte identical. Downstream is not affected.

## Side observations (non-blocking)

- The directive lacks `## Applied: F1 / F3 / F4` blocks because the orchestrator applied
  the edits directly rather than handing off to Codex. The destination state on disk is
  correct; this is a process note, not a verification failure.
- Per-stage exec log files `redteam/exec_logs/stage_075_sympy.log` and
  `..._mathematica.log` are absent on disk. The `.txt` transcripts serve as the post-fix
  evidence and confirm exit-0 + PASS for every assertion. If a persistent exec-log audit
  trail is desired, the orchestrator may want to copy or re-emit them on the next batch.
- The v1 verifier report previously living at this path classified the v1 patches as
  resolved with `findings_total: 3`. The v2 audit reopened all three (tautology shifted
  rather than discharged) and added F4. This v2 verifier report supersedes the v1
  classification.
- The pre-existing tautological identity checks (A1/A2, M1/M2) and round-trip checks
  (A3/A4, M3/M4) were intentionally retained per the directive ("may remain in place as
  informational consistency checks — they are tautological but harmless, and removing them
  is out of scope"). The Mathematica comment block now correctly labels them as
  informational and points at the asymptotic block as the substantive check. The SymPy
  comment around lines 68-82 likewise distinguishes the new substantive check from the
  prior tautological ones.

## Verdict justification

All four findings in the v2 directive are `resolved`. F1 is closed by the new
asymptotic-limit asserts in both engines, computed on free symbols via independent CAS
algorithms (`sp.limit` vs. `Limit[]`) — these are non-trivial consequences of the closed
forms that would fail if a factor were wrong, unlike the v1 identity / round-trip checks
that were algebraic-fraction cancellations. F3 is closed by the F1 fix (which provides a
genuine independent check) plus the rewritten Mathematica comment that no longer mislabels
the tautological identity as the independent leg. F2 is closed by user direction (a) on
the paper/notes side, with F4 added as the script-side companion that locks `alpha_r^2 =
100` against future drift. Both engines exit 0; the Mathematica transcript shows 15 PASS
lines; both `.txt` outputs were re-generated post-fix and are fresher than their scripts;
no downstream-affecting numeric changed. No regressions in the diff and no new findings
(verifier scope). Verdict: `verified`.
