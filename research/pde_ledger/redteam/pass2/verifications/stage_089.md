---
unit_id: 089
batch: III.5
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-05T18:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 2
findings_total: 2
material_change: false
---

# Verification — unit 089

**CHECKPOINT stage** — verified against the HIGHER BAR (both engines required,
substantive non-tautological assertions, exact paper alignment). It clears the
higher bar after the fix (see Verdict justification).

## Per-finding outcomes

### F1 — tautological_check (de-tautologize, not delete)

**Classification:** resolved

**What changed:**
- SymPy `scripts/moving_throat_pde_stage089_family1_minimal_isotropic_verdict_sympy_audit.py:67-70`:
  `eps_blk = sp.Integer(0)` → `eps_blk = sp.symbols("eps_blk", positive=True, real=True)`;
  introduced general `Q_gen = (1 + (1 - 2*eps_blk)*zeta)/(1 - eps_blk*zeta)`;
  `Q = sp.simplify(Q_gen.subs(eps_blk, 0))` (still the unblocked `1+zeta`).
- SymPy `:78`: the X−X assert `expect_zero(..., Q - (1+zeta), ...)` → `expect_zero("Stage-082 Q(zeta;0)=1+zeta reduction", Q_gen.subs(eps_blk, 0) - (1 + zeta), tol=1e-30)`.
- Mathematica `.wl:73-75`: after the general `q[zeta_,eps_]` definition, added
  `Clear[zetaRed]; expectZero["Q(zeta;0)=1+zeta reduction", q[zetaRed, 0] - (1 + zetaRed)]`.
- Mathematica `.wl:32-36`: added a minimal `expectZero` helper (sanctioned, recorded deviation; the file previously had only `expectTrue`/`expectApprox`).

**Assessment:**
Correct and non-tautological. The new check evaluates the GENERAL form
(`Q_gen` / `q[zetaRed, eps]`) at `eps=0` and compares to the claimed `1+zeta`,
so it exercises Q's algebraic structure rather than copying itself. A
transcription error in the general numerator/denominator (e.g. a wrong `zeta`
coefficient that does not vanish at `eps=0`) would make the difference nonzero
and FAIL. SymPy output line 8 prints `Stage-082 Q(zeta;0)=1+zeta reduction = 0`;
Mathematica output line 10 prints `Q(zeta;0)=1+zeta reduction = 0` + `PASS`.
The `expectZero` helper is mathematically correct:
`res = FullSimplify[Together[Expand[expr]]]`, pass iff `TrueQ[res === 0]` — an
exact symbolic-zero test (stricter than the numeric `expectApprox`), with the
common-denominator `Together` guarding against unsimplified rational forms.
Downstream `Q` is unchanged (`= 1 + zeta`), so `rho_suff/fail/max` are
byte-identical (confirmed against output). No collateral edits beyond F1/F2.

### F2 — tautological_check (de-tautologize)

**Classification:** resolved

**What changed:**
- SymPy `:125-135`: removed `Pe_req = sp.Integer(0); expect_zero("Pe_req ...", Pe_req)`
  (`0==0`). Replaced with `zero_bias_margin = sp.N(zeta_F1_at_zero - zeta_min, 40)`,
  a print, and `if not (zero_bias_margin > 0): raise AssertionError(...)`;
  then `Pe_req = sp.Integer(0)` constructed as the consequence.
- Mathematica `.wl:116-123`: removed `peReq = 0; expectApprox[..., peReq, 0, 10^-30]`.
  Replaced with `zeroBiasMargin = N[zetaF1AtZero - zetaMin, 40]`, a print, and
  `If[TrueQ[zeroBiasMargin > 0], pass[...], fail[...]]`; then `peReq = 0`.

**Assessment:**
Correct and genuinely can-fail. The self-check (`0==0` / `expectApprox[0,0]`) is
GONE in both engines. The replacement asserts positivity of the zero-bias
success margin `zeta_F1(0) - zeta_min` (≈ 0.6667), which is exactly the
precondition forcing `Pe_req = 0`: if `zeta_F1(0)` ever fell to/below
`zeta_min = 1/3` the margin would be ≤ 0 and the SymPy branch would raise /
the Mathematica branch would `fail` (Exit[1]). Note Codex used the directive's
explicit positivity-raise form, NOT the auditor's originally-proposed
`sp.Piecewise(...nan...)` gate — this correctly avoids the silent-pass trap
(`abs(complex(sp.nan)) > tol` evaluates False). The printed `Pe_req = 0` Output
line remains: SymPy output line 32 (`Pe_req = 0 -> paper Output app-stage089-Pe-zero verified`),
Mathematica output line 40 (PASS for the margin). SymPy output line 26 and
Mathematica output line 39 show the new margin = 0.6667185954688619953260007...
in both engines (agree).

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `Stage-082 Q(zeta;0)=1+zeta reduction = 0` (F1, general-form reduction)
- `zero-bias success margin zeta_F1(0) - zeta_min = 0.6667185954688619953260007503806538478887` (F2, can-fail)
- `Pe_req = 0     -> paper Output app-stage089-Pe-zero verified` (Output preserved)

**Mathematica:** exit=0. Notable lines:
- `Q(zeta;0)=1+zeta reduction = 0` + `PASS: Q(zeta;0)=1+zeta reduction` (F1)
- `zero-bias success margin zeta_F1(0) - zeta_min = 0.6667185954688619953260007503806538478887205567344657025096`40.` + `PASS: zero-bias success margin positive (=> Pe_req = 0)` (F2)
- `Stage 089 Mathematica audit passed.`

Both engines agree on every emitted deliverable (A_F1 = 1.00005192880219532865933408…,
zeta_max = 2.46752922945601223332958450…, rho_suff/fail/max, all Delta margins,
zero-bias margin 0.66671859546886199532600075…).

**Output freshness:** confirmed regenerated post-fix.
`.py` mtime 18:06:51, `.wl` mtime 18:07:41; both committed `.txt` outputs mtime
18:15:21 (newer). Committed `.txt` content is byte-consistent with the captured
exec logs (sympy 18:12:28, mathematica 18:13:41), both ending exit_code 0.

## Material-change assessment

`material_change`: **false**.

Every deliverable value is unchanged: `Q` is still
`sp.simplify(Q_gen.subs(eps_blk,0)) = 1 + zeta`, so `rho_suff = 3.46622291347846…`,
`rho_fail = 3.46752913273870…`, `rho_max = 3.46752922945601…`,
`zeta_min = 1/3`, `zeta_max = 2.46752922945601…`, `A_F1 = 1.00005192880220…`,
all `Delta_*` margins, and `Pe_req = 0` all match the audit-report reconciliation
table exactly. F1 only de-tautologizes a structural identity check (introduces no
new constant); F2 replaces a `0==0` echo with a positivity assertion on an
already-printed auxiliary margin (`Delta_AF1`, internal scaffolding) and keeps
`Pe_req = 0`. No downstream unit can be affected.

## Side observations (non-blocking)

- The `zero-bias success margin` and the previously-printed `Delta_AF1` are
  numerically the same quantity (`A_F1 - zeta_min`, since `zeta_F1(0) = A_F1`);
  the new F2 line is a deliberate restatement that ties the margin to the
  Pe_req=0 consequence with a can-fail guard. Not a defect.
- All "Stage 089" self-labels remain correct; the "Stage-075/082/62" references
  are cross-references (deferred per numbering policy). No stale self-label was
  reintroduced by the edits.
- `git status` shows only the four stage089 files touched for this unit (plus
  separate batch-III.5 sibling units 087/088); no paper/notes/prose edits.

## Verdict justification

Both findings are `resolved`, both engines exit 0 with all checks PASS, and the
diff contains no regressions. As a CHECKPOINT this clears the higher bar: the two
former tautologies are now substantive can-fail assertions present in BOTH engines
(F1 a general-form `Q(zeta;0)=1+zeta` structural reduction, F2 a positivity guard
on the zero-bias success margin that forces `Pe_req=0`); the engines agree to
displayed precision on every deliverable; the added `expectZero` helper is correct
(exact symbolic-zero, sanctioned deviation); paper alignment is preserved (the
boxed `Pe_req = 0` Output line still prints, no new constant introduced); and no
remaining tautological named assertion exists (A4 and A10/B9 from the audit
inventory are both eliminated). material_change is false.
