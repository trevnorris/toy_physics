---
unit_id: 125
batch: IV.3
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-29T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 1
findings_total: 1
material_change: false
---

# Verification — unit 125

## Per-finding outcomes

### F1 — insufficient_verification (R1)

**Classification:** resolved

**What changed:**

SymPy (`scripts/moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:85-92`):
the single weak smallness assertion
`expect_true("g[peaked@L proxy a=100] < 0.05", bool(abs(g_a_large) < sp.Rational(1, 20)))`
was deleted and replaced by three assertions:
- `expect_true("g[peaked@L proxy a=100] >= 0", bool(g_a_large >= 0))` (line 88)
- `expect_true("g[peaked@L proxy a=100] <= 1", bool(g_a_large <= 1))` (line 89)
- `expect_true("g[peaked@L proxy a=100] < 1/20", bool(g_a_large < sp.Rational(1, 20)))` (line 92)

Neither range bound is wrapped in `abs(...)`. The `< 1/20` line retains the original
smallness fact verbatim (`sp.Rational(1, 20)`) but is now subordinate to the genuine
range checks, exactly as the directive's claim manifest M1/M2/M3 specified.

Mathematica (`mathematica/moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:94-97`):
the genuine endpoint identity `expectZero["moment g[peaked@L] limit", Limit[gA, aSym -> Infinity]]`
(line 93) is kept unchanged, and two mirroring range checks were inserted after it at the
same `a = 100` finite proxy:
- `gApeaked = Chop[Re[Block[{$MaxExtraPrecision = 10000}, N[(gA /. aSym -> 100), 30]]]];`
- `expectTrue["g[peaked@L proxy a=100] >= 0", gApeaked >= 0]` (line 96)
- `expectTrue["g[peaked@L proxy a=100] <= 1", gApeaked <= 1]` (line 97)

The diff (`redteam/exec_logs/stage_125_diff.patch`) confirms these are the only two edits;
no collateral changes to any other check or derived quantity.

**Assessment:**

The edit correctly addresses the finding. The old assertion accepted small NEGATIVE
moments (a value of `-0.01` would have passed `abs(...) < 1/20`), so it could pass while
violating the lower-bound half of the paper theorem `0 <= g[sigma] <= 1`. The new
unwrapped `>= 0` check closes precisely that gap.

Non-tautology / can-it-fail: the new `>= 0` check is genuinely falsifiable. `g_a_large`
is `sp.N(g_a.subs(sigma_a, 100))` where `g_a` is the SymPy-DERIVED closed form
`hyper((sigma_param/2 + 1/2,), (1/2, sigma_param/2 + 3/2), -pi**2/16)` (script line 76,
echoed at SymPy output line 25) — i.e. the assertion evaluates a concrete numeric value
of a closed form (`0.0153964175481027`, output line 27), NOT the abstract integral of a
nonneg profile against a nonneg kernel re-run at assertion time. A sign error, wrong
kernel, wrong normalization, or mis-substituted `a` in that closed form would yield a
negative numeric `g_a_large`, which `bool(g_a_large >= 0)` would FAIL. It is therefore a
real implementation guard, not a structurally-forced-nonneg tautology. The Mathematica
side is equally substantive: `gApeaked` is the numeric (real-part, chopped) value of the
closed-form incomplete-gamma expression `gA` (output line 34), so a negative numeric value
there fails `gApeaked >= 0` as well.

The bounds `0` and `1` are the literal paper bounds per the directive/finding anchors
(paper/stages/stage_125.tex:16 boxed `0 <= g[sigma] <= 1`; notes :50-55) and match the
already-present uniform-source range checks (SymPy lines 94-95, Mathematica lines 99-100).
They are not fabricated. (Verifier is scripts-only and did not open the prose; relying on
the directive/finding anchors plus the parallel uniform-source bounds already in-script.)

The recorded deviation — Mathematica using `$MaxExtraPrecision = 10000` + `N[..., 30]`
inside `Re[...]` / `Chop[...]` on the incomplete-gamma closed form — is benign: it forces
a converged 30-digit numeric value and discards a spurious imaginary cancellation
artifact so the `>= 0` / `<= 1` comparisons resolve to explicit `True`. It does not weaken
the check (a truly negative real part would survive `Re`/`Chop` and fail the bound) and
does not alter any derived value (the printed `gApeaked = 0.0153964...` matches the SymPy
value to all shown digits, output line 41 vs 27).

## Exec log assessment

**SymPy:** exit=0. Notable lines (from `scripts/output/...sympy_audit.txt`):
- `g_a at sigma_param = 100 (peaked-at-L proxy) = 0.0153964175481027`
- `g[peaked@L proxy a=100] >= 0 = True`
- `g[peaked@L proxy a=100] <= 1 = True`
- `g[peaked@L proxy a=100] < 1/20 = True`
THREE new check lines present; the old single `g[peaked@L proxy a=100] < 0.05 = True`
line is GONE. The script prints its closing Conclusion lines, so no `expect_true`
raised → exit 0.

**Mathematica:** exit=0. Notable lines (from `mathematica/output/...mathematica_audit.txt`):
- `moment g[peaked@L] limit = 0` / `PASS: moment g[peaked@L] limit` (unchanged endpoint identity)
- `g_a at aSym -> 100 (peaked-at-L proxy) = 0.0153964175481027...`
- `PASS: g[peaked@L proxy a=100] >= 0`
- `PASS: g[peaked@L proxy a=100] <= 1`
(One benign `Limit::alimv` warning about assumptions ignored on the unchanged limit line;
not a failure.) The script reaches the Conclusion block and `Exit[0]`; every `expectTrue`/
`expectZero` printed `PASS`.

PASS-line tally for the strengthened check: SymPy = 3 new lines (`>= 0`, `<= 1`, `< 1/20`,
all `True`); Mathematica = 2 new `PASS:` lines (`>= 0`, `<= 1`). Matches the directive's
expected verification output exactly.

**Output freshness:** confirmed. Both `.txt` outputs have mtime 1780090018, newer than the
SymPy script (1780086665) and the Mathematica script (1780086955) — regenerated post-fix.

## Material-change assessment

`material_change`: false.

Only the verification surface changed: the smallness-only assertion was replaced by /
augmented with genuine range checks, and two mirroring range checks were added on the
Mathematica side. No derived or carried numerical result changed — `g_a_large` /
`gApeaked` is still `0.0153964175481027...`, the branch values `g_- = 0.758035...` and
`g_+ = 2.797951...` are untouched, and the parametric moment closed form is unchanged.
No downstream unit depends on a new or altered value.

## Side observations (non-blocking)

- The `r` symbol is defined in the SymPy script (line 49) via the `sqrt(12*(37/20)^2/pi^2 - 1)`
  form and reconciled with `R/(10 pi)` at line 54; this is upstream of the finding and was
  not part of F1. No action.
- The Mathematica `Limit::alimv` warning is cosmetic (assumptions involving the limit
  variable ignored on the `a -> Infinity` limit) and the limit still resolves to `0`.

## Verdict justification

The single finding R1 is fully resolved: the `abs(...) < 0.05` smallness test is gone,
replaced on the SymPy side by unwrapped `>= 0` / `<= 1` range checks plus a subordinate
`< 1/20` trend witness, and mirrored on the Mathematica side by two `a = 100` range checks
alongside the retained exact `a -> Infinity` limit. The strengthened `>= 0` check evaluates
a concrete numeric value of a SymPy/Mathematica-derived closed form (not the abstract
integral), so it genuinely fails on a negative-moment sign error and is non-tautological.
Both engines exit 0 with the expected new PASS lines, outputs were regenerated post-fix,
and no derived result changed (material_change = false).
