---
unit_id: 080
batch: III.4
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: insufficient
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 080 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.txt`

## What the script claims to verify

The unit converts four Stage-61 Family-1 Peclet thresholds `Pe_*(lambda_mu) = c_* lambda_mu^2`
(with `c_suff_chi=96.5285...`, `c_fail_chi=11220.5...`, `c_suff_J=22.0062...`,
`c_fail_J=2558.02...`) into quadrupole-demand thresholds `zeta_*(lambda_mu)`
through the Stage-62 Family-1 map `zeta = A_F1 * Omega(Pe)^2`, with
`A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2)`, `kappa_F1 = 12321/5`,
`y_F1*tan(y_F1) = 37`, and `Omega(Pe) = pi*Pe*(2*Pe*exp(Pe)+pi) /
((4*Pe^2 + pi^2)*(exp(Pe)-1))`. The SymPy script's only assertion is that
`lim_{lambda_mu -> oo} zeta_*(lambda_mu) = zeta_max^(F1)`. The Mathematica
script reproduces the same construction and additionally checks the four
numerical values `zeta_*(1)` against literal targets and a strict
monotonicity inequality `zetaSuffChi1 < zetaFailChi1 < zetaMax`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 71-72 | `raise AssertionError if abs(lim_{lam->oo} zeta_*(lam) - zeta_max) > 1e-10` (looped over the four thresholds) | tautological — see F1 |
| A2 | mathematica | 71 | `expectApprox["zeta_suff^(chi)(1) numeric check", zetaSuffChi1, 2.4662229134784638979, 1e-14]` | no — target is the SymPy output literal (F2) |
| A3 | mathematica | 72 | `expectApprox["zeta_fail^(chi)(1) numeric check", zetaFailChi1, 2.4675291327387028754, 1e-14]` | no — target is the SymPy output literal (F2) |
| A4 | mathematica | 73 | `expectApprox["zeta_suff^(J)(1) numeric check", zetaSuffJ1, 2.4425757147717912819, 1e-14]` | no — target is the SymPy output literal (F2) |
| A5 | mathematica | 74 | `expectApprox["zeta_fail^(J)(1) numeric check", zetaFailJ1, 2.4675273685505776147, 1e-14]` | no — target is the SymPy output literal (F2) |
| A6 | mathematica | 86 | `expectApprox["limit zeta_suff^(chi) -> zeta_max", limSuffChi, zetaMax, 1e-14]` | tautological — same shape as A1 |
| A7 | mathematica | 87 | `expectApprox["limit zeta_fail^(chi) -> zeta_max", ...]` | tautological |
| A8 | mathematica | 88 | `expectApprox["limit zeta_suff^(J) -> zeta_max", ...]` | tautological |
| A9 | mathematica | 89 | `expectApprox["limit zeta_fail^(J) -> zeta_max", ...]` | tautological |
| A10 | mathematica | 90 | `expectTrue["zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max"]` | partial — verifies monotonicity of `Omega` only |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:63-72`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:76-89`

**What's wrong:**

The SymPy script defines `zeta_max = sp.simplify(sp.limit(zeta, Pe, sp.oo))`
(line 33) and then, for each threshold `zeta_*` defined by `zeta_*(lam) =
zeta.subs(Pe, c_* * lam**2)` with `c_* > 0`, asserts that
`lim_{lam -> oo} zeta_*(lam) == zeta_max` (lines 63-72). Because
`c_* * lam**2 -> +oo` as `lam -> oo` for every positive constant `c_*`,
this limit equals `lim_{Pe -> oo} zeta(Pe) = zeta_max` by construction.
The assertion holds for **any** positive constant `c_*` regardless of
its numerical value, so it does not verify the Stage-61-derived numerical
thresholds `96.5285247264386`, `11220.5441626259`, `22.0062226330754`,
`2558.01892349205` in any way. Substituting bogus values like `0.001` or
`12345` for these constants would still PASS.

The Mathematica analogues (A6-A9) have the same defect: `Limit[zetaFn[c*lam^2],
lam -> Infinity]` evaluates to `aF1 * Pi^2/4 = zetaMax` independent of `c > 0`.

**Why this matters:**

The script's docstring promises to "convert the Stage-61 Family-1 Pe_req
thresholds into explicit quadrupole-demand thresholds zeta_req using the
Stage-62 Family-1 demand map", but the only assertions it carries do not
validate either the Stage-61 thresholds (the `c_*` constants) or the
Stage-62 map's evaluation at any specific Pe value. A reader concluding
"PASS therefore the Pe thresholds map correctly into zeta thresholds" is
mistaken: the PASS only confirms the saturation behavior of `Omega` for
large argument, which is property of `Omega` alone.

**Required change:**

Add at least one non-tautological assertion per threshold that ties the
numerical zeta value at `lambda_mu = 1` back to the explicit definition
`zeta(Pe) = A_F1 * Omega(Pe)^2`, computed by an independent evaluation of
`Omega` and `A_F1` (not via the same `zeta` symbolic chain). Concretely, at
`scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py`
after line 60, insert a block that:

1. Computes `omega_val(p) = pi*p*(2*p*exp(p)+pi)/((4*p**2 + pi**2)*(exp(p)-1))`
   as a fresh sympy expression (avoid reusing the `Omega` symbol bound on
   line 31; introduce a new helper name like `omega_explicit`).
2. Computes `A_F1_recomputed = (sp.Rational(12321,5) + sp.pi**2/4) /
   (sp.Rational(12321,5) + y_F1**2)`.
3. For each `(c, expected)` in
   `[(96.5285247264386, "2.4662229134784638979"),
     (11220.5441626259, "2.4675291327387028754"),
     (22.0062226330754, "2.4425757147717912819"),
     (2558.01892349205, "2.4675273685505776147")]`
   computes `pe_val = sp.Float(c, 30)`, `zeta_val = sp.N(A_F1_recomputed *
   omega_explicit.subs(Pe, pe_val)**2, 25)`, and asserts
   `abs(zeta_val - sp.Float(expected, 25)) < sp.Float('1e-14')`.
   On failure raise an AssertionError naming the failing threshold.

This makes the four numerical zeta values load-bearing — bogus Pe constants
no longer pass.

**Verification:**

Re-running the SymPy script should print four new lines of the form
`zeta numeric check c=<c>: diff=<small>` and exit 0; mutating any of the
four Pe constants by more than ~1e-12 should cause exit 1.

### F2 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:71-74`

**What's wrong:**

The Mathematica `expectApprox` targets

```
2.4662229134784638979   (line 71)
2.4675291327387028754   (line 72)
2.4425757147717912819   (line 73)
2.4675273685505776147   (line 74)
```

are exactly the SymPy script's printed output values
(`scripts/output/.../stage080..._sympy_audit.txt` lines 18-21). They are
not derived inside the Mathematica script and are not anchored to any
external authoritative source (no Stage-62 closed-form value, no
hand-computed cross-check). The Mathematica script's "numeric check"
therefore reduces to: "Mathematica's `aF1*omegaFn[c*1^2]^2` reproduces
SymPy's `A_F1*Omega(c)**2`", which only confirms that the two engines
share identical formulas and arithmetic — it does not confirm the result.

**Why this matters:**

The four `expectApprox` calls are the only substantive assertions in the
Mathematica file (A1-A5 from the inventory, before the tautological limit
checks A6-A9 and the monotonicity check A10). If the SymPy script's
formulas were wrong, the Mathematica targets — copied from that SymPy
output — would also be wrong, and the cross-engine check would still pass.
The two engines do not currently provide independent verification.

**Required change:**

Replace the literal targets on lines 71-74 with values derived inside the
Mathematica script from an independent re-expansion of the demand map.
Concretely, after line 53 in
`mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`,
add a block that recomputes each `zeta_*(1)` through a different algebraic
path, then uses those recomputed values as the targets:

1. Define `omegaIndep[p_] := Pi*p*(2*p*Exp[p] + Pi) / ((4*p^2 + Pi^2) *
   (Exp[p] - 1));` (a freshly named helper, deliberately structurally
   identical to `omegaFn` so the two evaluation paths can be compared, but
   used through a different code path).
2. Define `aF1Indep = (12321/5 + Pi^2/4) / (12321/5 + yRoot^2);` evaluated
   to `N[..., 50]` after the `FindRoot` call (i.e., the analytic form, not
   reusing the already-numerified `aF1`).
3. For each `c` in `{96.5285247264386, 11220.5441626259, 22.0062226330754,
   2558.01892349205}`, compute `zetaTarget_c = N[aF1Indep * omegaIndep[c]^2,
   40]`.
4. Replace the literal targets in `expectApprox` (lines 71-74) with the
   corresponding `zetaTarget_c` values.

This converts the check from "Mathematica matches SymPy's printed number"
to "Mathematica's main path matches Mathematica's independent path", and
makes the SymPy and Mathematica scripts cross-checkable a posteriori by
comparing their independent outputs (rather than by one transcribing the
other).

**Verification:**

The Mathematica output should still report all four checks PASS, but the
`expectApprox` source-code targets will no longer be literal numbers
copied from `scripts/output/.../stage080..._sympy_audit.txt`. The
verifier can confirm by grepping the `.wl` file for the literal strings
`2.4662229134784638979`, `2.4675291327387028754`, `2.4425757147717912819`,
`2.4675273685505776147` and seeing zero matches.

### F3 — insufficient_verification

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:23-78`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:32-90`

**What's wrong:**

The script's "FINAL LEDGER" prose (`.py` lines 74-78) claims:

```
Guaranteed success: zeta_req <= zeta_suff^(chi)(lambda_mu)
Guaranteed failure: zeta_req >= zeta_fail^(chi)(lambda_mu)
Conservative floor success/failure obtained by replacing (chi) with (J)
```

This claim presupposes (i) `zeta_*` is a monotonically *increasing* function
of `Pe`, so that `Pe_suff < Pe_fail` implies `zeta_suff < zeta_fail`, and
(ii) the (J) thresholds are conservatively narrower than the (chi)
thresholds (i.e., `zeta_suff^(J) <= zeta_suff^(chi)` and
`zeta_fail^(J) <= zeta_fail^(chi)`). Neither is asserted in the SymPy
script. The Mathematica script asserts only `zeta_suff^(chi)(1) <
zeta_fail^(chi)(1) < zeta_max` (line 90), covering one of the four
implied ordering relations at the single point `lambda_mu = 1`.

**Why this matters:**

The "ledger" interpretation depends on these ordering facts being true.
Without an assertion, a sign mistake in any of the four `Pe_*` constants
could silently put `zeta_suff > zeta_fail` and the script would still
exit 0. The current monotonicity check (A10) catches only the chi-pair at
`lambda_mu = 1`.

**Required change:**

Extend the strict-ordering assertion at
`mathematica/.../stage080..._mathematica_audit.wl:90` from one inequality
to the four that are actually implied by the FINAL LEDGER prose. Replace
line 90 with the four checks:

```
expectTrue["zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max",
  zetaSuffChi1 < zetaFailChi1 < zetaMax];
expectTrue["zeta_suff^(J)(1) < zeta_fail^(J)(1) < zeta_max",
  zetaSuffJ1 < zetaFailJ1 < zetaMax];
expectTrue["zeta_suff^(J)(1) <= zeta_suff^(chi)(1)",
  zetaSuffJ1 <= zetaSuffChi1];
expectTrue["zeta_fail^(J)(1) <= zeta_fail^(chi)(1)",
  zetaFailJ1 <= zetaFailChi1];
```

Do not extend the SymPy script for this finding — keep the change scoped
to the Mathematica side, since F1 already adds the substantive numerical
checks on the SymPy side.

**Verification:**

Mathematica output should print four `PASS:` lines for these inequalities
(replacing the one current PASS line for the chi-only ordering). The
verifier confirms by grepping the output for the four new check names.

## Independent-derivation check (Mathematica)

The Mathematica script is a near-line-by-line transliteration of the SymPy
script's construction:

SymPy lines 28-32:
```
y_F1 = sp.nsolve(y * sp.tan(y) - 37, 1.53, tol=1e-34, maxsteps=100)
kappa_F1 = sp.Rational(12321, 5)
A_F1 = (kappa_F1 + sp.pi**2 / 4) / (kappa_F1 + y_F1**2)
Omega = ... pi * Pe * (2 * Pe * exp(Pe) + pi) / ((4 * Pe**2 + pi**2) * (exp(Pe) - 1))
zeta = A_F1 * Omega**2
```

Mathematica lines 37-42:
```
kappaF1 = 12321/5;
etaF1 = 37;
yRoot = y /. FindRoot[y*Tan[y] == etaF1, {y, 1.53}, ...];
aF1 = N[(kappaF1 + Pi^2/4)/(kappaF1 + yRoot^2), 50];
omegaFn[p_] := Pi*p*(2*p*Exp[p] + Pi)/((4*p^2 + Pi^2)*(Exp[p] - 1));
zetaFn[p_] := aF1*omegaFn[p]^2;
```

These are the same algebraic recipe in two syntaxes. There is one small
piece of independent computation: `zetaMax` is computed analytically as
`aF1*Pi^2/4` on Mathematica line 43, whereas SymPy uses `sp.limit(zeta,
Pe, sp.oo)` on line 33. That is the only place the two engines exercise
genuinely different machinery on the same claim.

Both engines also reuse the same four Stage-61 numerical Pe thresholds
(SymPy lines 36-39 ↔ Mathematica lines 45-48), so the four `zeta_*(1)`
output values are guaranteed to agree to round-off without any independent
re-derivation — and indeed the Mathematica `expectApprox` targets are the
SymPy output strings themselves (see F2).

This pattern straddles `mathematica_transliteration` and
`hardcoded_result`; F2 captures the load-bearing instance (the literal
targets). I am not raising a separate `mathematica_transliteration`
finding because the directive in F2 already routes the fix toward
introducing an independent evaluation path inside Mathematica.

## Engine cross-check

Both engines completed exit 0 with matching numerical results to ~15
significant figures:

SymPy output:
```
zeta_suff^(chi)(1) = 2.4662229134784638979
zeta_fail^(chi)(1) = 2.4675291327387028754
zeta_suff^(J)(1)   = 2.4425757147717912819
zeta_fail^(J)(1)   = 2.4675273685505776147
zeta_max^(F1)      = 2.4675292294560112350
```

Mathematica output:
```
zeta_suff^(chi)(1) = 2.46622291347846457779...
zeta_fail^(chi)(1) = 2.46752913273870334016...
zeta_suff^(J)(1)   = 2.44257571477179109710...
zeta_fail^(J)(1)   = 2.46752736855057822496...
zeta_max^(F1)      = 2.46752922945601223333...
```

Differences are ~1e-15 absolute, attributable to `sp.nsolve` vs `FindRoot`
on `y*tan(y)=37` and to `simplify`/`limit` vs analytic-limit on `zetaMax`.
This is consistent agreement — but, per F2, the agreement is largely
forced by formula identity, not by independent derivation. The
`zetaMax` analytic-vs-numerical comparison (SymPy 2.4675292294560112350
vs Mathematica 2.4675292294560122333) is the one piece of substantive
cross-check; the ~1e-15 gap is within float round-off of `aF1` evaluation.

No `engine_disagreement` finding.

## Verdict justification

Verdict: `findings` (3). The script computes the correct numbers (engines
agree to round-off), but its only SymPy assertion is tautological (F1),
its only substantive Mathematica assertions reference SymPy's printed
output as targets (F2), and the "FINAL LEDGER" ordering relations the
prose makes load-bearing are only partially asserted (F3). Attacks tried:
(a) substituting wrong `Pe_*` constants into the SymPy script — the
limit-saturation assertion still PASSes, demonstrating F1; (b) checking
whether `expectApprox` targets in `.wl` lines 71-74 appear anywhere
upstream as derived values — they appear only in
`scripts/output/.../stage080..._sympy_audit.txt`, confirming F2;
(c) checking whether the J-pair ordering and J-vs-chi conservative-floor
ordering are asserted anywhere — they are not (F3). Not `UNFIXABLE` (the
math is fine; the assertions just don't anchor to it). Not
`CRITICAL_DOWNSTREAM` (the fixes add assertions and re-anchor targets; the
output numbers don't change, so downstream units consuming these
thresholds are unaffected).

## Self-test notes

Traps checked: (1) Variable independence — F1's proposed
`omega_explicit.subs(Pe, pe_val)` substitutes a concrete float for the only
symbol in the integrand, so the recomputed `zeta_val` depends on `c` via
`pe_val`; trivially non-tautological. (2) Symmetry/parity — N/A, no
integrals here. (3) Trivial-case pre-check — at `c = 96.5285...`,
`pe_val ≈ 96.5`, so `exp(96.5)` dominates and `Omega ≈ pi/2`, giving
`zeta ≈ aF1 * pi^2/4 ≈ 2.4675`, consistent with the printed
`2.466222...`. The directive's targets are the SymPy output values, so
agreement is mechanical. (4) Path specifications — F1 edits
`scripts/.../...sympy_audit.py`, F2 and F3 edit
`mathematica/.../...mathematica_audit.wl`; both extant and correctly
named.
