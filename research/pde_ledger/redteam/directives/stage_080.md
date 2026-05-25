---
unit_id: 080
batch: III.4
created_at: 2026-05-22T00:00:00Z
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-22T23:37:06-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 080

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:60-72`

**Issue:** The only SymPy assertion (lines 71-72) checks `lim_{lambda_mu -> oo}
zeta_*(lambda_mu) == zeta_max`. For every `c_* > 0`, `c_* * lambda_mu**2 -> +oo`
as `lambda_mu -> oo`, so the limit equals `lim_{Pe -> oo} zeta(Pe) = zeta_max`
regardless of the value of `c_*`. The Stage-61-derived numerical thresholds
(96.5285..., 11220.5..., 22.0062..., 2558.02...) are never exercised. Add
numerical-value assertions that tie `zeta_*(1)` back to an independent
evaluation of `A_F1 * Omega(Pe)**2`.

**Required change:**

In `scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py`,
**after line 60** (the existing `for k, v in vals.items(): print(k, '=', v)`
loop) and **before line 62** (the `# Large-lambda limits saturate...`
comment), insert the following block exactly:

```python

# Independent recomputation of zeta_*(1) via A_F1 * Omega(Pe)^2 with
# explicit float Pe values; the four reference targets are the
# corresponding output values printed above. This anchors the four
# Stage-61 numerical Pe thresholds (96.5285..., 11220.5..., 22.0062...,
# 2558.02...) to the four zeta values, breaking the tautology in the
# limit-saturation check below.
def _omega_explicit(p):
    return sp.pi * p * (2 * p * sp.exp(p) + sp.pi) / ((4 * p**2 + sp.pi**2) * (sp.exp(p) - 1))

A_F1_recomputed = (sp.Rational(12321, 5) + sp.pi**2 / 4) / (sp.Rational(12321, 5) + y_F1**2)
_numeric_targets = [
    ("zeta_suff^(chi)(1)", sp.Float('96.5285247264386', 30), sp.Float('2.4662229134784638979', 25)),
    ("zeta_fail^(chi)(1)", sp.Float('11220.5441626259', 30), sp.Float('2.4675291327387028754', 25)),
    ("zeta_suff^(J)(1)",   sp.Float('22.0062226330754', 30), sp.Float('2.4425757147717912819', 25)),
    ("zeta_fail^(J)(1)",   sp.Float('2558.01892349205', 30), sp.Float('2.4675273685505776147', 25)),
]
for _name, _pe_val, _expected in _numeric_targets:
    _zeta_val = sp.N(A_F1_recomputed * _omega_explicit(_pe_val)**2, 25)
    _diff = abs(sp.N(_zeta_val - _expected, 30))
    print(f"zeta numeric check {_name}: diff = {_diff}")
    if _diff > sp.Float('1e-14'):
        raise AssertionError(f"zeta numeric check {_name}: |{_zeta_val} - {_expected}| = {_diff} > 1e-14")

```

The leading blank line and trailing blank line in the inserted block are
intentional so the file's existing blank-line spacing is preserved.

Leave lines 62-72 (the large-lambda saturation block) unchanged.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 080` and
confirm (a) the new check prints four `zeta numeric check ...: diff = ...`
lines before the `limit zeta_*` lines, (b) the script exits 0, and
(c) grep of the script for `_numeric_targets` finds the new block.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py`
- summary: Added independent zeta numeric checks tying the four Stage-61 Pe thresholds to recomputed `A_F1 * Omega(Pe)^2` values.
- deviation: none

## F2 — hardcoded_result

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:53-74`

**Issue:** The four `expectApprox` targets on lines 71-74 are the SymPy
script's output values copied as literals, so the Mathematica "numeric
check" only confirms the two engines share identical formulas — it is not
an independent verification. Replace the literal targets with values
computed inside Mathematica via a freshly named independent evaluation
path (`omegaIndep`, `aF1Indep`).

**Required change:**

In
`mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`:

Step 1. **After line 53** (`zetaFailJ = zetaFn[peFailJ];`) and **before
line 55** (`Print["zeta_max^(F1) = ", fmt[zetaMax]];`), insert this block
exactly (preserve the surrounding blank line on line 54):

```
omegaIndep[p_] := Pi*p*(2*p*Exp[p] + Pi) / ((4*p^2 + Pi^2) * (Exp[p] - 1));
aF1Indep = N[(12321/5 + Pi^2/4) / (12321/5 + yRoot^2), 50];
zetaTargetSuffChi = N[aF1Indep * omegaIndep[ToExpression["96.5285247264386`30"]]^2, 40];
zetaTargetFailChi = N[aF1Indep * omegaIndep[ToExpression["11220.5441626259`30"]]^2, 40];
zetaTargetSuffJ   = N[aF1Indep * omegaIndep[ToExpression["22.0062226330754`30"]]^2, 40];
zetaTargetFailJ   = N[aF1Indep * omegaIndep[ToExpression["2558.01892349205`30"]]^2, 40];

```

Step 2. **Replace lines 71-74**

```
expectApprox["zeta_suff^(chi)(1) numeric check", zetaSuffChi1, ToExpression["2.4662229134784638979`25"], 10^-14];
expectApprox["zeta_fail^(chi)(1) numeric check", zetaFailChi1, ToExpression["2.4675291327387028754`25"], 10^-14];
expectApprox["zeta_suff^(J)(1) numeric check", zetaSuffJ1, ToExpression["2.4425757147717912819`25"], 10^-14];
expectApprox["zeta_fail^(J)(1) numeric check", zetaFailJ1, ToExpression["2.4675273685505776147`25"], 10^-14];
```

with

```
expectApprox["zeta_suff^(chi)(1) numeric check", zetaSuffChi1, zetaTargetSuffChi, 10^-14];
expectApprox["zeta_fail^(chi)(1) numeric check", zetaFailChi1, zetaTargetFailChi, 10^-14];
expectApprox["zeta_suff^(J)(1) numeric check", zetaSuffJ1, zetaTargetSuffJ, 10^-14];
expectApprox["zeta_fail^(J)(1) numeric check", zetaFailJ1, zetaTargetFailJ, 10^-14];
```

After these two steps, the four numeric targets are derived inside
Mathematica from `aF1Indep` and `omegaIndep` rather than from SymPy's
printed output.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 080`
and confirm (a) all four `numeric check` checks still PASS, (b) the
script exits 0, and (c) `grep -E '2\.46622291347846389|2\.46752913273870287|2\.44257571477179128|2\.46752736855057761' mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl` returns no matches (the literal SymPy-output targets are gone).

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`
- summary: Replaced literal numeric targets with independently computed Mathematica targets through `omegaIndep` and `aF1Indep`.
- deviation: none

## F3 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl:90`

**Issue:** The SymPy script's FINAL LEDGER prose (lines 74-78) implicitly
asserts four ordering facts: (i) `zeta_suff^(chi)(1) < zeta_fail^(chi)(1) <
zeta_max`, (ii) `zeta_suff^(J)(1) < zeta_fail^(J)(1) < zeta_max`,
(iii) `zeta_suff^(J)(1) <= zeta_suff^(chi)(1)`, (iv) `zeta_fail^(J)(1) <=
zeta_fail^(chi)(1)`. Only (i) is asserted (Mathematica line 90). Extend
to all four.

**Required change:**

In
`mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`,
**replace line 90**

```
expectTrue["zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max", zetaSuffChi1 < zetaFailChi1 < zetaMax];
```

with the following four lines:

```
expectTrue["zeta_suff^(chi)(1) < zeta_fail^(chi)(1) < zeta_max", zetaSuffChi1 < zetaFailChi1 < zetaMax];
expectTrue["zeta_suff^(J)(1) < zeta_fail^(J)(1) < zeta_max", zetaSuffJ1 < zetaFailJ1 < zetaMax];
expectTrue["zeta_suff^(J)(1) <= zeta_suff^(chi)(1)", zetaSuffJ1 <= zetaSuffChi1];
expectTrue["zeta_fail^(J)(1) <= zeta_fail^(chi)(1)", zetaFailJ1 <= zetaFailChi1];
```

Do not modify any other line; leave the trailing `Print[""];` /
`Print["Stage 080 Mathematica audit passed."];` / `Exit[0];` block intact.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 080`
and confirm four `PASS:` lines for the inequality checks (instead of the
current one) and the script still exits 0.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage080_family1_zeta_thresholds_mathematica_audit.wl`
- summary: Extended the Mathematica ordering assertions from the chi branch to all four ledger ordering facts.
- deviation: none
