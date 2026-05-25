---
unit_id: 079
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 079

## Per-finding outcomes

### F1 — hardcoded_result

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl:67-70` — the two `expectApprox` calls against the decimal constants `1.00005192880219492404747131934` and `2.46752922945601123498982913352` were removed and replaced with:

```
zeta0Sym = FullSimplify[Limit[aF1*omega^2, pe -> 0, Direction -> "FromAbove"], Assumptions -> $Assumptions];
zetaInfSym = FullSimplify[Limit[aF1*omega^2, pe -> Infinity], Assumptions -> $Assumptions];
expectApprox["zeta_F1(0+) - A_F1", N[zeta0Sym - aF1, 50], 0, 10^-40];
expectApprox["zeta_F1(inf) - A_F1 Pi^2/4", N[zetaInfSym - aF1*Pi^2/4, 50], 0, 10^-40];
```

**Assessment:**
Edit matches the directive exactly — Codex picked the explicitly-permitted `expectApprox`-with-target-0/tol-`10^-40` fallback (justified by the directive itself, since `aF1` is a high-precision numeric atom rather than a symbol). The decimal targets `1.00005…` and `2.46752…` no longer appear anywhere in the `.wl`. The new check exercises `Limit` taken on the product `aF1*omega^2`, so it is not a round-trip of SymPy's output. Captured `.txt` output (lines 29-32) shows:
```
zeta_F1(0+) - A_F1 diff = 0``49.69894745252931
PASS: zeta_F1(0+) - A_F1
zeta_F1(inf) - A_F1 Pi^2/4 diff = 0``49.306707698469
PASS: zeta_F1(inf) - A_F1 Pi^2/4
```
The residuals are exact zero at 50-digit precision (`0``49.7…` is mpmath/Mathematica's notation for "zero with ~49 digits of precision"), well inside the `10^-40` tolerance. Non-tautological: a wrong symbolic limit of `aF1*omega^2` would have produced a nonzero residual at 50-digit precision.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
Same edit as F1 (directive explicitly notes the two findings share lines 67-68). The pre-substituted scalar `zeta0 = N[aF1*omega0^2, 50]` and `zetaInf = N[aF1*omegaInf^2, 50]` (lines 61-62) and their `Print` calls (65-66) were preserved per the directive's instruction to retain them as diagnostic-only output. The assertion lines now compute `Limit[aF1*omega^2, ...]` directly, so Mathematica's `Limit` is exercised on the product rather than on a pre-evaluated `omega0` / `omegaInf`.

**Assessment:**
Resolves the tautology: the new assertion would catch a bug in Mathematica's `Limit[aF1*omega^2, pe -> 0]` (or `pe -> Infinity`) because the substitution `omega -> 1` / `omega -> Pi/2` no longer happens before the assertion runs. The diagnostic `zeta0` / `zetaInf` prints retained at lines 61-62 are no longer assertion-bearing. Outputs above confirm the `Limit` is exercised correctly (residuals identically zero).

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Three lines inserted after line 76 (the existing small-Pe linear coefficient check), matching the directive verbatim:

```
omegaPrime0 = FullSimplify[Limit[D[omega, pe], pe -> 0, Direction -> "FromAbove"], Assumptions -> $Assumptions];
Print["Omega'(0+) = ", fmt[omegaPrime0]];
expectApprox["Omega'(0+) - (4-Pi)/(2 Pi)", N[omegaPrime0 - (4 - Pi)/(2*Pi), 50], 0, 10^-40];
```

**Assessment:**
The script now contains `D[omega, pe]` (it did not before; confirmed via diff). The captured `.txt` output (lines 37-43):
```
Omega'(0+) = -1/2 + 2/Pi
…
Omega'(0+) - (4-Pi)/(2 Pi) diff = 0``99.43378662364455
PASS: Omega'(0+) - (4-Pi)/(2 Pi)
```
Mathematica's `FullSimplify` reduced the limit of the symbolic derivative to `-1/2 + 2/Pi`, which equals `(4-Pi)/(2*Pi)` (algebraically: `-1/2 + 2/Pi = (-Pi + 4)/(2*Pi) = (4-Pi)/(2*Pi)`). The residual is exact zero at extended precision. This forces Mathematica's symbolic-differentiation engine to corroborate the closed-form `Omega(Pe)` through a different code path than SymPy's `sp.series` (which underlies A5). Non-tautological: a wrong closed-form for `omega` would yield a `D[omega, pe]` whose `pe -> 0` limit disagreed with `(4-Pi)/(2*Pi)`.

The benign `N::meprec` warning (`$MaxExtraPrecision = 50. reached`) occurs because Mathematica's `N[]` of `-1/2 + 2/Pi - (4-Pi)/(2*Pi)` collapses to symbolic zero before numeric padding can complete; it does not invalidate the PASS.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script was modified for unit 079 (all three findings target the Mathematica `.wl` only); no SymPy exec log was captured this iteration. The pre-existing SymPy `.txt` output is unchanged and still shows all A1–A5 assertions consistent. No log required.

**Mathematica:** exit=n/a. The captured exec log `redteam/exec_logs/stage_079_mathematica.log` is missing from `exec_logs/`, but the refreshed canonical output `mathematica/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.txt` (mtime 2026-05-22 23:35, post-fix) contains the script's final line `Stage 079 Mathematica audit passed.` and all relevant `PASS:` lines:
```
PASS: Robin residual
PASS: Omega(0+) - 1
PASS: Omega(inf) - Pi/2
PASS: zeta_F1(0+) - A_F1
PASS: zeta_F1(inf) - A_F1 Pi^2/4
PASS: small-Pe constant coefficient check
PASS: small-Pe linear coefficient check
PASS: Omega'(0+) - (4-Pi)/(2 Pi)
```
The script ends with `Exit[0]` which is only reached when no `fail[]` triggers (each `fail` calls `Exit[1]`), so the captured output implies a clean exit. The two `Limit::alimv` warnings and the one `N::meprec` warning are benign (the auditor already flagged the `Limit::alimv` notes as non-blocking).

**Output freshness:** `mathematica/output/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.txt` mtime is 2026-05-22 23:35, newer than the script `mathematica/moving_throat_pde_stage079_family1_quadrupole_pe_map_mathematica_audit.wl` mtime 2026-05-22 23:34. SymPy `.txt` mtime 2026-05-22 23:35 (also refreshed, even though no SymPy edit was needed). Outputs are confirmed post-fix.

## Material-change assessment

`material_change`: false.

Rationale: F1/F2 changed only how Mathematica tests the bracket endpoints (numeric round-trip → symbolic-limit residual); the printed `zeta_F1(0+) = 1.0000519…` and `zeta_F1(inf) = 2.4675292…` values are unchanged (still come from the unchanged `zeta0 = N[aF1*omega0^2, 50]` and `zetaInf = N[aF1*omegaInf^2, 50]` diagnostic prints). F3 added a new assertion but did not alter any pre-existing printed value. No downstream unit's quoted constant is affected.

## Side observations (non-blocking)

- The SymPy script's `nsolve(tol=1e-34)` still does not actually deliver 34-digit precision — its `Robin residual = -1.42e-14` matches mpmath's default 15-digit precision, not the claimed tolerance. The auditor already noted this as a non-blocking note; nothing changed this iteration. Out of scope for this verification.
- The Mathematica output's `Omega'(0+)` print is shown as `-1/2 + 2/Pi`, which is algebraically `(4-Pi)/(2*Pi)` but in a different canonical form. `FullSimplify` chose the additive form; this is purely cosmetic and the residual check still confirms equality.

## Verdict justification

All three findings (F1, F2, F3) are `resolved`. Codex applied the directive verbatim with the explicitly-permitted `expectApprox(target=0, tol=10^-40)` fallback for F1/F2. The hardcoded decimal targets are gone, the new assertions force `Limit` on the product `aF1*omega^2` and `D[omega, pe]` (independent code paths), all eight PASS lines appear in the refreshed Mathematica output, and the script reaches `Exit[0]`. No collateral edits outside the target lines. No regressions in the diff. No material change to printed numeric values.

stage 079: verified
