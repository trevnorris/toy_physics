---
unit_id: 081
batch: III.4
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 4
scripts_checked:
  sympy: present
  mathematica: insufficient
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 081 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage081_family1_pi_thresholds_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.txt`

## What the script claims to verify

The docstring says the script translates Family-1 zeta-demand thresholds (numerically obtained at lambda_mu = 1, supplied as constants from Stage 63) into the product variable Pi_tr / C_mix using the Stage-35 inversion of zeta(Pi). Both scripts therefore (a) start from the support-demand law zeta = (Pi - C_mix) / (C_mix - eps_blk (2 C_mix - Pi)), (b) algebraically invert to get Pi(zeta) and hence Q(zeta; eps_blk) = Pi(zeta)/C_mix, (c) check the boundary values Q(0) = 1 and Q(1) = 2, (d) numerically evaluate Q at the five Stage-63 threshold values at eps_blk = 0, and (e) report the blocking ceiling eps_blk < 1/zeta_max^(F1).

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 40 | `expect_zero("Q(0)-1", Q.subs(zeta, 0) - 1)` | yes — exercises the derived Q at boundary |
| A2 | sympy | 41 | `expect_zero("Q(1)-2", Q.subs(zeta, 1) - 2)` | yes — exercises the derived Q at boundary |
| M1 | mathematica | 46 | `expectZero["Q(0)-1", (qq /. zeta -> 0) - 1]` | partial — `qq` is hardcoded, not derived from `piOfZeta` |
| M2 | mathematica | 47 | `expectZero["Q(1)-2", (qq /. zeta -> 1) - 2]` | partial — same |
| M3 | mathematica | 48 | `expectZero["dQ/dzeta exact formula", dqq - (1 - epsBlk)/(1 - epsBlk*zeta)^2]` | no — tautological calculus of hardcoded `qq` |
| M4-M8 | mathematica | 68-72 | `expectApprox[...]` of `qq /. zeta -> zeta_*` against literal target | no — targets are literally `1 + zeta_*` and `qq(zeta_*; 0) = 1 + zeta_*` by construction |
| M9 | mathematica | 76 | `expectApprox["blocking ceiling numeric check", epsCeiling, 0.40526368971137149977, 10^-14]` | no — checks `1/zetaMaxF1` against the typed digits of `1/zetaMaxF1` |

## Findings

### F1 — hardcoded_result

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:41`

**What's wrong:**
The Mathematica script derives `piOfZeta` from the premise via `Solve[zeta == zetaExpr, piTr]` at line 40, but then on line 41 simply assigns

```
qq = FullSimplify[(1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
```

This `qq` is hardcoded; it is never compared back to `FullSimplify[piOfZeta/cMix]`. All subsequent Mathematica checks (`Q(0)-1`, `Q(1)-2`, the eps=0 numeric checks, and the blocking ceiling) operate on `qq` only. The genuine inversion result `piOfZeta` is computed and printed but never assertively connected to `qq`, so the Mathematica side does not actually verify the inversion that the docstring claims to verify; it only verifies arithmetic properties of a literal rational function.

**Why this matters:**
If the inversion formula `piOfZeta` were ever modified (or the underlying `zetaExpr` definition mistyped) the Mathematica audit would still pass because `qq` is independent of `piOfZeta`. The cross-engine guarantee that the Mathematica side is exercising the physics, not just calculus of a written-down expression, is lost.

**Required change:**
Add an explicit `expectZero` between the derived inversion and the closed form: at line 41-42, replace the bare `qq = FullSimplify[...]` with a `qq` defined as `FullSimplify[piOfZeta / cMix]`, and append an `expectZero["Q matches piOfZeta/cMix", qq - (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)]` so the closed form is *tested* against the Solve-derived inversion rather than asserted. See directive for exact edit.

**Verification:**
The output transcript should display a new `PASS: Q matches piOfZeta/cMix` line and the existing `Q(0)-1` and `Q(1)-2` assertions should still pass against the Solve-derived `qq`.

### F2 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:42,48`

**What's wrong:**
Line 42 defines `dqq = FullSimplify[D[qq, zeta], Assumptions -> $Assumptions]`. Line 48 then asserts

```
expectZero["dQ/dzeta exact formula", dqq - (1 - epsBlk)/(1 - epsBlk*zeta)^2];
```

Both sides are nothing but symbolic calculus of the already-known closed form `qq = (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)`. The check cannot fail unless Mathematica's symbolic differentiation engine itself is broken; it does not test any physical claim, nor does it test that `qq` correctly represents the inversion. The output line `dQ/dzeta exact formula = 0` is the result of this trivially true identity.

**Why this matters:**
This is reported as a non-trivial verification in the output transcript (`PASS: dQ/dzeta exact formula`), giving false confidence that the Q function's slope has been physically validated. In reality, no physical content is checked.

**Required change:**
Either delete the `dQ/dzeta exact formula` assertion (line 48), or replace it with an assertion that genuinely couples calculus and physics. The minimally invasive fix is to delete line 48 and the `dqq` definition on line 42, since the Q matches inversion check added by F1 already ties `qq` back to the physics. See directive.

**Verification:**
The output should no longer contain the `dQ/dzeta exact formula = 0` or `PASS: dQ/dzeta exact formula` lines.

### F3 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:68-72`

**What's wrong:**
The five `expectApprox` calls compare `qq /. zeta -> zeta_* /. epsBlk -> 0`, which equals `1 + zeta_*` by direct substitution into the hardcoded `qq`, against literal targets that are exactly `1 + zeta_*` to high precision:

| zeta_* | 1 + zeta_* | target literal in script |
|---|---|---|
| `2.46622291347846` | `3.46622291347846…` | `3.4662229134784601214` |
| `2.46752913273870` | `3.46752913273870…` | `3.4675291327386998930` |
| `2.44257571477179` | `3.44257571477179…` | `3.4425757147717899187` |
| `2.46752736855058` | `3.46752736855058…` | `3.4675273685505798582` |
| `2.46752922945601` | `3.46752922945601…` | `3.4675292294560100537` |

The trailing digits in the targets (e.g. `…01214`, `…98930`) are just the high-precision arithmetic noise of `1 + zeta_*` in 50-digit `N[]`. The check is `1 + zeta_* ≈ 1 + zeta_*` at tolerance `10^-14`. It does not cross-check the result against an independent source; the SymPy script does not produce these targets either — they are just printed.

**Why this matters:**
These checks masquerade as numeric cross-validation between Mathematica's evaluation of Q and an external high-precision target. They are in fact comparing a quantity to itself. If `qq` were defined wrong, the same wrong value would be the target. (Note: this is the same root cause as F1; the targets were generated from the script's own formula, not from an independent computation.)

**Required change:**
Replace the literal target list (lines 68-72) with assertions of the symbolic form `qq /. zeta -> zeta_* /. epsBlk -> 0 == 1 + zeta_*`, i.e. assert that at eps_blk = 0 the explicit formula reduces to `1 + zeta` for each of the five threshold values. This converts the check into an exercise of `qq`'s functional form rather than a check of constants against themselves. See directive.

**Verification:**
The transcript will show five new `expectZero`-style lines confirming `qq(zeta_*; 0) - (1 + zeta_*) == 0` instead of `expectApprox diff` lines against magic numbers.

### F4 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage081_family1_pi_thresholds_mathematica_audit.wl:74-76`

**What's wrong:**
Line 74 computes `epsCeiling = N[1/zetaMaxF1, 40]`. Line 76 then checks `epsCeiling` against `0.40526368971137149977`, which is exactly the decimal expansion of `1/2.46752922945601`. The assertion `expectApprox["blocking ceiling numeric check", epsCeiling, 0.40526368971137149977, 10^-14]` therefore reduces to `1/zetaMaxF1 ≈ digits(1/zetaMaxF1)` — a check of `N[]`'s reciprocal against its own decimal expansion.

**Why this matters:**
The blocking ceiling has physical meaning (denominator positivity requires `eps_blk < 1/zeta_max^(F1)`), but the numeric check does not test any cross-derivation. It is informational only and gives no protection against a wrong `zetaMaxF1` input.

**Required change:**
Replace the literal target with an explicit symbolic check, e.g. `expectZero["blocking ceiling reciprocal", epsCeiling * zetaMaxF1 - 1]` (with appropriate numerical tolerance via Abs). This exercises the reciprocal relationship rather than asserting `x == digits(x)`. See directive.

**Verification:**
The transcript should show a new `Abs[…] < tol` (or `expectZero`) line confirming `epsCeiling * zetaMaxF1 == 1` instead of the current `diff = …e-17` line against a magic number.

## Independent-derivation check (Mathematica)

The Mathematica script is *partially* a transliteration of the SymPy script: it shares the same `zetaExpr` premise (line 39 vs SymPy line 33), the same `Solve` step (line 40 vs SymPy line 34), the same five Stage-63 numerical thresholds, the same eps_blk = 0 evaluation pattern, and the same blocking ceiling. Quoted side-by-side:

SymPy line 33-36:
```
zeta_expr = (Pi - Cmix) / (Cmix - eps_blk * (2 * Cmix - Pi))
Pi_of_zeta = sp.solve(sp.Eq(zeta, zeta_expr), Pi)[0]
Pi_of_zeta = sp.simplify(Pi_of_zeta)
Q = sp.simplify(Pi_of_zeta / Cmix)
```

Mathematica line 39-41:
```
zetaExpr = FullSimplify[(piTr - cMix)/(cMix - epsBlk*(2*cMix - piTr)),...];
piOfZeta = FullSimplify[piTr /. First[Solve[zeta == zetaExpr, piTr]], ...];
qq = FullSimplify[(1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta), Assumptions -> $Assumptions];
```

The Mathematica script departs from a pure transliteration only at line 41, where instead of `qq = FullSimplify[piOfZeta/cMix]` (which would parallel the SymPy `Q = sp.simplify(Pi_of_zeta / Cmix)`), it writes a hand-encoded closed form. This departure is what creates F1 and the cascade F3/F4. The departure is *less* independent, not more — it weakens the second-engine guarantee rather than providing an independent derivation path. I do not file a separate `mathematica_transliteration` finding because the substantive issue is captured by F1; if F1's fix is applied (making `qq = FullSimplify[piOfZeta/cMix]` and then asserting it equals the closed form) the second engine still does the same Solve step as SymPy, but at least both engines then trace back to the same premise without an unverified hand-encoded shortcut.

## Engine cross-check

Both engines report the same final closed form `Q(zeta; eps_blk) = (1 + zeta - 2 eps_blk zeta)/(1 - eps_blk zeta)`:

- SymPy output line 14: `Q(zeta;eps_blk) = (2*eps_blk*zeta - zeta - 1)/(eps_blk*zeta - 1)` — algebraically identical (numerator and denominator both negated).
- Mathematica output line 14: `Q(zeta;eps_blk) = (1 + zeta - 2*epsBlk*zeta)/(1 - epsBlk*zeta)`.

Both engines report the same five eps_blk = 0 values to ~15 digits (e.g. `3.4662229134784601214` from both). The blocking ceiling `0.40526368971137...` matches between SymPy and Mathematica outputs. No cross-engine disagreement.

## Verdict justification

The two engines agree on the closed-form Q and on the numerical evaluations. The SymPy side genuinely derives `Q` from the premise via `solve` and substitutes the (hardcoded-from-Stage-63) zeta thresholds into it. The Mathematica side computes the inversion but does *not* algebraically tie it to the assertions it makes; its assertions instead exercise a hand-written closed form against itself (F1, F2, F3, F4). All four findings are confined to the Mathematica script and are mechanically fixable. Attacks attempted: I checked that `Solve` returns a unique branch (linear in Pi, so yes), that the assumptions (`piTr > 0`, `epsBlk >= 0, epsBlk < 1`) don't silently kill a branch (they don't — both engines arrive at the same Q), that the `Q(1) = 2` identity is non-trivial (it requires the `(2 - 2 eps)/(1 - eps)` cancellation, so yes, it does exercise the formula's structure at one point), and that the five eps=0 targets are not from an independent source (they are not — they are literally `1 + zeta_*`). Verdict: `findings`, four findings, no stop-cold.

## Self-test notes

Traps checked: (1) `Solve` branch ambiguity — none, zetaExpr is linear in piTr so the inversion is unique; (2) `simplify` hiding branch errors under aggressive assumptions — `$Assumptions` includes `zeta >= 0`, `epsBlk < 1` which keeps the denominator nonzero on the physical domain, so no hidden simplification; (3) parity/integration traps — not applicable, no integrals in this unit; (4) the proposed directive replaces hardcoded-constant comparisons with `expectZero` on symbolic differences which I mentally verified reduce to 0 (`qq /. zeta -> zeta_* /. epsBlk -> 0 - (1 + zeta_*) == 0` is immediate from the closed form, and `epsCeiling * zetaMaxF1 - 1 == 0` is immediate by definition of `epsCeiling = 1/zetaMaxF1`).
