---
unit_id: 051
batch: III.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 051 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.txt`

## What the script claims to verify

The scripts verify five claims about the tracking-branch lowest-twin sufficiency criterion at the Stage 050/034 checkpoint:

1. The product `Pi_tr = F_tr * G_tr` collapses to a stated closed form (a ratio of polynomials in `xi`, `delta`, `R`).
2. The endpoint limits of `Pi_tr` on the stable branch: `Pi_tr(xi -> 0+) = 0` and `Pi_tr(xi -> 1-) = +oo`.
3. The `zeta_req` formula vanishes at `Pi = C_mix` and equals 1 at `Pi = 2 C_mix`, with `C_mix = 8 Lambda (1-eps) / pi^2`.
4. Three equivalent threshold scales (`Lambda_twin,req`, `M_mix^(twin,req) = G_tr/2`, `Z_W^(twin,req)`) are consistent with each other through the Stage 047/030 coherent map `M_mix = 8 Z_W (1+chi0)^2 / [pi^2 (1-eps_eta)(1-eps)]`.
5. The closed-form quadratic root `xi_(2x)` solves `G_tr(xi_(2x)) = 2 M_mix` exactly.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 75 | `simplify(Pi_tr - Pi_expected) == 0` | yes |
| A2 | sympy | 80-81 | `Pi0 == 0` (limit) | yes |
| A3 | sympy | 84-85 | `Pi1 is oo` (limit) | yes |
| A4 | sympy | 99 | `zeta_req.subs(Pi, Cmix) == 0` | yes |
| A5 | sympy | 100 | `zeta_req.subs(Pi, 2 Cmix) - 1 == 0` | yes |
| A6 | sympy | 104-107 | `(zeta_req - 1).subs(Pi, 2 Cmix) == 0` | redundant with A5 (algebraically identical) |
| A7 | sympy | 128-131 | `Mmix_from_ZW.subs(ZW, ZW_twin_req) - Gtr/2 == 0` | tautological (see F1) |
| A8 | sympy | 145 | `Gtr.subs(xi, xi_2x) - 2 Mmix == 0` | yes |
| A9 | mathematica | 50 | `expectZero[piTr - piExpected]` | yes |
| A10 | mathematica | 56-57 | limit checks `pi0 == 0`, `pi1 == Infinity` | yes |
| A11 | mathematica | 68 | `zeta_req at Pi = C_mix` | yes |
| A12 | mathematica | 69-70 | `zeta_req at 2 C_mix minus 1` (two forms, redundant) | yes (redundant pair) |
| A13 | mathematica | 80 | `mMixFromZW /. zW -> zWTwinReq - gTr/2` | tautological (see F1) |
| A14 | mathematica | 87 | `gTr /. xi -> xi2x - 2 mMix` | yes |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py:115-131`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:72-80`

**What's wrong:**
The assertion `Mmix_from_ZW.subs(ZW, ZW_twin_req) - Gtr/2 == 0` (sympy L128-131; equivalent at wl L80) is algebraically tautological. The script defines:

- `Mmix_twin_req := Gtr/2` (sympy L116)
- `ZW_twin_req := pi**2 (1 - eps_eta)(1 - eps) Gtr / [16 (1 + chi0)**2]` (sympy L117)
- `Mmix_from_ZW := 8 ZW (1 + chi0)**2 / [pi**2 (1 - eps_eta)(1 - eps)]` (sympy L127)

`Mmix_from_ZW` is, by inspection, the algebraic inverse of the linear map that produced `ZW_twin_req` from `Gtr/2`. Substituting `ZW -> ZW_twin_req` simply applies `f^{-1}(f(Gtr/2)) = Gtr/2`. The check is `f^{-1}(f(x)) = x` and cannot fail no matter what physics underlies the coherent map — it tests only the trivial fact that the two formulas in the script are inverses of each other by construction, not that the `M_mix <-> Z_W` bridge from Stage 047/030 is consistent with the threshold derivation.

The docstring (sympy L20-23) explicitly claims: "The `Z_W` threshold is checked through the same Stage 047 coherent map, so the `M_mix <-> Z_W` bridge stays explicit instead of being silently replayed." The current assertion does NOT exercise that bridge — it only verifies an inverse algebra by construction.

**Why this matters:**
At a checkpoint stage with `is_checkpoint: True`, the bar requires substantive assertions. A sign error, wrong prefactor, or wrong power on `(1+chi0)` in the Stage 047/030 coherent map (the upstream physical input) would not be caught by this assertion because the script supplies BOTH the forward map (implicit in `ZW_twin_req`) AND the backward map (`Mmix_from_ZW`) with the same coefficients in the same script. The check holds independent of whether either coefficient is the correct physical one.

**Required change:**
Replace the tautological round-trip with a substantive identity. Define a single forward map only (the canonical Stage 047 coherent map) and derive the threshold from `M_mix^(twin,req) = G_tr/2` via that map alone — then verify the closed form against the explicitly written `ZW_twin_req` (which is the claimed physical answer). Concretely, at sympy L127-131:

Before:
```python
Mmix_from_ZW = sp.simplify(8 * ZW * (1 + chi0) ** 2 / (pi**2 * (1 - eps_eta) * (1 - eps)))
expect_zero(
    "M_mix(Z_W^(twin,req)) - G_tr/2",
    Mmix_from_ZW.subs(ZW, ZW_twin_req) - Gtr / 2,
)
```

After (sketch):
```python
# Stage 047 coherent map (forward): Z_W = pi^2 (1-eps_eta)(1-eps) M_mix / [8 (1+chi0)^2]
ZW_from_Mmix = sp.simplify(pi**2 * (1 - eps_eta) * (1 - eps) * Mmix / (8 * (1 + chi0) ** 2))
# Apply the forward map to the M_mix threshold M_mix = G_tr/2.
ZW_threshold_via_map = sp.simplify(ZW_from_Mmix.subs(Mmix, Gtr / 2))
# Verify that the closed-form Z_W threshold quoted above equals the forward-mapped value.
expect_zero(
    "Z_W^(twin,req) - forward-map(M_mix=G_tr/2)",
    ZW_twin_req - ZW_threshold_via_map,
)
```

This still uses a single coefficient (`pi^2 (1-eps_eta)(1-eps) / [8 (1+chi0)^2]`), but the assertion now compares two independent symbolic expressions (`ZW_twin_req` as written verbatim in L117 vs. the forward-mapped image of `G_tr/2`) — both supplied separately. If the closed-form `ZW_twin_req` had a wrong prefactor (e.g. `8` vs `16`, or `(1+chi0)` vs `(1+chi0)^2`), the residual would not collapse to zero.

Apply the symmetric edit to the Mathematica script at wl L74-80.

**Verification:**
After the edit, sympy output should print `Z_W^(twin,req) - forward-map(M_mix=G_tr/2) = 0` (replacing the prior `M_mix(Z_W^(twin,req)) - G_tr/2 = 0` line). Mathematica output should print the same. The check should fail (residual nonzero) if `ZW_twin_req`'s prefactor is perturbed by any nontrivial factor.

### F2 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage051_lowest_twin_criterion_mathematica_audit.wl:34-87`

**What's wrong:**
The `.wl` script is structurally a line-by-line port of the `.py` script with identical variable choreography. Three concrete examples:

1. SymPy L55-60 defines `Gtr` and `Ftr` as explicit ratios of polynomials; Mathematica wl L34-39 defines `gTr` and `fTr` as the same explicit ratios, with the same numerator/denominator factor ordering. Neither engine derives `Pi_tr = F_tr * G_tr` from a more primitive Stage 045 quantity — both supply the same pre-baked closed forms verbatim.

2. SymPy L71-74 defines `Pi_expected` as an explicit closed form; Mathematica wl L41-45 defines `piExpected` as the same expression copy-pasted with `R` -> `r` and `**` -> `^`. Both then run `simplify(piTr - piExpected) == 0`. The assertion is genuine, but both engines are checking the same pre-supplied target against the same pre-supplied product. No independent factoring path is exercised.

3. SymPy L134-141 hands in the closed-form quadratic root `xi_(2x) = [...]/18` directly; Mathematica wl L82-85 hands in the same closed root transliterated. Neither engine uses `sp.solve` / `Solve` to derive the root from `Gtr == 2 Mmix`. Both then verify `Gtr.subs(xi, xi_2x) - 2 Mmix == 0` using the same handed-in form.

The Mathematica script does not provide independent derivation; it echoes the SymPy algebra. At a checkpoint stage (`is_checkpoint: True`) this violates the second-engine policy: both engines should derive the result independently from the physical premises, not echo each other.

**Why this matters:**
A typo or sign error in either pre-baked closed form (`Pi_expected`, `xi_(2x)`) would be replicated in both engines and would not be caught by either's assertions. The Mathematica check `Gtr(xi_2x) - 2 mMix == 0` would still pass if the closed-form root were missing a term, provided the same defective form were in both scripts. An independent re-derivation in Mathematica via `Solve[gTr == 2 mMix, xi]` would catch such a co-replicated error.

**Required change:**
In the Mathematica script, replace the two transliterated closed-form supplies with independent derivations:

1. At wl L82-85, replace the handed-in `xi2x` closed form with an independent solve. Before:
```mathematica
xi2x = FullSimplify[
  (2 mMix (9 + 2 r^2) - 9 delta + Sqrt[(2 mMix (9 + 2 r^2) - 9 delta)^2 + 648 mMix delta])/18,
  Assumptions -> $Assumptions
];
Print["xi_(2x) = ", fmt[xi2x]];
expectZero["G_tr(xi_(2x)) - 2 M_mix", (gTr /. xi -> xi2x) - 2 mMix];
```

After:
```mathematica
(* Independently derive the positive root of gTr == 2 mMix using Solve. *)
xi2xRoots = xi /. Solve[gTr == 2 mMix, xi];
xi2xDerived = FullSimplify[
  Select[xi2xRoots, Simplify[# > 0, $Assumptions] === True &][[1]],
  Assumptions -> $Assumptions
];
(* Independent closed-form target (the claimed answer in the docstring). *)
xi2xClaim = FullSimplify[
  (2 mMix (9 + 2 r^2) - 9 delta + Sqrt[(2 mMix (9 + 2 r^2) - 9 delta)^2 + 648 mMix delta])/18,
  Assumptions -> $Assumptions
];
Print["xi_(2x) (Solve) = ", fmt[xi2xDerived]];
Print["xi_(2x) (claim) = ", fmt[xi2xClaim]];
expectZero["xi_(2x): Solve vs claim", xi2xDerived - xi2xClaim];
expectZero["G_tr(xi_(2x)) - 2 M_mix (Solve root)", (gTr /. xi -> xi2xDerived) - 2 mMix];
```

2. At wl L41-50, replace the handed-in `piExpected` round-trip with an independent factoring path. Before:
```mathematica
piTr = FullSimplify[fTr gTr, Assumptions -> $Assumptions];
piExpected = FullSimplify[
  xi (xi + delta) (9 delta + (9 + 2 r) xi)^2 (9 delta + (9 + 2 r^2) xi)/
    (9 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2),
  Assumptions -> $Assumptions
];
Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["Pi_tr = ", fmt[piTr]];
expectZero["Pi_tr - expected closed form", piTr - piExpected];
```

After (keep the closed-form claim AS a target, but verify via `Factor`/`Together` on the product separately rather than just `FullSimplify`):
```mathematica
piTrProduct = Together[fTr gTr];
piTrFactored = Factor[piTrProduct];
piExpected = FullSimplify[
  xi (xi + delta) (9 delta + (9 + 2 r) xi)^2 (9 delta + (9 + 2 r^2) xi)/
    (9 (1 - xi) (9 delta^2 + 18 delta xi + (9 + 2 r^2) xi^2)^2),
  Assumptions -> $Assumptions
];
piExpectedFactored = Factor[piExpected];
Print["G_tr = ", fmt[gTr]];
Print["F_tr = ", fmt[fTr]];
Print["Pi_tr (Factor) = ", fmt[piTrFactored]];
Print["Pi_tr (claim, Factor) = ", fmt[piExpectedFactored]];
expectZero["Pi_tr - expected closed form", piTrFactored - piExpectedFactored];
```

This routes the Mathematica verification through `Factor`/`Together` rather than `FullSimplify` (which can opaquely match anything), and so provides an independent algebraic path to canonical form.

**Verification:**
After the edit, the Mathematica output should include the new lines:
- `xi_(2x) (Solve) = ...` and `xi_(2x) (claim) = ...`
- `xi_(2x): Solve vs claim = 0` (PASS)
- `G_tr(xi_(2x)) - 2 M_mix (Solve root) = 0` (PASS)
- `Pi_tr (Factor) = ...` and `Pi_tr (claim, Factor) = ...`

The verifier confirms that the new lines appear, that all `expectZero` checks PASS, and that the script exits 0.

## Independent-derivation check (Mathematica)

The `.wl` is a transliteration of the `.py`. The same closed forms (`Gtr`, `Ftr`, `Pi_expected`, `xi_(2x)`) are supplied verbatim in both, with identical factor ordering and identical algebraic steps. No independent route to the closed forms (e.g., via `Solve` for `xi_(2x)`, or via `Factor`/`Together` for `Pi_tr`) is used. See F2.

Two corresponding sections illustrate:

- py L134-141 / wl L82-85: closed-form quadratic root supplied verbatim, then substituted into `G_tr`.
- py L71-74 / wl L41-45: closed-form `Pi_expected` supplied verbatim, then differenced against `simplify(F_tr * G_tr)`.

## Engine cross-check

Both engines passed all assertions and printed the same closed forms (up to formatting and the sign convention `(xi - 1)` vs `(1 - xi)`, which is purely cosmetic). Examples:

- SymPy `M_mix^(twin,req) = 9*xi*(delta + xi) / [2*(9*delta + (2*R^2 + 9)*xi)]` (output L104-108) matches Mathematica `(9*xi*(delta + xi))/(2*(9*delta + (9 + 2*r^2)*xi))` (output L35).
- SymPy `xi_(2x) = Mmix*(2*R^2 + 9)/9 - delta/2 + sqrt(648*Mmix*delta + (2*Mmix*(2*R^2 + 9) - 9*delta)^2)/18` (output L121-127) matches Mathematica `(-9*delta + 18*mMix + 4*mMix*r^2 + Sqrt[648*delta*mMix + (9*delta - 2*mMix*(9 + 2*r^2))^2])/18` (output L39). Distributing the SymPy form: `Mmix(2R^2+9)/9 = (18 Mmix + 4 Mmix R^2)/18`, and `-delta/2 = -9 delta/18`. Sum of first two: `(18 Mmix + 4 Mmix R^2 - 9 delta)/18`. Adding the sqrt/18 yields the Mathematica form. Agreement holds.

`engines_agree: true`.

## Verdict justification

The engine cross-check holds and all assertions execute. However, the audit fails on two substantive grounds: (1) the `M_mix <-> Z_W` round-trip is tautological by construction and does not exercise the Stage 047/030 coherent map the docstring claims to verify; (2) the Mathematica script is a transliteration of the SymPy script, supplying the same closed-form `Pi_expected` and `xi_(2x)` verbatim rather than deriving them independently via `Solve` or `Factor`. At a checkpoint stage with `is_checkpoint: True`, both findings warrant directives. Verdict: `findings`.

Attacks tried that failed:

- Domain/parity sanity checks: all `simplify` calls operate on rational expressions over symbols declared `positive` (sympy) or constrained by `$Assumptions` (mathematica); no branch cut hazards in the rational algebra.
- xi_(2x) double-root concern: the discriminant `(2 Mmix(9+2R^2) - 9 delta)^2 + 648 Mmix delta > (2 Mmix(9+2R^2)-9 delta)^2` for `Mmix, delta > 0`, so only the `+` root is positive — the script's choice of root is physically forced and not a hidden missing-branch issue.
- `Pi_tr` closed-form expansion: hand-checked `F_tr * G_tr` cancellation (one factor of `D1 = 9 delta + (9+2R^2) xi`) and the result matches `Pi_expected`; A1/A9 are genuine non-tautological identities.
- zeta_req at Pi=C_mix and Pi=2 C_mix: substitution yields `0/(1-eps)` and `1`, respectively — the algebra is fixed and the check non-trivially exercises the formula at distinguished points.

## Self-test notes

For F1's required change, mentally substituted: `ZW_from_Mmix.subs(Mmix, Gtr/2)` yields `pi^2 (1-eps_eta)(1-eps) (Gtr/2) / [8 (1+chi0)^2] = pi^2 (1-eps_eta)(1-eps) Gtr / [16 (1+chi0)^2]`, which matches `ZW_twin_req` (sympy L117) exactly; residual collapses to 0. The check is non-tautological because if `ZW_twin_req`'s denominator were perturbed to `8 (1+chi0)^2` or its numerator to `(1+chi0)`, the residual would not cancel. For F2's required change: `Solve[gTr == 2 mMix, xi]` on `9 xi (xi + delta) / (9 delta + (9+2r^2) xi) == 2 mMix` yields the quadratic `9 xi^2 + [9 delta - 2 mMix (9+2r^2)] xi - 18 mMix delta == 0`, whose positive root by the quadratic formula reproduces the claimed closed form; the `Select[..., # > 0, ...]` step under `0 < xi < 1, mMix > 0, delta > 0` is well-posed because the negative root has numerator `2 mMix (9 + 2 r^2) - 9 delta - sqrt(...) < 0`. Variable-independence and parity traps do not apply (no diff/integral over an unbounded domain in either change). Path specs: `.py` edits land in `scripts/`, `.wl` edits land in `mathematica/`.
