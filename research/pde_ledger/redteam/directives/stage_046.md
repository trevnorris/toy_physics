---
unit_id: 046
batch: III.1
created_at: 2026-05-22T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-22T12:52:41-06:00
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 046

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl:55-89`

**Issue:**

The Mathematica script is a line-by-line transliteration of the SymPy script: identical hand-typed polynomial coefficients (`pR`, `p1`, `p2`) in the same ordering as the SymPy `P_R`, `P1`, `P2`; identical hand-typed `dGExpected`, `dFExpected`, `gDiffExpected`, `fDiffExpected` as the SymPy `dG_expected`, `dF_expected`, `G_diff_expected`, `F_diff_expected`; mirrored helper functions (`expectZero`, `expectPositiveCoefficients`); identical variable choreography. This defeats the two-engine policy. A typo in any hand-typed expected form lands in both files and passes both engines undetected.

**Required change:**

Replace lines 55-89 of `moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl` so that the Mathematica script performs independent CAS reasoning rather than re-asserting the SymPy script's hand-typed expected forms. Keep lines 1-54 (helpers, banner, symbols, `gTr`/`fTr`/`gFlat`/`fFlat` declarations, strong-split-endpoint checks). Replace lines 55-89 with:

```mathematica
(* Derivatives obtained directly from gTr/fTr; no hand-typed comparison polynomial *)
dGdR = Together[D[gTr, r]];
dFdR = Together[D[fTr, r]];
Print["dG_tr/dR = ", fmt[Factor[dGdR]]];
Print["dF_tr/dR = ", fmt[Factor[dFdR]]];

(* Sign of dG/dR on the bound domain: must be <= 0 on 0 <= r <= 1.
   Use Reduce on the closure (here we test strict negativity on the open interval). *)
dGSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, dGdR < 0], Reals];
Print["Reduce[dG/dR < 0 on (0,1)^3] = ", fmt[dGSignClaim]];
If[TrueQ[dGSignClaim], pass["dG/dR < 0 on bound domain"],
  fail["dG/dR < 0 on bound domain", dGSignClaim]];

dFSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, dFdR > 0], Reals];
Print["Reduce[dF/dR > 0 on (0,1)^3] = ", fmt[dFSignClaim]];
If[TrueQ[dFSignClaim], pass["dF/dR > 0 on bound domain"],
  fail["dF/dR > 0 on bound domain", dFSignClaim]];

(* Branch-difference forms derived directly; no hand-typed p1/p2/gDiffExpected/fDiffExpected *)
deltaG = Together[gTr - gFlat];
deltaF = Together[fFlat - fTr];
Print["G_tr - G_flat = ", fmt[Factor[deltaG]]];
Print["F_flat - F_tr = ", fmt[Factor[deltaF]]];

(* Verify the (1 - r^2) factor in deltaG via polynomial division *)
{deltaGQuotient, deltaGRemainder} =
  PolynomialQuotientRemainder[Together[Numerator[Factor[deltaG]]], 1 - r^2, r];
Print["deltaG numerator / (1 - r^2) quotient = ", fmt[Factor[deltaGQuotient]]];
Print["deltaG numerator / (1 - r^2) remainder = ", fmt[Factor[deltaGRemainder]]];
expectZero["(1 - r^2) divides numerator of G_tr - G_flat", deltaGRemainder];

(* Verify the (1 - r) factor in deltaF via polynomial division *)
{deltaFQuotient, deltaFRemainder} =
  PolynomialQuotientRemainder[Together[Numerator[Factor[deltaF]]], 1 - r, r];
Print["deltaF numerator / (1 - r) quotient = ", fmt[Factor[deltaFQuotient]]];
Print["deltaF numerator / (1 - r) remainder = ", fmt[Factor[deltaFRemainder]]];
expectZero["(1 - r) divides numerator of F_flat - F_tr", deltaFRemainder];

(* Sign of branch difference on the bound domain. Both should be > 0 on (0,1). *)
deltaGSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, deltaG > 0], Reals];
Print["Reduce[G_tr - G_flat > 0 on (0,1)^3] = ", fmt[deltaGSignClaim]];
If[TrueQ[deltaGSignClaim], pass["G_tr > G_flat on bound domain"],
  fail["G_tr > G_flat on bound domain", deltaGSignClaim]];

deltaFSignClaim =
  Reduce[ForAll[{r, xi, delta}, 0 < r < 1 && 0 < xi < 1 && delta > 0, deltaF > 0], Reals];
Print["Reduce[F_flat - F_tr > 0 on (0,1)^3] = ", fmt[deltaFSignClaim]];
If[TrueQ[deltaFSignClaim], pass["F_flat > F_tr on bound domain"],
  fail["F_flat > F_tr on bound domain", deltaFSignClaim]];
```

Do NOT include any hand-typed `pR = 4*r^4*xi^3 + ...`, `p1 = 18*r^2*delta^2*xi + ...`, `p2 = 18*r^3*delta^2*xi^2 + ...`, `dGExpected = -36*r*xi^2*...`, `dFExpected = ...`, `gDiffExpected = 18*xi^2*(1 - r^2)...`, or `fDiffExpected = 4*xi*(1 - r)*...` literals.

If Mathematica's `Reduce[ForAll[...]]` for any of the four sign claims returns a non-`True` result (e.g., it returns `Reduce`'s reduced predicate), keep that result as the printed output, but the corresponding `If[TrueQ[...], pass, fail]` will trigger `fail` and `Exit[1]`. In that case, report the unreduced predicate in the `## Applied: F1` block under `deviation` and stop — do NOT mask the failure.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 046`. The new output must contain `PASS: dG/dR < 0 on bound domain`, `PASS: dF/dR > 0 on bound domain`, `PASS: G_tr > G_flat on bound domain`, `PASS: F_flat > F_tr on bound domain`, `PASS: (1 - r^2) divides numerator of G_tr - G_flat`, `PASS: (1 - r) divides numerator of F_flat - F_tr`. Exit code 0. The verifier will also grep the script for the strings `pR = 4*r^4*xi^3`, `p1 = 18*r^2*delta^2*xi`, `p2 = 18*r^3*delta^2*xi^2`, `dGExpected = -36*r*xi^2`, `dFExpected = ` and confirm all are absent.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl`
- summary: Replaced the hand-typed comparison-polynomial checks with direct Mathematica derivative, factor-divisibility, and Reduce-based sign checks.
- deviation: none

## F2 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py:142-148`

**Issue:**

The script's filename and section 4 heading advertise "tracking-branch bounds". The bound `G_flat <= G_tr <= xi` and `1/(1-xi) <= F_tr <= F_flat` for `0 <= R <= 1` is currently stated only via `print` (sympy:143-147). The monotonicity ingredients (positivity of `P_R`, `P1`, `P2`) and the strong-split endpoints are checked, but no inequality assertion ties them to the bound. A regression that inverted the sign of `(1 - R)` in `F_diff_expected` would still pass the existing `expect_zero("F_flat - F_tr formula", delta_F - F_diff_expected)` check because both sides would cancel the wrong sign.

**Required change:**

In `moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`, insert the following block AFTER line 141 (after `expect_positive_coefficients("P2", P2, R, delta, xi)`) and BEFORE the existing `banner("4. Endpoint bounds")` on line 143:

```python
banner("3b. Sign verification of branch-difference factors")

# Boundary values: at R=1 the tracking and flat branches must coincide; at R=0
# the tracking branch must hit the strong-split endpoint.
expect_zero(
    "G_tr - G_flat vanishes at R=1",
    sp.simplify((G_tr - G_flat).subs(R, 1)),
)
expect_zero(
    "F_flat - F_tr vanishes at R=1",
    sp.simplify((F_flat - F_tr).subs(R, 1)),
)
expect_zero(
    "G_tr at R=0 equals xi",
    sp.simplify(G_tr.subs(R, 0) - xi),
)
expect_zero(
    "F_tr at R=0 equals 1/(1-xi)",
    sp.simplify(F_tr.subs(R, 0) - 1 / (1 - xi)),
)

# Numerical sign sampling of delta_G = G_tr - G_flat and delta_F = F_flat - F_tr
# at three representative interior points (R, xi, delta). All must be strictly
# positive for the bound G_flat <= G_tr and F_tr <= F_flat to hold on the open
# interval R in (0,1).
delta_G_value = sp.simplify(G_tr - G_flat)
delta_F_value = sp.simplify(F_flat - F_tr)
sample_points = [
    {R: sp.Rational(1, 4), xi: sp.Rational(1, 3), delta: sp.Rational(1, 2)},
    {R: sp.Rational(1, 2), xi: sp.Rational(1, 2), delta: sp.Rational(1, 4)},
    {R: sp.Rational(3, 4), xi: sp.Rational(1, 5), delta: sp.Rational(2, 3)},
]
for idx, pt in enumerate(sample_points, start=1):
    g_sample = sp.simplify(delta_G_value.subs(pt))
    f_sample = sp.simplify(delta_F_value.subs(pt))
    print(f"sample {idx}: G_tr - G_flat = {g_sample}, F_flat - F_tr = {f_sample}")
    if g_sample <= 0:
        raise AssertionError(
            f"sample {idx} violates G_tr > G_flat: got {g_sample}"
        )
    if f_sample <= 0:
        raise AssertionError(
            f"sample {idx} violates F_tr < F_flat: got {f_sample}"
        )
```

Apply the analogous additions to `moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl` AFTER the corresponding sign-check block introduced by F1 (i.e., at the end of the section F1 inserts, before the final `Print[""]; Print["Stage 046 Mathematica audit passed."]; Exit[0];`). Use Mathematica idioms:

```mathematica
banner["3b. Boundary-value sign checks"];
expectZero["G_tr - G_flat vanishes at r=1", FullSimplify[(gTr - gFlat) /. r -> 1, Assumptions -> $Assumptions]];
expectZero["F_flat - F_tr vanishes at r=1", FullSimplify[(fFlat - fTr) /. r -> 1, Assumptions -> $Assumptions]];
expectZero["G_tr at r=0 equals xi", FullSimplify[(gTr /. r -> 0) - xi, Assumptions -> 0 < xi < 1 && delta > 0]];
expectZero["F_tr at r=0 equals 1/(1-xi)", FullSimplify[(fTr /. r -> 0) - 1/(1 - xi), Assumptions -> 0 < xi < 1 && delta > 0]];

samplePts = {
  {r -> 1/4, xi -> 1/3, delta -> 1/2},
  {r -> 1/2, xi -> 1/2, delta -> 1/4},
  {r -> 3/4, xi -> 1/5, delta -> 2/3}
};
Do[
  Module[{pt, gs, fs},
    pt = samplePts[[k]];
    gs = FullSimplify[(gTr - gFlat) /. pt];
    fs = FullSimplify[(fFlat - fTr) /. pt];
    Print["sample ", k, ": G_tr - G_flat = ", fmt[gs], ", F_flat - F_tr = ", fmt[fs]];
    If[!TrueQ[gs > 0], fail["sample " <> ToString[k] <> " violates G_tr > G_flat", gs]];
    If[!TrueQ[fs > 0], fail["sample " <> ToString[k] <> " violates F_tr < F_flat", fs]];
  ],
  {k, 1, Length[samplePts]}
];
```

Keep the existing `banner("4. Endpoint bounds")` print block (sympy:143-148) and the existing `Print[""]; Print["Stage 046 Mathematica audit passed."]` (wl:91-93) intact — they document the bound for the reader.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-sympy 046` and `redteam exec-mathematica 046`. Each new output must contain:

SymPy:
- `G_tr - G_flat vanishes at R=1 = 0`
- `F_flat - F_tr vanishes at R=1 = 0`
- `G_tr at R=0 equals xi = 0`
- `F_tr at R=0 equals 1/(1-xi) = 0`
- Three `sample N: G_tr - G_flat = <positive value>, F_flat - F_tr = <positive value>` lines

Mathematica:
- `PASS: G_tr - G_flat vanishes at r=1`
- `PASS: F_flat - F_tr vanishes at r=1`
- `PASS: G_tr at r=0 equals xi`
- `PASS: F_tr at r=0 equals 1/(1-xi)`
- Three `sample N: G_tr - G_flat = <positive value>, F_flat - F_tr = <positive value>` lines

Both must exit 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage046_tracking_branch_bounds_mathematica_audit.wl`
- summary: Added boundary-value zero checks and three interior positive sample checks for the tracking branch differences in both audit scripts.
- deviation: none
