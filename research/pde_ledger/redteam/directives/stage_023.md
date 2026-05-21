---
unit_id: 023
batch: I.2
created_at: 2026-05-21T00:00:00Z
findings_count: 4
stop_cold: null
applied: true
applied_at: 2026-05-21T13:48:45-06:00
findings_applied: 4
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 023

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:247-252`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:130-132`

**Issue:**
The three "Nbar0 / aN0 / bN0 formula" assertions are tautological: `(Nbar0, aN0, bN0)` comes from `grouped_parts(N020, N021, N022)`, which by definition returns the exact RHS the assertion compares against. The checks reduce to `x - x == 0` and cannot fail.

**Required change:**
Replace the three tautological assertions with a single additivity check that exercises the linearity of `grouped_parts` for the N-bundle.

In the SymPy script, replace lines 247-252 with:

Before (sympy lines 247-252):
```
    expect_zero(
        "Nbar0 formula",
        Nbar0 - (N020 + 2 * N021 + 2 * N022) / 5,
    )
    expect_zero("aN0 formula", aN0 - (2 * N020 - N021 - N022) / 10)
    expect_zero("bN0 formula", bN0 - (N021 - N022) / 2)
```

After:
```
    # Non-tautological linearity test: grouped_parts(N0 + N2 lane sum)
    # equals the sum of the individual grouped_parts outputs.
    Nbar02, aN02, bN02 = grouped_parts(N020 + N220, N021 + N221, N022 + N222)
    expect_zero("Nbar0 + Nbar2 additivity", Nbar02 - (Nbar0 + Nbar2))
    expect_zero("aN0 + aN2 additivity", aN02 - (aN0 + aN2))
    expect_zero("bN0 + bN2 additivity", bN02 - (bN0 + bN2))
```

In the Mathematica script, replace lines 130-132 with:

Before (mathematica lines 130-132):
```
expectZero["Nbar0 formula", nbar0 - (n020 + 2*n021 + 2*n022)/5];
expectZero["aN0 formula", aN0 - (2*n020 - n021 - n022)/10];
expectZero["bN0 formula", bN0 - (n021 - n022)/2];
```

After:
```
(* Non-tautological linearity test: grouped_parts of a lane sum
   equals the sum of the individual grouped_parts outputs. *)
{nbar02, aN02, bN02} = groupedParts[n020 + n220, n021 + n221, n022 + n222];
expectZero["Nbar0 + Nbar2 additivity", nbar02 - (nbar0 + nbar2)];
expectZero["aN0 + aN2 additivity", aN02 - (aN0 + aN2)];
expectZero["bN0 + bN2 additivity", bN02 - (bN0 + bN2)];
```

Note: the Mathematica script does NOT currently declare `n220, n221, n222` symbols in its `Clear[...]` and `$Assumptions` lines (line 95-96). After replacement, add them to both. Specifically:

In Mathematica line 95 (the `Clear[...]` call), change:
```
Clear[k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022];
```
to:
```
Clear[k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022, n220, n221, n222];
```

In Mathematica line 96 (the `$Assumptions = ...` call), change:
```
$Assumptions = Element[{k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022}, Reals];
```
to:
```
$Assumptions = Element[{k20, k21, k22, m20, m21, m22, b020, b021, b022, b220, b221, b222, b420, b421, b422, z020, z021, z022, z220, z221, z222, z420, z421, z422, n020, n021, n022, n220, n221, n222}, Reals];
```

Also: the SymPy script declares `N220, N221, N222` symbols (line 192) but does not destructure them into `Nbar2, aN2, bN2` (lines 225 already does this). The Mathematica script declares `b220, ..., z220, ...` but does NOT currently declare `n220, n221, n222`. The above additions handle that. After the additions, the Mathematica script must also compute `{nbar2, aN2, bN2} = groupedParts[n220, n221, n222];` before the new assertions. Insert this line immediately before the new `nbar02, aN02, bN02` line. The SymPy script already has the analogous `Nbar2, aN2, bN2 = grouped_parts(N220, N221, N222)` at its line 225 (which currently is left unused), so no addition is needed on the SymPy side.

**Verification command:**
After Codex applies, the verifier runs `redteam exec-sympy 023` and `redteam exec-mathematica 023`. The new check names ("Nbar0 + Nbar2 additivity", "aN0 + aN2 additivity", "bN0 + bN2 additivity") should appear in both saved outputs, each evaluating to 0, and the scripts should exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- summary: Replaced the tautological N-bundle formula assertions with grouped-parts additivity checks and added the missing Mathematica N2 lane symbols/destructuring.
- deviation: none

## F2 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:307-314`

**Issue:**
The assertion `expect_zero("P0 normalization target", P0 - N0/D0)` at lines 311-314 is a duplicate of the assertion at line 290 (also `expect_zero("P0 - N0/D0", P0 - N0/D0)`). It is misnamed: it does not test any normalization target. The intermediate `P0_target = sp.solve(sp.Eq(mhat**2 * P0, 54 * G * c_s**5 / (5 * a**5 * c**5)), N0)[0]` at line 309 is computed and then never used.

**Required change:**
Replace lines 307-314 with a substantive check that actually exercises the normalization target.

Before (sympy lines 307-314):
```
    subbanner("III.3 — Universal normalization product")
    P0_target = sp.simplify(sp.solve(sp.Eq(mhat**2 * P0, 54 * G * c_s**5 / (5 * a**5 * c**5)), N0)[0])
    # This is the target on N0, but the more natural statement is on P0.
    expect_zero(
        "P0 normalization target",
        P0 - N0/D0,
    )
```

After:
```
    subbanner("III.3 — Universal normalization product")
    N0_target = sp.simplify(sp.solve(sp.Eq(mhat**2 * (N0/D0), 54 * G * c_s**5 / (5 * a**5 * c**5)), N0)[0])
    # Substantive check: N0_target, when substituted into mhat^2 * P0, must
    # reproduce the universal normalization 54 G c_s^5 / (5 a^5 c^5).
    expect_zero(
        "N0_target reproduces universal normalization",
        (mhat**2 * (N0_target/D0)) - 54 * G * c_s**5 / (5 * a**5 * c**5),
    )
```

This (a) renames the unused-and-misleading `P0_target` symbol to the accurate `N0_target` (the solver is solving for N0), and (b) replaces the duplicate `P0 - N0/D0` check with a check that actually substitutes `N0_target` into `mhat^2 * P0` and verifies the universal value drops out. Substituting concrete `D0=2`, this yields `mhat^2 * (54·G·c_s^5/(5·a^5·c^5))·D0·mhat^(-2)·D0^(-1) - 54·G·c_s^5/(5·a^5·c^5)` — let me check: `N0_target = 54·G·c_s^5·D0/(5·a^5·c^5·mhat^2)`, then `(N0_target/D0) = 54·G·c_s^5/(5·a^5·c^5·mhat^2)`, then `mhat^2 · (N0_target/D0) = 54·G·c_s^5/(5·a^5·c^5)`, then the residual is 0. Correct.

**Verification command:**
After Codex applies, `redteam exec-sympy 023` saved output no longer contains a `P0 normalization target = 0` line that duplicates the earlier `P0 - N0/D0 = 0` line; instead it contains a new line `N0_target reproduces universal normalization = 0`. Script exits 0.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
- summary: Replaced the duplicate P0 normalization assertion with an N0_target substitution check against the universal normalization.
- deviation: none

## F3 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py:305-306`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:157-158`

**Issue:**
`P2.subs(N2, N2_target)` is 0 by construction because `N2_target` was solved from `P2 == 0`. Same for `P4` under both substitutions. The checks verify solver correctness, not physics.

**Required change:**
Replace the substitute-and-check assertions with closed-form-comparison assertions.

In the SymPy script, replace lines 305-306:

Before (sympy lines 305-306):
```
    expect_zero("P2 under N2_target", P2.subs(N2, N2_target))
    expect_zero("P4 under N2_target,N4_target", P4.subs({N2: N2_target, N4: N4_target}))
```

After:
```
    # Independent closed-form derivation of N2_target and N4_target.
    # If the solver disagrees with these closed forms, the substitution check
    # would still pass tautologically; comparing to the closed form catches it.
    N2_target_closed = 2 * D2 * N0 / D0
    N4_target_closed = sp.simplify(N0 * (2 * D0 * D4 + D2**2) / D0**2)
    expect_zero("N2_target closed form", sp.simplify(N2_target - N2_target_closed))
    expect_zero("N4_target closed form", sp.simplify(N4_target - N4_target_closed))
```

In the Mathematica script, replace lines 157-158:

Before (mathematica lines 157-158):
```
expectZero["P2 under N2 target", p2 /. n2 -> n2Target];
expectZero["P4 under N2,N4 targets", p4 /. {n2 -> n2Target, n4 -> n4Target}];
```

After:
```
(* Independent closed-form derivation of n2Target and n4Target. *)
n2TargetClosed = 2*d2*n0/d0;
n4TargetClosed = FullSimplify[n0*(2*d0*d4 + d2^2)/d0^2, Assumptions -> $Assumptions];
expectZero["N2_target closed form", FullSimplify[n2Target - n2TargetClosed, Assumptions -> $Assumptions]];
expectZero["N4_target closed form", FullSimplify[n4Target - n4TargetClosed, Assumptions -> $Assumptions]];
```

The closed forms come from the saved output (sympy lines 167-173 of the .txt show `N2 = 2*D2*N0/D0` and `N4 = N0*(2*D0*D4 + D2**2)/D0**2`), so the closed forms are not invented — they are the printed values of the solver, written out explicitly so the assertion has independent content.

**Verification command:**
After Codex applies, the saved outputs contain `N2_target closed form = 0` and `N4_target closed form = 0` instead of the previous `P2 under N2 target / P4 under N2,N4 targets` lines. Scripts exit 0.

## Applied: F3

- files_changed:
  - `scripts/moving_throat_pde_stage023_full_grouped_bundle_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- summary: Replaced target-substitution checks with explicit closed-form comparisons for N2_target and N4_target.
- deviation: none

## F4 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:83-91` (Section II.0 one-port)
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl:165-178` (Section III.3 Gamma5_port)

**Issue:**
The Mathematica script reproduces the SymPy script's algebraic path step-for-step (same `Delta/S/Q/H/P` expressions, same `Series[..., {omega, 0, n}]` choreography, same `omega·D[h2,z]/h2` Hankel construction). Two engines must derive results via algebraically independent paths so that shared algebraic errors are caught.

**Required change:**
Add two independent numerical-substitution cross-checks to the Mathematica script that use distinct mechanics (numerical evaluation at concrete parameters and direct small-z Bessel expansion) from the SymPy script's symbolic series path. Do NOT remove the existing symbolic checks; add the new ones as additional assertions.

(a) Section II.0 numerical cross-check. After Mathematica line 91 (end of the symbolic Z/N one-port assertions), insert:

```
(* Independent-path verification: numerical substitution at fixed
   parameter values to verify Z_n, N_n closed forms agree with the
   rational function expansion. This breaks structural correspondence
   with the SymPy symbolic series approach. *)
zRule = {omegaU -> 2, omegaW -> 3, rMix -> 1, gU -> 1, gW -> 2};
deltaNum = deltaExpr /. zRule;
sNum = sExpr /. zRule;
qNum = qExpr /. zRule;
hNum = hExpr /. zRule;
pNum = pExpr /. zRule;
zNum0 = qNum/deltaNum;
zNum2 = (qNum*sNum - hNum*deltaNum)/deltaNum^2;
zNum4 = (qNum*(sNum^2 - deltaNum) - sNum*hNum*deltaNum)/deltaNum^3;
nNum0 = pNum^2/deltaNum^2;
nNum2 = 2*pNum*(pNum*sNum - deltaNum*gW)/deltaNum^3 /. zRule;
nNum4 = (deltaNum^2*gW^2 - 2*deltaNum*pNum^2 - 4*deltaNum*pNum*sNum*gW + 3*pNum^2*sNum^2)/deltaNum^4 /. zRule;
(* Compare to direct numerical Taylor coefficients of the rational function. *)
zNumFromSeries = SeriesCoefficient[(qExpr - hExpr*omega^2)/(deltaExpr - sExpr*omega^2 + omega^4) /. zRule, {omega, 0, #}] & /@ {0, 2, 4};
nNumFromSeries = SeriesCoefficient[((pExpr - gW*omega^2)^2/(deltaExpr - sExpr*omega^2 + omega^4)^2) /. zRule, {omega, 0, #}] & /@ {0, 2, 4};
expectZero["Z0 numerical cross-check", zNum0 - zNumFromSeries[[1]]];
expectZero["Z2 numerical cross-check", zNum2 - zNumFromSeries[[2]]];
expectZero["Z4 numerical cross-check", zNum4 - zNumFromSeries[[3]]];
expectZero["N0 numerical cross-check", nNum0 - nNumFromSeries[[1]]];
expectZero["N2 numerical cross-check", nNum2 - nNumFromSeries[[2]]];
expectZero["N4 numerical cross-check", nNum4 - nNumFromSeries[[3]]];
```

This evaluates Z_n, N_n at concrete parameters by two independent mechanics: (1) closed-form plug-in, and (2) `SeriesCoefficient` on the rational function. If the closed forms were wrong, the two would disagree.

(b) Section III.3 Bessel small-z cross-check. After Mathematica line 175 (`expectZero["Stage-5 Gamma5_port anchor", ...]`), insert:

```
(* Independent-path verification: compute the 5th-order coefficient of
   Y2/Y2_static by direct small-z Taylor expansion of j2 + I y2 (and its
   derivative), bypassing the omega*D[h2,z]/h2 ratio path used above. *)
h2SmallZ = Normal[Series[j2 + I y2, {z, 0, 9}]];
h2DerivSmallZ = Normal[Series[D[j2 + I y2, z], {z, 0, 8}]];
lambdaSmallZ = Normal[Series[(z*h2DerivSmallZ/h2SmallZ) /. z -> omega*a/cS, {omega, 0, 6}]];
ySmallZ = Normal[Series[1/lambdaSmallZ, {omega, 0, 5}]];
yStatSmallZ = ySmallZ /. omega -> 0;
yHatSmallZ = Expand[ySmallZ/yStatSmallZ];
gamma5Bessel = FullSimplify[Coefficient[yHatSmallZ, omega, 5]/I, Assumptions -> $Assumptions];
expectZero["Gamma5_port via direct Bessel small-z expansion", gamma5Bessel - a^5/(27*cS^5)];
```

This computes the same 5th-order coefficient by Taylor-expanding `j2 + I·y2` directly in `z` to high enough order, then composing with `z = omega·a/cS`, rather than the closed-form `omega·D[h2,z]/h2 /. z -> omega·a/cS` path. The two paths share the input functions but apply Series at different stages — independent enough to catch a series-order or coefficient-extraction error.

**Verification command:**
After Codex applies, `redteam exec-mathematica 023` saved output contains new pass lines: `Z0 numerical cross-check`, `Z2 numerical cross-check`, `Z4 numerical cross-check`, `N0 numerical cross-check`, `N2 numerical cross-check`, `N4 numerical cross-check`, and `Gamma5_port via direct Bessel small-z expansion`. All evaluate to 0. The Mathematica script exits 0. The pre-existing symbolic checks remain present and still pass.

## Applied: F4

- files_changed:
  - `mathematica/moving_throat_pde_stage023_full_grouped_bundle_mathematica_audit.wl`
- summary: Added the required independent Mathematica numerical-substitution and direct small-z Bessel cross-checks.
- deviation: none
