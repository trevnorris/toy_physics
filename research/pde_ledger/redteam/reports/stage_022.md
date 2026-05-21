---
unit_id: 022
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 2
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 022 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage022_grouped_p2_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.txt`

## What the script claims to verify

The scripts certify a five-step algebraic bridge that takes Stage-3/4 reduced-sector data and lands on the 2.5PN grouped-quadrupole normalization product. Step I derives the omega^2 / omega^4 coefficients of the normalized response `Y_resp(omega) = D0 / D_cons(omega)` as `u2 = -D2/D0` and `u4 = (D2^2 - D0 D4)/D0^2`. Step II computes the outgoing internal prefactor `Pref(omega) = D0 N(omega) / D_cons(omega)^2` and the full branch coefficients `K0, K2, K4, Gamma5` after multiplying by an `l=2` outgoing fingerprint `Y2_out = 1 + A omega^2 + B omega^4 + i G5 omega^5`. Step III defines a grouped-real P2 forward map `(x20, x21, x22) -> (xbar, ax, bx)` and asserts the algebraic inverse, plus an isotropy collapse `xbar = xQ, ax = 0, bx = 0` when all three lane coefficients coincide. Step IV expands the one-lane Stage-4 prototype transfer factor `N_proto = (P0_proto - gW omega^2)^2 / (Delta0 - S2 omega^2 + omega^4)^2` in omega and confirms its `N0, N2, N4` coefficients, then round-trips into Stage-4 Maxwell/mixed symbols via `Delta0 = Omega_A^2 Omega_W^2 - R^2`, `S2 = Omega_A^2 + Omega_W^2`, `P0_proto = Omega_A^2 gW + R g_A`. Step V derives the `l=2` outgoing fingerprint `A_stage4 = a^2/(9 c_s^2)`, `B_stage4 = 4 a^4/(81 c_s^4)`, `G5_stage4 = a^5/(27 c_s^5)` from the spherical-Hankel-function DtN combinatorics, solves `mhat^2 P0 Gamma5_port = gamma_GR = 2G/(5c^5)` for the invariant product, and asserts the three `mhat = 1` literal targets `K0_target = 54 G c_s^5/(5 a^5 c^5)`, `K2_target = 6 G c_s^3/(5 a^3 c^5)`, `K4_target = 8 G c_s/(15 a c^5)`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1  | sympy        | 80      | `u2 + D2/D0 == 0` (from series of `D0/Dcons`) | yes |
| A2  | sympy        | 81      | `u4 - (D2^2 - D0 D4)/D0^2 == 0` | yes |
| A3  | sympy        | 118     | `P0 - N0/D0 == 0` | yes |
| A4  | sympy        | 119     | `P2 - (D0 N2 - 2 D2 N0)/D0^2 == 0` | yes |
| A5  | sympy        | 120-123 | `P4 - (D0^2 N4 - 2 D0(D2 N2 + D4 N0) + 3 D2^2 N0)/D0^3 == 0` | yes |
| A6  | sympy        | 133     | `K0 - P0 == 0` | partial (`K0 = pref[0] * Y2_out[0] = P0 * 1 = P0` by construction) |
| A7  | sympy        | 134     | `K2 - (P2 + A P0) == 0` | yes |
| A8  | sympy        | 135     | `K4 - (P4 + A P2 + B P0) == 0` | yes |
| A9  | sympy        | 136     | `Gamma5 - G5 P0 == 0` | yes |
| A10 | sympy        | 181-183 | inverse round-trip `x2N_back - x2N == 0` | yes |
| A11 | sympy        | 193-195 | isotropy `xbar_iso - xQ = ax_iso = bx_iso = 0` when all lanes equal `xQ` | yes |
| A12 | sympy        | 217-222 | `N0, N2, N4` prototype formulas | yes |
| A13 | sympy        | 250-261 | round-trip `Nseries[Nproto.subs(dict_back)].coeff(omega,k) - Nseries[Nproto].coeff(omega,k).subs(dict_back) == 0` | NO (tautological — see F2) |
| A14 | sympy        | 293-295 | Stage-4 outgoing fingerprint `A, B, G5` from spherical-Hankel `Lambda2` | yes |
| A15 | sympy        | 315-318 | `mhat=1` invariant-product target literal | yes |
| A16 | sympy        | 331-338 | `mhat=1` `K2, K4` target literals | yes |
| M1  | mathematica  | 40-41   | mirror of A1, A2 via `Series[d0/dCons]` | yes (same claim, non-independent route — F1) |
| M2  | mathematica  | 57-66   | mirror of A3-A9 via `pref = d0*nFac/dCons^2` series | yes (same claim, non-independent route — F1) |
| M3  | mathematica  | 75-80   | mirror of A10, A11 inverse / isotropy | yes (same claim, non-independent route — F1) |
| M4  | mathematica  | 94-99   | mirror of A12 prototype | yes (same claim, non-independent route — F1) |
| M5  | mathematica  | 107-109 | mirror of A13 round-trip | NO (tautological — see F2) |
| M6  | mathematica  | 134-136 | mirror of A14 outgoing fingerprint | yes (same claim, non-independent route — F1) |
| M7  | mathematica  | 147-149 | mirror of A15, A16 invariant-product literals | yes (same claim, non-independent route — F1) |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl:33-149`
- compared against `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py:65-338`

**What's wrong:**

The Mathematica audit script is a section-by-section transliteration of the SymPy audit script. Every algebraic step matches the SymPy choreography (same intermediate variables, same series-expansion approach, same coefficient extraction order, same target formulas). The only differences are syntactic — `sp.series` vs `Series[...]`, `coeff` vs `Coefficient`, `simplify` vs `FullSimplify`, snake_case vs camelCase. Three matched excerpts:

1. **Section I — `Y_resp` series and coefficient extraction.**

   SymPy (`scripts/.../stage022_..._sympy_audit.py:71-74`):
   ```
   Dcons = D0 + D2 * omega**2 + D4 * omega**4
   Yresp = sp.expand(sp.series(D0 / Dcons, omega, 0, 6).removeO())
   u2 = sp.simplify(Yresp.coeff(omega, 2))
   u4 = sp.simplify(Yresp.coeff(omega, 4))
   ```

   Mathematica (`mathematica/.../stage022_..._mathematica_audit.wl:33-36`):
   ```
   dCons = d0 + d2*omega^2 + d4*omega^4;
   yResp = Expand[Normal[Series[d0/dCons, {omega, 0, 4}]]];
   u2 = FullSimplify[Coefficient[yResp, omega, 2], Assumptions -> $Assumptions];
   u4 = FullSimplify[Coefficient[yResp, omega, 4], Assumptions -> $Assumptions];
   ```

   Identical algebraic path: build `Dcons`, take the Taylor series of `D0/Dcons` around omega=0, read off the omega^2 and omega^4 coefficients. The Mathematica version does not, for instance, verify `u2 = -D2/D0` via the inverse-relation `Y_resp * Dcons - D0 = 0` and a coefficient match on the *product* (which is a genuinely distinct route).

2. **Section IV — Stage-4 prototype `N_proto` series.**

   SymPy (lines 210-214):
   ```
   Nproto = (P0_proto - gW * omega**2) ** 2 / (Delta0 - S2 * omega**2 + omega**4) ** 2
   Nseries = sp.expand(sp.series(Nproto, omega, 0, 6).removeO())
   N0 = sp.simplify(Nseries.coeff(omega, 0))
   N2 = sp.simplify(Nseries.coeff(omega, 2))
   N4 = sp.simplify(Nseries.coeff(omega, 4))
   ```

   Mathematica (lines 87-91):
   ```
   nProto = (p0proto - gW*omega^2)^2/(delta0 - s2*omega^2 + omega^4)^2;
   nSeries = Expand[Normal[Series[nProto, {omega, 0, 4}]]];
   n0Proto = FullSimplify[Coefficient[nSeries, omega, 0], Assumptions -> $Assumptions];
   n2Proto = FullSimplify[Coefficient[nSeries, omega, 2], Assumptions -> $Assumptions];
   n4Proto = FullSimplify[Coefficient[nSeries, omega, 4], Assumptions -> $Assumptions];
   ```

   Same `N_proto` factored form, same coefficient-extraction order. An independent Mathematica derivation would, for example, factor `N_proto = (P0_proto - gW omega^2)^2 * (1 / D(omega))^2`, compute `(1/D)` via `InverseSeries`/polynomial-inversion, and convolve the two factors; or use `SeriesData` operations on the factored form directly. The current code matches SymPy's pattern.

3. **Section V — `Lambda2` from spherical Hankel and coefficient match.**

   SymPy (lines 280-291):
   ```
   j2 = (sp.Rational(3, 1) / z**3 - sp.Rational(1, 1) / z) * sp.sin(z) - 3 * sp.cos(z) / z**2
   y2 = -((sp.Rational(3, 1) / z**3 - sp.Rational(1, 1) / z) * sp.cos(z) + 3 * sp.sin(z) / z**2)
   h2 = sp.simplify(j2 + I * y2)
   Lambda2 = sp.simplify((omega := sp.symbols("omega", real=True)) * sp.diff(h2, z) / h2).subs(z, omega * a / c_s)
   Lambda2_series = sp.series(Lambda2, omega, 0, 7).removeO()
   Y2 = sp.simplify(sp.series(1 / Lambda2_series, omega, 0, 6).removeO())
   ```

   Mathematica (lines 121-126):
   ```
   j2 = ((3/z^3) - 1/z) Sin[z] - 3 Cos[z]/z^2;
   y2 = -((3/z^3) - 1/z) Cos[z] - 3 Sin[z]/z^2;
   h2 = FullSimplify[j2 + I y2, Assumptions -> $Assumptions];
   lambda2 = FullSimplify[(omega D[h2, z]/h2) /. z -> omega a/cS, Assumptions -> $Assumptions];
   lambda2Series = Normal[Series[lambda2, {omega, 0, 6}]];
   y2Resp = Normal[Series[1/lambda2Series, {omega, 0, 5}]] // FullSimplify;
   ```

   Same explicit `j2`, `y2`, `h2 = j2 + i y2`, same `Lambda2 = omega * h2'(z)/h2(z)` then substitute, same series-then-invert recipe. Mathematica has a built-in `SphericalHankelH1[n, z]` which gives an independent route to `h2`; the script does not use it. This makes the Mathematica check a transliteration of SymPy's explicit polynomial-rational form rather than an independent derivation from a standard-library special function.

**Why this matters:**

The two-engine policy exists so that a SymPy bug or convention quirk does not silently propagate into "verified". When both scripts walk the same algebraic path, a bug in the shared *idea* (e.g. wrong choice of branch, wrong normalization of `Y2_out`, sign error in `h2 = j2 + i y2` convention) cannot be caught by the cross-check. Both engines pass, both engines are wrong.

This stage is the load-bearing bridge from the Stage-4 reduced-sector prototype to the 2.5PN normalization target, so the second-engine policy applies at full strength.

**Required change:**

Rewrite the Mathematica audit script so that each section computes the LHS of its `expectZero` calls via a route that is *not* "build the same intermediate, expand the same series, read the same coefficient". Concretely:

- Section I — replace the `Series[d0/dCons]` route with the inverse-relation route: form `yRespCand = 1 + u2Sym * omega^2 + u4Sym * omega^4` with symbolic unknowns `u2Sym, u4Sym`, multiply by `dCons`, expand, and use `SolveAlways` (or `Solve` on the coefficient list) to determine `u2Sym, u4Sym` from `yRespCand * dCons - d0 = 0 mod omega^6`. Then assert `u2Sym + d2/d0 == 0` etc. from the solution.

- Section II — same inverse-relation strategy for `Pref = d0 * nFac / dCons^2`. Form `prefCand` as an unknown polynomial of degree 4, set `prefCand * dCons^2 - d0 * nFac = 0 mod omega^6`, solve for the coefficients, compare to the target formulas. For the branch step, replace `pref * y2Out` series with a direct polynomial expansion in omega (a single Expand, no Series) and read coefficients from the polynomial. Adjust `expectZero` LHS to be `(K0_via_polynomial - target) == 0`.

- Section IV — replace `Series[nProto]` with the inverse-relation route: assume `nSeriesCand = n0Cand + n2Cand omega^2 + n4Cand omega^4` and require `nSeriesCand * (delta0 - s2 omega^2 + omega^4)^2 - (p0proto - gW omega^2)^2 = 0 mod omega^6`. Solve, compare to the target formulas. Also, drop the IV.2 round-trip (it is tautological — see F2).

- Section V — replace explicit `j2, y2, h2` with `SphericalHankelH1[2, z]` (Mathematica's built-in) and verify `Lambda2 = omega * D[h, z]/h` matches; then proceed with the series. This makes Mathematica's derivation of `A_stage4, B_stage4, G5_stage4` independent of the explicit polynomial-rational form used by SymPy.

The targets (`u2 = -d2/d0`, `P0 = N0/D0`, `N0 = P0^2/Delta0^2`, `A_stage4 = a^2/(9 c_s^2)`, etc.) are the *claim*; do not change them. Change the *route* used to compute the LHS.

**Verification:**

The verifier runs `redteam exec-mathematica 022`. The new `.wl` file must (a) exit code 0, (b) contain `SolveAlways` or `Solve` calls on coefficient equations in Sections I, II, IV, (c) contain `SphericalHankelH1[2, ...]` in Section V, (d) keep all existing `expectZero` target literals unchanged on the RHS of each check.

### F2 — tautological_check

**Severity:** low

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage022_grouped_p2_sympy_audit.py:245-261`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage022_grouped_p2_normalization_bridge_mathematica_audit.wl:104-109`

**What's wrong:**

Section IV.2 in both engines runs a "round-trip" check that is algebraically guaranteed by construction. The SymPy form is:

```
Nproto_stage4 = (
    (dict_back[P0_proto] - gW * omega**2) ** 2
    / (dict_back[Delta0] - dict_back[S2] * omega**2 + omega**4) ** 2
)
Nseries_stage4 = sp.expand(sp.series(Nproto_stage4, omega, 0, 6).removeO())
expect_zero(
    "N0 round-trip into Stage-4 symbols",
    Nseries_stage4.coeff(omega, 0) - sp.simplify(N0.subs(dict_back)),
)
```

`Nproto_stage4` is defined by literally substituting `dict_back` into the same template that produced `Nproto`. The check then claims `Series[Nproto.subs(dict_back)].coeff(omega, k) == Series[Nproto].coeff(omega, k).subs(dict_back)`. Since `dict_back` introduces no `omega`-dependence (`Delta0 = Omega_A^2 Omega_W^2 - R^2`, `S2 = Omega_A^2 + Omega_W^2`, `P0_proto = Omega_A^2 gW + R g_A` are all `omega`-independent), Taylor expansion in `omega` commutes with substitution into those coefficients. The residual is therefore identically zero — the check cannot fail unless SymPy itself has a bug in `series` or `subs`, neither of which is the physics this audit is supposed to certify.

The Mathematica mirror (lines 104-109) has the same structure (`nStage4 = Series[((p0Back - gW*omega^2)^2)/(...)^2]`, then assert `Coefficient[nStage4, omega, k] - (nKProto /. {delta0 -> deltaBack, s2 -> s2Back, p0proto -> p0Back}) == 0`) and the same tautology.

**Why this matters:**

A passing check that cannot fail is dead weight in a verification ledger. It contributes nothing to confidence in the dictionary mapping itself (does `P0_proto = Omega_A^2 gW + R g_A` actually reproduce the Stage-4 forward result?), because both LHS and RHS are built from the same template. If the dictionary entries were wrong, both LHS and RHS would be equally wrong, and the residual would still be zero.

The actual claim that needs verifying — that the Stage-4 Maxwell/mixed expressions for `Delta0, S2, P0_proto` are the ones that come out of the Stage-3/4 reduced-sector derivation — is not in scope for this script (that's an upstream-unit claim). So the right answer is to delete the redundant check, not to expand it.

**Required change:**

Delete the Section IV.2 round-trip block in both engines. The dictionary print-outs at SymPy lines 238-243 and Mathematica lines 101-103 can remain (they are documentation, not assertions), but the three `expect_zero`/`expectZero` calls and the construction of `Nseries_stage4` / `nStage4` should be removed.

SymPy edit: delete lines 245-261 (the `Nproto_stage4` and `Nseries_stage4` constructions, plus the three `expect_zero` round-trip calls).

Mathematica edit: delete lines 104-109 (the `nStage4` construction plus the three `expectZero` round-trip calls).

**Verification:**

The verifier runs `redteam exec-sympy 022` and `redteam exec-mathematica 022`. Both must exit 0. The output files must no longer contain lines `N0 round-trip into Stage-4 symbols`, `N2 round-trip into Stage-4 symbols`, `N4 round-trip into Stage-4 symbols`. The remaining assertions (IV.1 prototype formulas, V fingerprint coefficients, V invariant target) must still appear and still print residual `0`.

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. See F1 for three matched code excerpts. The script uses the same algebraic recipe in each of the five sections: build the same intermediate (`Dcons`, `Pref`, `Nproto`, explicit `j2/y2/h2`), take the same Taylor series in `omega` around 0, extract the same coefficients, compare to the same target formula. There is no section in which Mathematica reaches the claim by a path that diverges from SymPy's algebra.

## Engine cross-check

Both engines pass all assertions (SymPy output `EXIT_CODE: 0` with all `expect_zero` residuals reported as `0`; Mathematica output `EXIT_CODE: 0` with all `PASS` lines). The two engines report identical final forms for `Gamma5_port = a^5/(27*c_s^5)` and `mhat^2 * P0 = 54 G c_s^5 / (5 a^5 c^5 mhat^2)`. They agree on result; the disagreement issue is that they agree via the *same* derivation, which is what F1 flags.

## Verdict justification

`findings` (two findings: one `mathematica_transliteration`, one `tautological_check`). The SymPy script holds up under attack on its own terms: the series-coefficient checks in Sections I-II non-trivially verify the closed-form expressions for `u2, u4, P0, P2, P4, K0..K4, Gamma5` (a typo in any target formula would flip the residual sign and fail the assertion); the inverse-map round-trip in Section III is not tautological because the inverse formulas `x20 = xbar + 4 ax`, etc., are independent algebraic statements that could fail if a coefficient were wrong; the Stage-4 prototype expansion in Section IV.1 non-trivially verifies the `N0, N2, N4` formulas; the spherical-Hankel derivation in Section V genuinely verifies `A_stage4 = a^2/(9 c_s^2)`, `B_stage4 = 4 a^4/(81 c_s^4)`, `G5_stage4 = a^5/(27 c_s^5)` from the explicit `h2 = j2 + i y2` form, and the literal mhat=1 target checks in V.2/V.3 verify that the constant arithmetic `(2/5) / (1/27) = 54/5`, `(54/5)(1/9) = 6/5`, `(54/5)(4/81) = 8/15` is correct. Adversarial attacks I tried that failed: (i) checking whether positivity of `z` in V causes hidden simplification issues — the series is in `omega` around 0, not in `z`, and the substitution `z -> omega*a/c_s` happens before series so no `sqrt(z^2)` artifacts; (ii) checking whether `h2 = j2 + i y2` has the right small-`omega` parity to produce a purely-imaginary `omega^5` coefficient — yes, `j2` is even in `z`, `y2` is odd in `z`, so `Lambda2 = omega h2'/h2` has even-in-omega real part and odd-in-omega imaginary part, and `Y2 = 1/Lambda2` inherits that parity pattern, putting `G5` at `omega^5/I`; (iii) checking `nonzero` vs `positive` assumption mismatches between `D0, D2, D4` (declared nonzero, real) and the simplification path — no positivity is assumed in I-IV, and Section V's positive declarations are consistent with the stated physical scope. The two real defects are structural: the Mathematica engine walks the same algebraic path as SymPy (F1), and the IV.2 round-trip is a substitution-commutes-with-series identity that cannot fail (F2).

## Self-test notes

- I checked that F1's required changes (`SolveAlways` on coefficient equations, `SphericalHankelH1` for `h2`) produce LHS expressions in the *same* variables already used by the script, so the existing `expectZero` target literals continue to match. The Mathematica `SolveAlways` returns rules `{u2Sym -> -d2/d0, u4Sym -> (d2^2 - d0 d4)/d0^2}`, which after `/.` substitution into `u2Sym + d2/d0` gives `0` — trap (1) variable-independence: `Solve` produces rules whose RHS still depends on `d0, d2, d4` (since those are the symbolic coefficients of `dCons`), so the target subtraction is non-trivial. Trap (2) parity: doesn't apply (no unbounded integrals). Trap (3) trivial-case: for `nFac = 1` (i.e. `n2 = n4 = 0, n0 = 1`), the inverse-relation route gives `P0 = 1/d0`, which matches `N0/D0 = 1/D0`; for `nFac = 0` (all `n0 = n2 = n4 = 0`), both routes return `0`, consistent.
- I checked that F2's deletion doesn't strand any later assertion. The IV.2 block defines `Nproto_stage4` and `Nseries_stage4`, which are not used after line 261 in SymPy or line 109 in Mathematica. Section V independently derives `A_stage4, B_stage4, G5_stage4` from the spherical Hankel function, not from `Nseries_stage4`. So the deletion is safe.
- Path specifications for F1 use the existing Mathematica file path — no missing-script path-guessing risk here, since both engine scripts already exist.
