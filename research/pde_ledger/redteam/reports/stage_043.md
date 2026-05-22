---
unit_id: 043
batch: III.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-22T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 043 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage043_support_direction_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage043_support_direction_mathematica_audit.txt`

## What the script claims to verify

The audit claims that integrating out the split-U doublet from a continuum extension with one symmetry-allowed `U/phi` bilinear (a) produces an effective support vector `y = gB v + gU gS D_U v` with a closed-form direction factor `R_phi = (1 + sigma0/(1+delta_U)) / (1 + sigma0)`; (b) that the splitting invariant `D_phi = kappa0 y1 - kappa1 y0` reduces to the expected single product; (c) that the U-overlap contraction `v.D_U.v` under the overlap relation `kappa1^2 = (2/11) sigma` gives the split blocking factor `(1 - (2/11) delta_U/(1+delta_U))`; (d) that the physical support baseline reduces to `8/(pi^2) * cB^2 * (1+c_etaU c_Uphi/(KU cB))^2 / [Keta_eff Kphi_eff (1-eps_eta)(1-eps_phi^split)]`; (e) that the support and mixed vectors are aligned iff `gB gR = gW gS`, with the closed-form mismatch `delta_U(rho0-sigma0) / [(1+delta_U)(1+rho0)(1+sigma0)]`.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 74 | `expect_zero(R_phi - Rphi_expected)` | yes |
| A2 | sympy | 81 | `expect_zero(D_phi - Dphi_expected)` | yes |
| A3 | sympy | 93 | `expect_zero(SU_sub - SU_expected)` | yes |
| A4 | sympy | 101 | `expect_zero(Aphi_eff - Aphi_eff_expected)` | no (tautological) |
| A5 | sympy | 120 | `expect_zero(Msupp_cont_eval - Msupp_expected)` | no (tautological) |
| A6 | sympy | 130 | `expect_zero(Dzphi - Dzphi_expected)` | yes |
| A7 | sympy | 134 | `expect_zero(tracking via gS=gB gR/gW)` | yes |
| A8 | sympy | 142 | `expect_zero(mismatch - mismatch_expected)` | yes |
| B1 | mathematica | 58 | `expectZero[rPhi - rPhiExpected]` | yes |
| B2 | mathematica | 59 | `expectZero[dPhi - dPhiExpected]` | yes |
| B3 | mathematica | 79 | `expectZero[sUSub - sUExpected]` | yes |
| B4 | mathematica | 80 | `expectZero[aPhiEff - aPhiEffExpected]` | no (tautological) |
| B5 | mathematica | 98 | `expectZero[mSuppContEval - mSuppExpected]` | no (tautological) |
| B6 | mathematica | 113 | `expectZero[dPhiZ - dPhiZExpected]` | yes |
| B7 | mathematica | 114 | `expectZero[tracking via gS->gB gR/gW]` | yes |
| B8 | mathematica | 125 | `expectZero[mismatch - mismatchExpected]` | yes |

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:43-125`

**What's wrong:**

The `.wl` script is a line-by-line port of the `.py` script with identifier renaming (snake_case → camelCase) and no independent derivation path. Compare the choreography:

SymPy (lines 59–81):
```
DU = sp.diag(1 / KU, 1 / (KU * (1 + delta_U)))
v = sp.Matrix([kappa0, kappa1])
y = sp.simplify(gB * v + gU * gS * DU * v)
sigma0 = sp.simplify(gU * gS / (KU * gB))
Rphi = sp.simplify((y1 / y0) / (kappa1 / kappa0))
Rphi_expected = sp.simplify((1 + sigma0 / (1 + delta_U)) / (1 + sigma0))
Dphi = sp.factor(sp.expand(kappa0 * y1 - kappa1 * y0))
Dphi_expected = sp.simplify(-kappa0 * kappa1 * gB * sigma0 * delta_U / (1 + delta_U))
```

Mathematica (lines 43–53):
```
dU = DiagonalMatrix[{1/kU, 1/(kU (1 + deltaU))}];
v = {kappa0, kappa1};
y = FullSimplify[gB v + gU gS (dU.v), ...];
sigma0 = FullSimplify[gU gS/(kU gB), ...];
rPhi = FullSimplify[(y1/y0)/(kappa1/kappa0), ...];
rPhiExpected = FullSimplify[(1 + sigma0/(1 + deltaU))/(1 + sigma0), ...];
dPhi = FullSimplify[kappa0 y1 - kappa1 y0, ...];
dPhiExpected = FullSimplify[-kappa0 kappa1 gB sigma0 deltaU/(1 + deltaU), ...];
```

Same matrix definition, same vector, same `y` formula, same `sigma0` formula, same `R_phi` ratio, same expected closed form, same `D_phi` definition, same expected closed form — all in the same order. The same pattern continues through sections 2–5: the substitution choices (`kappa1^2 → (2/11) sigma`, `kappa0^2 → 8/Pi^2`), the construction of `Aphi_eff = Kphi_eff - cUphi^2 * SU_expected`, and the `Msupp_cont` quotient with the `mu_eta * mu_phi` cancellation all appear identically in both engines.

The second-engine policy requires the Mathematica script to derive the result from physical premises independently — e.g., starting from a different parameterization, building `y` from a linear solve / nullspace computation, or computing `R_phi` as a tangent ratio rather than as `(y1/y0)/(kappa1/kappa0)`. None of that is present.

**Why this matters:**

A transliterated second engine cannot catch algebra mistakes in the first — both will reproduce the same error in lock-step. Cross-engine agreement under transliteration is no stronger than running the same script twice.

**Required change:**

Re-cast the Mathematica script so the same five results are arrived at by an independent computational path. Concrete suggestions (Codex picks one consistent route):
- Solve for `y` symbolically as `LinearSolve` on `(gB^{-1} I - gU gS gB^{-1} DU) y = v` form, or equivalently express `y` as the residue of a 1×1 effective Green's function `gB(1 + sigma0 D_U_dimensionless)`, and read off `R_phi` from `y[[2]]/(kappa1) ÷ y[[1]]/kappa0` AFTER simplification — not via the same intermediate `Rphi = (y1/y0)/(kappa1/kappa0)` expression.
- Derive `D_phi` as the determinant of the `2×2` matrix `[[y0, y1],[kappa0, kappa1]]` via `Det[...]` rather than by writing the expression `kappa0 y1 - kappa1 y0` by hand.
- Derive the overlap `v.D_U.v` AND the closed `(1 - (2/11) delta_U/(1+delta_U))` factor by expanding `(v.D_U.v) / (sigma/kU)` and evaluating the limit `deltaU → 0` (should give 1) and `deltaU → infinity` (should give 9/11), then closed-form-fit; verify both endpoints agree.
- Compute `M_supp` by first computing the propagator quotient symbolically with no `kappa0^2 → 8/Pi^2` hand-substitution, then verifying the prefactor at `kappa0^2 → 8/Pi^2` separately. The current `Msupp_cont` cancels `mu_eta * mu_phi` against `mu_eta * mu_phi` — that step adds no content.
- Compute the mismatch via `Series[rPhiExpected - rU, {deltaU, 0, 2}]` and confirm the leading nonzero coefficient matches `(rho0 - sigma0)/((1+sigma0)(1+rho0))`, then extend to exact closed form by `Together[...]`.

**Verification:**

After Codex applies, the new `.wl` script should NOT contain the lines `kappa0 y1 - kappa1 y0` (replaced by `Det`), `kappa0^2 -> 8/Pi^2` as a single-step substitution into the constructed quotient (replaced by per-prefactor verification), and the variable choreography should not appear as a Mathematica image of the Python lines. The output `.txt` should still report `PASS` for each of the five claims, with at least one new diagnostic line printed (e.g., a limit-check or a series-expansion-check) that has no SymPy counterpart.

### F2 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:96-101`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:69-80`

**What's wrong:**

The `A_phi^(eff) - expected` check is algebraically trivial because both sides are constructed from the same `SU_expected`.

SymPy (lines 96–101):
```
eps_phi_split = sp.simplify(eps_phi * (1 - sp.Rational(2, 11) * delta_U / (1 + delta_U)))
Aphi_eff = sp.simplify(Kphi_eff - cUphi**2 * SU_expected)
Aphi_eff_expected = sp.simplify(Kphi_eff * (1 - eps_phi_split)).subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))
...
expect_zero("A_phi^(eff) - expected", Aphi_eff - Aphi_eff_expected)
```

Substituting `eps_phi = cUphi^2 * sigma/(KU * Kphi_eff)` into `Kphi_eff * (1 - eps_phi_split) = Kphi_eff - eps_phi * Kphi_eff * (1 - (2/11) delta_U/(1+delta_U))` gives `Kphi_eff - cUphi^2 * sigma/KU * (1 - (2/11) delta_U/(1+delta_U)) = Kphi_eff - cUphi^2 * SU_expected`. That is precisely the definition of `Aphi_eff` on line 97. So `Aphi_eff - Aphi_eff_expected = 0` by algebraic substitution alone — independent of whether the *physical* identification `eps_phi = cUphi^2 sigma/(KU Kphi_eff)` is the right pole-shift, and independent of whether the `(1 - (2/11) delta_U/(1+delta_U))` factor is the right overlap reduction.

The same shape appears in the Mathematica script at lines 69–80 (`aPhiEff = kPhiEff - cUphi^2 * sUExpected`, `aPhiEffExpected = kPhiEff*(1-epsPhiSplit) /. epsPhi -> cUphi^2 sigma/(kU kPhiEff)`).

Both engines re-confirm that `Kphi_eff - cUphi^2 * S = Kphi_eff - cUphi^2 * S`. This is not a check; it's a re-statement.

**Why this matters:**

This is the claim that links the overlap reduction `(1 - (2/11) delta_U/(1+delta_U))` to the *physical* split blocking ratio `eps_phi^(split)`. Without an assertion that genuinely exercises that link, a sign error or factor of 2 in `eps_phi^(split)` could not be caught here — both sides would still match because both are constructed from the same factor.

**Required change:**

Add a genuine assertion that anchors `eps_phi` *separately from* the construction of `Aphi_eff`. Specifically, in both `.py` and `.wl`:

In the SymPy script (`/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py`), between current lines 101 and 102 (before subbanner 26.3), insert an independent check that the *minimal-overlap* limit `delta_U = 0` reduces `Aphi_eff` to `Kphi_eff * (1 - eps_phi)` with `eps_phi = cUphi^2 sigma/(KU Kphi_eff)`:

```python
Aphi_eff_min = sp.simplify(Aphi_eff.subs(delta_U, 0))
Aphi_eff_min_expected = sp.simplify(Kphi_eff - cUphi**2 * sigma / KU)
expect_zero("A_phi^(eff) at delta_U=0 (minimal)", Aphi_eff_min - Aphi_eff_min_expected)
```

Then add the genuinely non-tautological identity that the *ratio* `(Kphi_eff - Aphi_eff) / (Kphi_eff - Aphi_eff_min)` equals exactly `(1 - (2/11) delta_U/(1+delta_U))`:

```python
overlap_ratio = sp.simplify((Kphi_eff - Aphi_eff) / (Kphi_eff - Aphi_eff_min))
overlap_ratio_expected = sp.simplify(1 - sp.Rational(2, 11) * delta_U / (1 + delta_U))
expect_zero("split-vs-minimal overlap ratio", overlap_ratio - overlap_ratio_expected)
```

This check is non-tautological because `Aphi_eff` and `Aphi_eff_min` are each set up from the contracted `SU` quantity, and the ratio specifically tests that `kappa1^2 / kappa0^2 = 2/9` (equivalently `kappa1^2/sigma = 2/11`) produces the claimed `2/11` weight on the `delta_U` correction. A sign or factor error in `SU_expected` would make this ratio fail.

Apply the same shape to the Mathematica script after line 80, mirroring the limit `deltaU -> 0` and the ratio check, but using the independent path called out in F1 (e.g., evaluate `aPhiEffMin = Limit[aPhiEff, deltaU -> 0]` rather than `aPhiEff /. deltaU -> 0`).

**Verification:**

After Codex applies, the SymPy output should contain a new line `A_phi^(eff) at delta_U=0 (minimal) = 0` and a new line `split-vs-minimal overlap ratio = 0`. The Mathematica output should contain `PASS: A_phi^(eff) at deltaU=0 (minimal)` and `PASS: split-vs-minimal overlap ratio`. Verifier confirms both new checks appear AND scripts exit 0.

### F3 — tautological_check

**Severity:** high
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage043_support_direction_sympy_audit.py:107-120`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage043_support_direction_mathematica_audit.wl:84-98`

**What's wrong:**

The `M_supp - expected` check is structurally tautological: `Msupp_cont` is set up as `(numer)/(denom)` and `Msupp_expected` is set up as the same `(numer)/(denom)` post-cancellation of `mu_eta * mu_phi`.

SymPy (lines 107–120):
```python
Msupp_cont = sp.simplify(
    (kappa0**2 * cB**2 * (1 + ceU * cUphi / (KU * cB))**2 / (mu_eta * mu_phi))
    / ((Keta_eff * (1 - eps_eta) / mu_eta) * (Kphi_eff * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))) / mu_phi))
)
Msupp_expected = sp.simplify(
    (sp.Rational(8, 1) / sp.pi**2)
    * (cB**2 / (Keta_eff * Kphi_eff))
    * (1 + ceU * cUphi / (KU * cB))**2
    / ((1 - eps_eta) * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))))
)
Msupp_cont_eval = sp.simplify(Msupp_cont.subs(kappa0**2, sp.Rational(8, 1) / sp.pi**2))
expect_zero("M_supp - expected", Msupp_cont_eval - Msupp_expected)
```

The `mu_eta * mu_phi` in the numerator of `Msupp_cont` cancels the `mu_eta * mu_phi` in its denominator algebraically. After that cancellation and the hand-substitution `kappa0^2 = 8/pi^2`, the residual is identical to `Msupp_expected` term-for-term: `(8/pi^2)(cB^2/(Keta_eff Kphi_eff))(1+ceU cUphi/(KU cB))^2 / [(1-eps_eta)(1-eps_phi_split_phys)]`. The check thus verifies `(a/b)*(b/a) * X == X`. Nothing physical.

Critically, the value `kappa0^2 = 8/pi^2` (which the docstring claim labels the "Z_phi-driven baseline") is *substituted by hand* into `Msupp_cont`, not derived from any earlier result. Whatever number is plugged in, both sides absorb it identically.

The same shape appears in the Mathematica script lines 84–98 (`mSuppCont` over `mSuppExpected`, with `kappa0^2 -> 8/Pi^2` substitution).

**Why this matters:**

This is claim #5 ("the actual physical support baseline"). The only piece that touches new physics is the value `kappa0^2 = 8/pi^2` — but it is treated as a given, not derived, and not used in any non-cancelling combination. A wrong baseline factor would propagate unnoticed.

**Required change:**

In both engines, split the `M_supp` claim into two genuinely independent assertions:

(a) The *mu-independence* / cancellation structure: verify that `Msupp_cont` written WITHOUT the `kappa0^2` substitution is independent of `mu_eta` and `mu_phi` (currently obvious by construction but absent as a check). In SymPy:

```python
Msupp_mu_independent = sp.simplify(sp.diff(Msupp_cont, mu_eta)) + sp.simplify(sp.diff(Msupp_cont, mu_phi))
expect_zero("M_supp independent of mu_eta, mu_phi", Msupp_mu_independent)
```

(b) The *baseline-value identification*: assert separately that the substitution `kappa0^2 = 8/pi^2` produces the leading `8/pi^2` numerical prefactor — by writing `Msupp_expected` symbolically with a free baseline `B`, then asserting `B = 8/pi^2` only at the final comparison:

```python
B = sp.symbols("B_baseline", positive=True, real=True)
Msupp_expected_sym = sp.simplify(
    B
    * (cB**2 / (Keta_eff * Kphi_eff))
    * (1 + ceU * cUphi / (KU * cB))**2
    / ((1 - eps_eta) * (1 - eps_phi_split.subs(eps_phi, cUphi**2 * sigma / (KU * Kphi_eff))))
)
Msupp_cont_in_B = sp.simplify(Msupp_cont.subs(kappa0**2, B))
expect_zero("M_supp structural form (free baseline)", Msupp_cont_in_B - Msupp_expected_sym)
expect_zero("baseline value B = 8/pi^2", (B - sp.Rational(8, 1) / sp.pi**2).subs(B, sp.Rational(8, 1) / sp.pi**2))
```

Now the structural form check uses a *symbolic* `B`, so it actually exercises whether `Msupp_cont` reduces to the claimed product structure for an arbitrary baseline; the second check separately records the hard-coded baseline value. If the `8/pi^2` is wrong, the second assertion is the place to fix it; the structural assertion remains valid.

Apply the same split to the Mathematica script lines 84–98, using `D[mSuppCont, muEta]` and `D[mSuppCont, muPhi]` for the mu-independence check (verified by `FullSimplify[#, $Assumptions]&`), and a free symbol `bBaseline` for the structural form. Reach the structural form via an independent path (per F1) — for instance, build `mSuppExpected_sym` from the residue of the 2×2 effective propagator at the support pole rather than copying the SymPy quotient form.

**Verification:**

After Codex applies, the SymPy output should contain:
- `M_supp independent of mu_eta, mu_phi = 0`
- `M_supp structural form (free baseline) = 0`
- `baseline value B = 8/pi^2 = 0`

The Mathematica output should contain the same three diagnostic lines (named identically). All five claims should still pass.

## Independent-derivation check (Mathematica)

The Mathematica script is structurally a line-by-line port of the SymPy script with the variable-renaming convention `delta_U → deltaU`, `K_U → kU`, `c_etaU → cEtaU`, etc. Every intermediate expression and every assertion appears in the same order with the same algebraic form. There is no independent derivation path. This is filed as F1 (`mathematica_transliteration`).

## Engine cross-check

Both engines produce identical residuals (`0`) for all eight assertions. The pre-residual displays also match:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `R_phi` | `(K_U δ_U g_B + K_U g_B + g_S g_U) / [(δ_U+1)(K_U g_B + g_S g_U)]` | `1 - (deltaU gS gU)/((1+deltaU)(gS gU + gB kU))` |
| `D_phi` | `-δ_U g_S g_U κ_0 κ_1 / (K_U(δ_U+1))` | `-(deltaU gS gU kappa0 kappa1)/(kU + deltaU kU)` |
| `v.D_U.v` | `σ(9 δ_U + 11) / [11 K_U (δ_U+1)]` | `(sigma + (9 deltaU sigma)/11)/(kU + deltaU kU)` |
| `A_phi^(eff)` | `(11 K_U K_φ_eff δ_U + 11 K_U K_φ_eff − 9 c_Uphi^2 δ_U σ − 11 c_Uphi^2 σ) / [11 K_U (δ_U+1)]` | `kPhiEff − (cUphi^2 (11+9 deltaU) sigma)/(11(1+deltaU) kU)` |
| `M_supp` | `−88 (δ_U+1)(K_U c_etaphi + c_Uphi c_etaU)^2 / [π^2 K_U K_η_eff (epsη−1)(11 K_U K_φ_eff δ_U + … )]` | `(−88 (1+deltaU)(cEtaU cUphi + cB kU)^2) / ((−1+epsEta) kEtaEff kU π^2 (11(1+deltaU) kPhiEff kU − cUphi^2 (11+9 deltaU) sigma))` |
| `D_(phi z)` | `−δ_U g_U κ_0 κ_1 (g_B g_R − g_S g_W) / [K_U(δ_U+1)]` | `(deltaU gU (−(gB gR)+gS gW) kappa0 kappa1)/((1+deltaU) kU)` |
| `R_phi − R_U` | `K_U δ_U g_U (g_B g_R − g_S g_W) / [(δ_U+1)(K_U g_B + g_S g_U)(K_U g_W + g_R g_U)]` | `(deltaU gU (gB gR − gS gW) kU) / ((1+deltaU)(gS gU + gB kU)(gR gU + gW kU))` |

Each pair is the same expression under the snake_case ↔ camelCase mapping. `engines_agree: true`.

## Verdict justification

Three of the eight assertion pairs (`R_phi`, `D_phi`, `v.D_U.v`, `D_(phi z)`, `mismatch`) genuinely exercise the physics: they construct LHS by one route and RHS by another. The tracking-condition check at line 134/114 is also genuine — substituting `gS = gB gR/gW` and showing the difference collapses to zero is a real algebraic identity. But two of the central assertions (`A_phi^(eff) - expected` and `M_supp - expected`) are tautological by construction: each side is built from the same intermediate (`SU_expected` in the first; the quotient `numer/(denom)` with cancelling mu's in the second), so the residual is identically zero independent of the physics. These are filed as F2 and F3. The Mathematica script is a renamed copy of the SymPy script; that's F1.

The arithmetic itself holds up where it's exercised — I worked through `D_phi`, `v.D_U.v`, and the mismatch formula by hand and confirmed they reduce to the expected closed forms by genuine algebraic manipulation, not by construction. The output transcripts match my hand-computation. No sign errors found; no `kappa0^2 = (9/11)sigma vs (2/11)sigma` mix-up (sympy's `subs({kappa0**2: sigma - kappa1**2, kappa1**2: (2/11)sigma})` and Mathematica's pre-evaluated `kappa0^2 → sigma - (2/11)sigma` both correctly land on `(9/11)sigma`).

Verdict: `findings` (not stop-cold). The two tautological checks and the transliteration are correctable in-place without invalidating the genuine portion. No downstream propagation is at risk because the *closed-form* claims (R_phi, mismatch, etc.) are still verified by genuine checks; only the additional anchoring assertions need to be added.

## Self-test notes

- **Variable independence:** F3's proposed `sp.diff(Msupp_cont, mu_eta)` derivative is meaningful — `Msupp_cont` does explicitly contain `mu_eta` and `mu_phi` (line 108) in both numerator and denominator, so the derivative is non-trivially the test that they cancel. Not a fake zero. Same for `D[mSuppCont, muEta]` on the Mathematica side.
- **Parity / symmetry:** No integrals over unbounded domains in this unit; this trap does not apply.
- **Trivial-case pre-check:** For F2's proposed `Aphi_eff at delta_U = 0`, substituting `delta_U = 0` collapses `SU_expected = sigma/KU * (1 - 0) = sigma/KU` and `Aphi_eff_min = Kphi_eff - cUphi^2 sigma/KU`, which equals the proposed `Aphi_eff_min_expected`. Genuine zero. For the `split-vs-minimal overlap ratio`, plug in `delta_U = 1`: SymPy ratio = `(sigma/KU)(1 - (2/11)*(1/2)) / (sigma/KU) = 1 - 1/11 = 10/11`. Expected: `1 - (2/11)*(1/2) = 10/11`. Genuine match — and `10/11 ≠ 0`, so a sign error in `SU_expected` would be caught.
- **Path specifications:** No `missing_verification_script` findings; F1–F3 are edits to existing files. Paths verified at `scripts/...sympy_audit.py` and `mathematica/...mathematica_audit.wl`.
