---
unit_id: 024
batch: II.1
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 024 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.txt`

## What the script claims to verify

The scripts verify, on the canonical real STF l=2 basis {E20, E21c, E21s, E22c, E22s} normalized by NORM = sqrt(15/(8π)): (i) orthonormality of the basis under the quadratic overlap with the unit-sphere I4 isotropy measure (4π/15)(δδ+δδ+δδ); (ii) the corresponding angular source-map identity that gram·svec = svec; (iii) a grouped (20/21/22) bundle collapse on equal lanes and the two unequal-lane witness values for the (a_x, b_x) defect coordinates; (iv) the closed forms of B0/B2/B4 (two-mode BdG geometric series), Z0/Z2/Z4 and N0/N2/N4 (one-pair Maxwell/mixed coupled-pair series in ω with denominator Δ - Sω² + ω⁴), and D0/D2/D4 (sum law); (v) the exact axisymmetric Y20 triple-overlap matrix M = diag(κ*, κ*/2, κ*/2, -κ*, -κ*) with κ* = √5/(7√π) computed via the I6 = (4π/105)·sum-of-15-pairings measure; (vi) the grouped axisymmetric splitting weights (λ20, λ21, λ22) = (1, 1/2, -1) and the law b_x = 3·a_x; and (vii) the first-order ε-expansion of P_A = N_A/D_A and the matching law b_P = 3·a_P. The assertion that mhat_ang = 1 exactly in Section I.2 follows by construction from Section I.1's gram = I5 result.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 164 | `expect_zero("Gram - I5", gram - sp.eye(5))` | yes |
| A2 | sympy | 171 | `expect_zero("projected coefficients - source coefficients", projected - svec)` | partial (corollary of A1) |
| A3 | sympy | 193-195 | `xbar - x0`, `a_x on equal lanes`, `b_x on equal lanes` | yes (defines + checks) |
| A4 | sympy | 200-203 | unequal-lane witnesses for a_x, b_x | yes |
| A5 | sympy | 225-228 | `C_alpha - lambda_B,alpha I_etaalpha` | partial (definition consistency) |
| A6 | sympy | 237-239 | B0/B2/B4 sum formulas | yes |
| A7 | sympy | 275-286 | Z0/Z2/Z4, N0/N2/N4, D0/D2/D4 formulas | yes |
| A8 | sympy | 306 | `M - M_target` (κ* = √5/(7√π) diagonal) | yes |
| A9 | sympy | 323-326 | xbar, a_x, b_x, b_x - 3 a_x for grouped axisymmetric coords | yes |
| A10 | sympy | 359-371 | P20/P21/P22 expansions and grouped defects | yes |
| B1 | mathematica | 85 | `Gram - I5` | yes |
| B2 | mathematica | 90 | `projected coefficients - source coefficients` | partial (corollary) |
| B3 | mathematica | 98-104 | grouped-bundle collapse + witnesses | yes |
| B4 | mathematica | 114 | C_alpha definition consistency | partial |
| B5 | mathematica | 121-123 | B0/B2/B4 | yes |
| B6 | mathematica | 151-159 | Z0/Z2/Z4, N0/N2/N4, D0/D2/D4 | yes |
| B7 | mathematica | 167 | `m20 - mtarget` | yes |
| B8 | mathematica | 177-180 | xbar, a_x, b_x, b_x - 3 a_x | yes |
| B9 | mathematica | 191-200 | P20/P21/P22, grouped P-defects | yes |

The two "partial" rows (A2/B2 and A5/B4) are technically valid (their assertions can fail if their definitions are misspecified) but add little adversarial signal beyond the directly-anchored rows.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl:1-211`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage024_overlap_isotropy_sympy_audit.py:1-387`

**What's wrong:**
The `.wl` is a line-by-line port of the `.py`. The two engines do not derive the unit's results from independent formulations of the angular algebra; they execute the same algorithmic recipe in two languages. Concretely:

1. Identical recursive `pairings` helper. SymPy (lines 71-81) and Mathematica (lines 31-41) both define the same recursion (`a = first`, loop over remaining indices, prepend `(a, b)` pair, recurse on the rest). The Mathematica version's `Take[list, {2, i - 1}]` followed by `Drop[list, i]` mirrors the Python `lst[1:i] + lst[i+1:]` slice-merge exactly.

2. Identical I4 prefactor (4π/15) and three-pair sum, identical I6 prefactor (4π/105) and 15-pair sum. The hand-decomposed three pairings in `I4` are written in the same `{(0,1),(2,3)}, {(0,2),(1,3)}, {(0,3),(1,2)}` order (sympy line 105; mathematica lines 43-47).

3. Identical Section III intermediate symbols. SymPy lines 249-253 define `Delta`, `S`, `Q`, `H`, `P` with the formulas `OmegaU²·OmegaW² - Rr²`, `OmegaU² + OmegaW²`, `GU²·OmegaW² + 2·GU·GW·Rr + GW²·OmegaU²`, `GU² + GW²`, `OmegaU²·GW + Rr·GU`. Mathematica lines 132-136 define `deltaPair`, `sPair`, `qPair`, `hPair`, `pPair` with the same formulas in the same order. These are not standard physics names; they are the SymPy author's working shorthands. Their reuse verbatim is the strongest single tell.

4. Identical target formulas in identical order: `Z0 - Q/Δ`, `Z2 - (Q·S - H·Δ)/Δ²`, `Z4 - (Q·(S² - Δ) - S·H·Δ)/Δ³` (sympy 275-277 vs mathematica 151-153); `N0 - P²/Δ²`, `N2 - 2·P·(P·S - Δ·GW)/Δ³`, `N4 - (Δ²·GW² - 2·Δ·P² - 4·Δ·P·S·GW + 3·P²·S²)/Δ⁴` (sympy 278-283 vs mathematica 154-156).

5. Identical Section IV.2 reparameterization. SymPy lines 316-318 write `x20 = x0 + eps*lam20*x1` (with lam20 = 1, lam21 = 1/2, lam22 = -1), and the .wl writes `x20ax = x0 + eps*x1`, `x21ax = x0 + eps*x1/2`, `x22ax = x0 - eps*x1` (lines 171-173) — the same construction with the hardcoded λ values pre-substituted.

6. Identical Section V `lane_ratio` formulation. SymPy line 348 returns `expand(series((N0 + eps·lam·N1)/(D0 + eps·lam·D1), eps, 0, 2).removeO())`; Mathematica line 185 returns `Expand[Normal[Series[(n0 + eps·lam·n1)/(d0 + eps·lam·d1), {eps, 0, 1}]]]` — same numerator, same denominator, same expansion order, same target `(P0 + eps·P1)`, `(P0 + eps·P1/2)`, `(P0 - eps·P1)`.

7. Section ordering, banner strings, and witness substitutions match (e.g., the two witnesses `{x20 -> x0 + 1, x21 -> x0, x22 -> x0}` and `{x20 -> x0, x21 -> x0 + 1, x22 -> x0}` are the same two substitutions in the same order).

**Why this matters:**
The second-engine policy exists so that a coding bug or convention error in one engine's setup is caught by an independently formulated cross-check. When the `.wl` repeats the `.py`'s symbol naming, intermediate shorthand variables, and target expressions verbatim, both engines fail or pass together regardless of whether the underlying claim is correct. The current state offers only a "Mathematica's `FullSimplify` agrees with SymPy's `simplify` on identical input" guarantee, not an independent verification of the physics.

**Required change:**
Rewrite the Mathematica audit so the angular-integration layer is set up by a structurally different route, and the coupled-pair series uses Mathematica-native polynomial machinery rather than a port of the SymPy shorthand. Edit `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage024_overlap_isotropy_mathematica_audit.wl` in place (preserving the existing assertion list and PASS/FAIL contract), with the following structural changes:

1. Build the I4 and I6 isotropy moments by direct symbolic integration over the 2-sphere using Mathematica's `Integrate[..., {theta, 0, Pi}, {phi, 0, 2 Pi}]` against `Sin[theta]` weight and Cartesian-component products `n[i]` defined by `n[1] = Sin[theta] Cos[phi]; n[2] = Sin[theta] Sin[phi]; n[3] = Cos[theta]`. Replace the hand-written `pairings`/`KroneckerDelta` Wick sum with these explicit integrals. Wrap as `i4Direct[i_, j_, k_, l_]` and `i6Direct[i_, j_, k_, l_, m_, n_]`. The two definitions must coincide on all 81 (resp. 729) index choices, but should be encoded independently of the SymPy `pairings` recursion.

2. Replace the working symbols `deltaPair`, `sPair`, `qPair`, `hPair`, `pPair` with direct in-place expressions in each `Z_n`/`N_n` target. Do not introduce per-engine shorthand symbols whose names mirror the SymPy script; expand the closed forms in terms of `omegaU, omegaW, rPair, gU, gW` only and let Mathematica's `Together`/`Cancel` produce the rational form. The assertion targets remain unchanged in mathematical content but must not name `Q, H, P, S, Delta`.

3. Replace the recursive `pairings` helper (lines 31-41) with a Mathematica-native `Subsets`/`Permutations`-based enumerator (e.g., `pairings6Native[]` that enumerates the 15 perfect matchings of 6 elements via `Subsets[Range[6], {2}]` and recursion over remaining elements), or remove it entirely once the direct integration approach in change (1) is in place. In either case, the line-by-line correspondence to the SymPy recursion must be broken.

4. In Section IV.2, replace the pre-substituted forms `x20ax = x0 + eps*x1`, `x21ax = x0 + eps*x1/2`, `x22ax = x0 - eps*x1` with a generic `x[lam_] := x0 + eps*lam*x1` and call `x[1]`, `x[1/2]`, `x[-1]`. This matches the structural intent of the section (lambda-parameterized lanes) without echoing the SymPy hardcoded substitution.

5. Keep the existing assertion list (the names printed, the residuals checked, the PASS line ordering) intact. The verifier reads `expectZero` lines; do not alter the names or the mathematical content of the residuals.

**Verification:**
After Codex applies, the verifier confirms:
- `redteam exec-mathematica 024` exits 0 with the same PASS lines as the current output (all 28 PASS rows).
- Inspection of the rewritten `.wl` shows: (a) `Integrate[..., {theta, 0, Pi}, {phi, 0, 2 Pi}]` calls present in the angular-moment definitions, (b) no `Delta`, `Q`, `H`, `P`, `S` (or `deltaPair`, `qPair`, `hPair`, `pPair`, `sPair`) bindings, (c) no recursive `pairings` helper structurally parallel to the SymPy `pairings(lst)` (or that the helper has been removed in favor of the direct-integral route), (d) a single `x[lam_]` parameterizer instead of three pre-substituted axisymmetric lanes.

## Independent-derivation check (Mathematica)

The Mathematica script does not derive the angular overlap algebra independently from first principles. The clearest single piece of evidence is the verbatim reuse of the SymPy author's working symbols `Delta, S, Q, H, P` (renamed only by suffix `Pair`) for the Maxwell/mixed coupled-pair shorthand — these are not standard physics names but the SymPy script's local notation. A second tell is the byte-for-byte parallel `pairings` helper. A genuinely independent Mathematica derivation would build the angular moments via `Integrate` over the 2-sphere (Mathematica's native idiom), not by writing a Python-style index-pair recursion. This finding is filed as F1.

## Engine cross-check

Both engines produce all-zero residuals across the full assertion list. The .py output exits 0 with every `expect_zero` printing `= 0`; the .wl output exits 0 with every `expectZero` printing `= 0` followed by `PASS:`. Because the two scripts share the same intermediate algebra, the agreement is expected but not informative as an independent check (see F1).

## Verdict justification

The math holds up: the gram matrix is genuinely computed via I4 contraction and reduces to I5; the BdG and Maxwell/mixed series coefficients in ω match the rational-function targets at orders 0, 2, 4; the triple-overlap matrix M reduces to the claimed diag(κ*, κ*/2, κ*/2, -κ*, -κ*); the first-order P_A expansion follows from the standard ε-series of N_A/D_A. Attacks I tried that failed: (i) testing the I.2 "source-map identity" as a tautology — it is technically a corollary of I.1 but its `gram·svec - svec` form can fail if `gram ≠ I5`, so it is redundant rather than tautological; (ii) testing the κ* = √5/(7√π) constant as a hardcoded result — it is the claim being verified, and the I6-based LHS computes it independently, so the match is real; (iii) testing the IV.2 hardcoded λ values (1, 1/2, -1) as a `hardcoded_result` — they are the diagonal entries of the M matrix derived in IV.1 (modulo κ*) and act as parameter declarations rather than answer-against-answer comparisons; (iv) testing the Section V series for branch-cut issues — the expansion is in formal ε with D0 ≠ 0 explicitly assumed (nonzero=True), so no hidden division-by-zero. The single real attack that lands is the second-engine independence: F1 (`mathematica_transliteration`).

## Self-test notes

I checked: (1) no `sp.diff` or `D[expr, var]` calls appear in either script (no derivative-zero traps to verify), (2) the integrals in Section III/IV are over the unit 2-sphere with a measure that's manifestly even under n → -n only at even total power (I4, I6 are even-power moments) so the moments are nonzero by parity — consistent with the claimed nonzero diagonal entries; the Section IV.2 grouped reparameterization is a finite linear-algebra identity and contains no integrals to parity-check, (3) trivial-case pre-check of the F1 directive's "Required change" is structural (rebuild integrals, drop shorthand symbols, replace helper) and introduces no new `expectZero` rows to mentally evaluate — the existing assertion list is preserved. Path target in the directive is the `.wl` under `mathematica/`, which is the correct directory for the Mathematica engine.
