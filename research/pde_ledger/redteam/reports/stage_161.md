---
unit_id: 161
batch: IV.6
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage161_dn_similarity_slippage.md]
  paper_appendix: present
---

# Audit unit 161 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_161.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage161_dn_similarity_slippage.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_161}` on line 1356, no narrative row)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.txt`

## What the paper claims

The card's quoted output line states:

> "Last defect is `Delta_Q = -sigma_* Xi_slip dPi_tan / (1 - sigma_*)`, with `Xi_slip = Xi_gamma - 2 Xi_L`."

The notes flesh this out across seven sections:
1. Exact bare slippage `B_W = ((1+r_c)/9)(eps_gamma - eps_kappa)`.
2. Linearized slippage `dB_W = ((1+r_c*)/9)(d eps_gamma - d eps_kappa)`.
3. Exact D/N-tube even defect `eps_kappa = 12 L_W^2/(pi^2 a^2 (1+r_c)) - 1` and the difference identity `d eps_gamma - d eps_kappa = d ln gamma_0 - 2 d ln(L_W/a)` (in particular the `dr_c` cancellation).
4. Susceptibility decomposition `Upsilon_Pi = ((1+r_c*)/9)(Xi_gamma - 2 Xi_L)` and the collapse `Delta_Q = -(sigma_*/(1-sigma_*)) Xi_slip dPi_tan`, `N_Q - 1 = +(sigma_*/(1-sigma_*)) Xi_slip dPi_tan`.
5. D/N similarity-preservation theorem: if `Xi_gamma = 2 Xi_L` then `Xi_slip = 0` and `Delta_Q = 0`, `N_Q - 1 = 0`.
6. Numeric `(1+r_F1^2)/9 ~= 0.462362334687869` using `r_F1 ~= 1.77799353547498`.

The "Checks" sub-block of the card adds three meta-checks (deviations about the renormalized canonical point, even-preservation imposed before reading the odd defect, tangent motion on the parent family gives `delta_perp = 0`).

## What the script claims to verify

Both `.py` and `.wl` scripts (transcript-confirmed PASS) check the same six items in the same order:
1. `B_W - (1+r_c)*(eps_gamma - eps_kappa)/9 == 0`.
2. Coefficient of `eps^1` in a series of `((1+r_c* + eps*drc)/9)*(eps*deps_g - eps*deps_k)` equals `(1+r_c*)*(deps_g - deps_k)/9`.
3. `eps_kappa` from `kappa_0 = 4 L_W^2 / (pi^2 a^2)` reduces to `12 L_W^2 / (pi^2 a^2 (1+r_c)) - 1`; first variation on the branch reduces to `2 dL_W/L_W - 2 da/a - dr_c/(1+r_c)`.
4. `d eps_gamma = 9 dgamma_0/(1+r_c) - dr_c/(1+r_c)` rewritten as `d ln gamma_0 - d ln(1+r_c)` (after substituting `dln_gamma0 -> 9 dgamma0/(1+r_c)`); difference identity using the branch constraint.
5. `Upsilon_Pi = ((1+r_c*)/9)(Xi_gamma - 2 Xi_L)` and the collapse of `Delta_Q`, `N_Q - 1` to their `sigma_*`-only forms.
6. Substitution `Xi_gamma -> 2 Xi_L` makes `Delta_Q` and `N_Q - 1` vanish; numeric `(1+r_F1^2)/9 ~= 0.4623623346878688`. Also reduces `Delta_Q` in `(dSigma_0, dS)` and `(dThat, dS)` using carry-forward coefficients `0.832409471081635`, `-1.16275838754222`, and `6.42981496203006`.

The card's "Checks" 1-3 (deviations about canonical point, even-preservation, `delta_perp = 0`) are framing prerequisites for the entire moving-throat sub-block; they are not separate identities Stage 161 must symbolically verify, but the script does enforce them implicitly by linearizing about `eps_kappa = eps_gamma = 0` (Check 1) and only working with `dB_W` and odd defects (Check 2). `delta_perp = 0` is a prior-stage result not exercised here.

## Paper -> script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| `B_W = ((1+r_c)/9)(eps_gamma - eps_kappa)` | sympy line 44 / wl line 35 | match |
| `dB_W = ((1+r_c*)/9)(d eps_g - d eps_k)` | sympy line 52 / wl line 43 | partial (constructed-by-Taylor — algebraic tautology; see F2) |
| `eps_kappa = 12 L_W^2/(pi^2 a^2 (1+r_c)) - 1` | sympy line 58 / wl line 50 | match (printed, then differentiated) |
| `d eps_kappa = 2 dL_W/L_W - 2 da/a - dr_c/(1+r_c)` | sympy line 68 / wl line 63 | match (uses branch constraint) |
| `d eps_gamma = d ln gamma_0 - d ln(1+r_c)` | sympy line 77 / wl line 67 | mismatch (tautological by construction; see F1) |
| `d eps_gamma - d eps_kappa = d ln gamma_0 - 2 d ln(L_W/a)` | sympy line 81 / wl line 71 | match (uses branch constraint) |
| `Upsilon_Pi = ((1+r_c*)/9)(Xi_gamma - 2 Xi_L)` | sympy line 95 / wl line 90 | match (definition) |
| `Delta_Q = -(sigma_*/(1-sigma_*)) Xi_slip dPi_tan` | sympy line 102 / wl line 97 | match |
| `N_Q - 1 = +(sigma_*/(1-sigma_*)) Xi_slip dPi_tan` | sympy line 106 / wl line 101 | match |
| `Xi_gamma = 2 Xi_L => Delta_Q = N_Q - 1 = 0` | sympy lines 121-128 / wl lines 113-114 | match |
| `(1+r_F1^2)/9 ~= 0.462362334687869` | sympy line 131 / wl line 117 | match (numeric) |
| Stage 159 transport `dPi_tan ~ 0.8324 dSigma_0 - 1.1628 dS` | sympy line 112 / wl line 106 | extra (paper notes mention but this is downstream Stage 159 carry-forward, not a Stage 161 deliverable — informational, not a finding) |

Dominant pattern: aligned, with one tautology block and one Taylor-extraction tautology that should be hardened against the actual physical derivation.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero(BW - (1+rc)*(eps_g - eps_k)/9)` | claim 1 (B_W decomposition) | yes |
| A2 | sympy | 52 | `expect_zero(dBW - (1+rc_star)*(deps_g - deps_k)/9)` | claim 2 (dB_W linearization) | partial — checks series-coeff of constructed polynomial, not a derived quantity |
| A3 | sympy | 68 | `expect_zero(depsk_diff)` after branch substitution | claim 3 (`d eps_kappa` identity) | yes |
| A4 | sympy | 77 | `expect_zero(depsg_branch - (dln_gamma0 - drc/(1+rc)).subs(dln_gamma0, 9 dgamma0/(1+rc)))` | claim 4 (`d eps_gamma` rewrite) | no — tautological by construction |
| A5 | sympy | 88 | `expect_zero(diff_identity)` after branch substitution | claim 5 (difference identity) | yes |
| A6 | sympy | 102 | `expect_zero(DeltaQ + sigma * (Xi_g - 2 Xi_L) dPi / (1 - sigma))` | claim 6 (`Delta_Q` collapse) | yes (tests `(1+rc*)/9 <-> 9/((1+rc*))` cancellation) |
| A7 | sympy | 106 | `expect_zero(NQm1 - sigma * (Xi_g - 2 Xi_L) dPi / (1 - sigma))` | claim 7 (`N_Q - 1` collapse) | yes |
| A8 | sympy | 121 | `expect_zero(DeltaQ.subs(Xi_gamma, 2 Xi_L))` | claim 8 (preservation) | partial — substitution into a manifestly-proportional expression |
| A9 | sympy | 125 | `expect_zero(NQm1.subs(Xi_gamma, 2 Xi_L))` | claim 8 | partial |
| B1 | wl | 35 | mirror of A1 | claim 1 | yes |
| B2 | wl | 43 | mirror of A2 | claim 2 | partial |
| B3 | wl | 63 | mirror of A3 | claim 3 | yes |
| B4 | wl | 67 | mirror of A4 | claim 4 | no — tautological |
| B5 | wl | 83 | mirror of A5 | claim 5 | yes |
| B6 | wl | 97 | mirror of A6 | claim 6 | yes |
| B7 | wl | 101 | mirror of A7 | claim 7 | yes |
| B8 | wl | 113 | mirror of A8 | claim 8 | partial |
| B9 | wl | 114 | mirror of A9 | claim 8 | partial |

## Findings

### F1 — tautological_check

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py:72-80`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:65-70`

**What's wrong:**
In the sympy script the variable `depsg_branch` is defined on line 72 as

```
depsg_branch = sp.simplify(9 * dgamma0 / (1 + rc) - drc / (1 + rc))
```

and the very next assertion (lines 77-80) is

```
expect_zero(
    "d eps_gamma rewritten as d ln gamma0 - d ln(1+r_c)",
    depsg_branch - (dln_gamma0 - drc / (1 + rc)).subs(dln_gamma0, 9 * dgamma0 / (1 + rc)),
)
```

After the `.subs(dln_gamma0, 9*dgamma0/(1+rc))` the RHS of the subtraction becomes `9*dgamma0/(1+rc) - drc/(1+rc)`, which is exactly `depsg_branch`. The assertion therefore reduces to `depsg_branch - depsg_branch == 0`, algebraically guaranteed regardless of physics. The paper's claim is `d eps_gamma = d ln gamma_0 - d ln(1+r_c)`, which uses the branch-point relation `gamma_{0,*} = (1+r_c)/9`. The script never asserts that relation; it just renames `9*dgamma0/(1+rc)` to `dln_gamma0` and substitutes it back, so the identity `d ln gamma_0 = 9*dgamma_0/(1+r_c)` is never tested.

The Mathematica file does the same thing more nakedly at lines 67-70:

```
expectZero[
  "d eps_gamma rewritten as d ln gamma0 - d ln(1+r_c)",
  depsGBranch - (9*dgamma0/(1 + rc) - drc/(1 + rc))
];
```

Here `depsGBranch` was set on line 65 to `FullSimplify[9*dgamma0/(1 + rc) - drc/(1 + rc)]`, so the residual is identically zero by construction (no `subs` even needed).

**Why this matters:**
The notes describe this rewrite (notes lines 168-189) as the place where the odd defect's hybridization piece (`drc/(1+rc)`) is identified, justifying the cancellation that yields the boxed difference identity `d eps_gamma - d eps_kappa = d ln gamma_0 - 2 d ln(L_W/a)`. The script claims to "verify the rewrite" but actually verifies nothing — a sign error in the definition of `depsg_branch` would still produce zero residual here.

**Required change:**
Replace the tautology by a non-trivial check anchored to the branch-point relation `gamma_{0,*} = (1+r_c)/9`. Concretely: introduce `dln_gamma0` as a symbol, eliminate `dgamma0` via `dgamma0 = gamma_{0,*} * dln_gamma0 = (1+r_c)*dln_gamma0/9`, and assert that under this substitution `depsg_branch` reduces to `dln_gamma0 - drc/(1+rc)`. That exercises the physics step (`d gamma_0 = gamma_{0,*} d ln gamma_0` on the branch) instead of just renaming a sub-expression. In the Mathematica script make the same change, ensuring it is reproduced from the physics rather than mirroring the sympy edit.

**Verification:**
After the patch, the assertion at sympy line 77 should be checking `depsg_branch.subs(dgamma0, (1+rc)*dln_gamma0/9) - (dln_gamma0 - drc/(1+rc))` (or equivalent). A deliberate sign flip in `depsg_branch`'s definition should now make the assertion FAIL; pre-patch, it would still PASS. The Mathematica analog should fail in the same way under the same sabotage.

### F2 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py:49-52`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:40-43`

**What's wrong:**
The "linearized slippage law" check expands a hand-crafted polynomial in `eps`:

```
BW_lin = sp.expand((((1 + rc_star + eps * drc) / 9) * (eps * deps_g - eps * deps_k)).series(eps, 0, 2).removeO())
dBW = sp.expand(BW_lin.coeff(eps, 1))
expect_zero("linearized slippage law", dBW - (1 + rc_star) * (deps_g - deps_k) / 9)
```

The polynomial `((1 + rc_star + eps*drc)/9)*(eps*deps_g - eps*deps_k)` expands to `(1+rc_star)*eps*(deps_g - deps_k)/9 + drc*eps^2*(deps_g - deps_k)/9`. Extracting the `eps^1` coefficient mechanically returns `(1+rc_star)*(deps_g - deps_k)/9`, which is exactly what the assertion compares against. The check verifies sympy's `series`/`coeff` machinery, not that the linearization of `B_W = (1+r_c)(eps_g - eps_k)/9` about `eps_g = eps_k = 0` actually gives that form. The Mathematica `Series`/`Coefficient` analog (lines 40-43) is identical.

**Why this matters:**
A sign error or factor-of-two slip in the linearization step would not be caught here, because the script does not derive `dB_W` from the exact `B_W` formula on line 42 — it constructs a stand-in polynomial whose first-order coefficient is fixed by hand. The notes (lines 92-100) derive the linearization by perturbing the bare expression around the basepoint; the script should mirror that derivation.

**Required change:**
Replace the constructed polynomial with the actual linearization of the exact `BW` from line 42. Introduce small symbols `deps_k`, `deps_g`, define `BW_pert = BW.subs({eps_kappa: eps*deps_k, eps_gamma: eps*deps_g})`, compute `dBW = sp.diff(BW_pert, eps).subs(eps, 0)`, then assert `dBW - (1 + rc)*(deps_g - deps_k)/9 == 0` (or simplify with `rc -> rc_star` if a star-evaluated form is desired — but note the script then needs to be consistent about the basepoint). Apply the analogous derivation in Mathematica.

**Verification:**
A sabotage that changes `BW` on line 42 to `gamma0 - kappa0/2` (wrong factor) would, post-patch, make the line 52 assertion FAIL because the linearization residual would no longer match `(1+rc)*(deps_g - deps_k)/9`. Pre-patch it would still PASS since `dBW` is computed independently of `BW`.

### F3 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:1-131`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py:1-141`

**What's wrong:**
The Mathematica file is a near-literal port of the SymPy file. Corresponding evidence:

(a) Both open with an identical mislabelled banner `"STAGE 144 — D/N SIMILARITY SLIPPAGE DECOMPOSITION"` (sympy line 37, wl line 26). The wrong stage number "144" is preserved verbatim — an independent derivation would not reproduce a typo.

(b) Both scripts perform the same artificial `eps`-series construction for the linearization check (see F2). In Mathematica:
```
bWLin = Normal[Series[((1 + rcStar + eps*drc)/9)*(eps*depsG - eps*depsK), {eps, 0, 1}]];
dBW = Expand[Coefficient[bWLin, eps, 1]];
```
in SymPy:
```
BW_lin = sp.expand((((1 + rc_star + eps * drc) / 9) * (eps * deps_g - eps * deps_k)).series(eps, 0, 2).removeO())
dBW = sp.expand(BW_lin.coeff(eps, 1))
```
Same artificial polynomial, same variable choreography, same `eps^1` extraction. Mathematica would naturally use `D[expr, eps] /. eps -> 0` or `Series` over a physically meaningful quantity.

(c) Both scripts implement the exact same `depsk_branch` derivation using `subs(12*LW**2, pi**2 a**2 (1+rc))` / `/. 12*lW^2 -> Pi^2*a^2*(1 + rc)`; the multiplied-and-reduced residual `depsk_diff` / `depsKDiff` is recovered via the SAME multiplier `pi^2 * a^2 * (1+rc) * LW` (sympy line 66) and `PolynomialRemainder[..., -12*lW^2 + a^2*Pi^2*(1 + rc), lW]` (wl lines 57-61). An independent Mathematica derivation would more naturally reduce on the branch using `Eliminate`/`GroebnerBasis` over the branch ideal, or substitute `lW -> Pi*a*Sqrt[(1+rc)/12]` and `dLW -> ...`.

(d) The carry-forward print block (sympy lines 134-140, wl lines 120-126) is reproduced line-for-line, comment-for-comment.

This violates the second-engine policy: the Mathematica script does not independently re-derive any of the identities; it mirrors the SymPy choreography in Mathematica syntax.

**Why this matters:**
The point of having two engines is to catch transcription or simplification errors that any single CAS might make. A line-for-line port reproduces whatever defect (or whatever taut construction, cf. F1 and F2) the first engine had. The mislabelled banner is a tell that the wl file was generated/copied from the py file without re-derivation.

**Required change:**
Re-derive the Mathematica side from the underlying expressions rather than from the SymPy script structure. Concretely:
- For the bare slippage check (wl line 35): keep, this is the only structurally inevitable identity.
- For the linearization (wl lines 40-43): replace the `Series` of the artificial polynomial with `D[(gamma0 - kappa0/3) /. {epsKappa -> t*depsK, epsGamma -> t*depsG}, t] /. t -> 0`, then `FullSimplify` against `(1+rcStar)*(depsG - depsK)/9`. (See F2 — both engines change the same way.)
- For the `d eps_kappa` branch reduction (wl lines 52-62): replace the polynomial-remainder hack with a direct algebraic reduction, e.g. substitute `lW -> Pi*a*Sqrt[(1 + rc)/12]/Sqrt[1]` and `dLW -> ...` from the differentiated branch condition, then verify against `2*dLW/lW - 2*da/a - drc/(1 + rc)` using `FullSimplify`.
- For the difference identity (wl lines 71-83): same — derive `depsG - depsK` directly and reduce, do not mirror the SymPy multiplier-and-substitute trick.
- Correct the banner "STAGE 144" to "STAGE 161" at wl line 26 (and also at sympy line 37) — the wrong stage number indicates the file was duplicated without review.

**Verification:**
After patch, a side-by-side `diff` between the .py and .wl algebraic structure should show distinct simplification strategies (Series vs. D, subs vs. /. Eliminate, multiplier-trick vs. direct substitution). The banners should both read "STAGE 161 — D/N SIMILARITY SLIPPAGE DECOMPOSITION". Both scripts should still PASS with exit 0.

## Independent-derivation check (Mathematica)

The Mathematica `.wl` is a transliteration of the SymPy `.py`. Beyond the three evidence items quoted in F3, the assertion order, variable names, banner text (including the "STAGE 144" typo), and even the comment-only carry-forward block at the end match the SymPy script line for line. This is not a re-derivation; it is a translation. See F3.

## Engine cross-check

Both transcripts report exit 0 and "PASS". All printed expressions agree symbolically:

| Quantity | SymPy | Mathematica |
|---|---|---|
| `B_W` | `(eps_gamma - eps_kappa)*(r_c + 1)/9` | `((epsGamma - epsKappa)*(1 + rc))/9` |
| `dB_W` | `(deps_g + deps_g r_c_star - deps_k - deps_k r_c_star)/9` (rearranged) | `depsG/9 - depsK/9 + (depsG*rcStar)/9 - (depsK*rcStar)/9` |
| `eps_kappa` | `(12 LW^2 - pi^2 a^2 (r_c + 1))/(pi^2 a^2 (r_c + 1))` | `-1 + (12 lW^2)/(a^2 Pi^2 (1 + rc))` |
| `d eps_kappa` (on branch) | `(24 LW dLW - pi^2 a^2 dr_c - 2 pi^2 a da (r_c + 1))/(pi^2 a^2 (r_c + 1))` | `(-12 lW (a drc lW - 2 a dLW (1 + rc) + 2 da lW (1 + rc)))/(a^3 Pi^2 (1 + rc)^2)` — same up to multiplier `1/(12 LW^2) * (12 LW^2)` over the branch |
| `Upsilon_Pi` | `-(2 Xi_L - Xi_gamma)(r_c_star + 1)/9` | `((1 + rcStar)*(xiGamma - 2 xiL))/9` |
| `Delta_Q` | `dPi_tan sigma_* (-2 Xi_L + Xi_gamma)/(sigma_* - 1)` | `(dPiTan sigmaStar (xiGamma - 2 xiL))/(-1 + sigmaStar)` |
| `(1+r_F1^2)/9` | `0.462362334687868748` | `0.46236233468786880105466350900749...` |

No engine disagreement.

## Verdict justification

The script verifies every paper-side deliverable to the level claimed (B_W decomposition, linearization, even defect's D/N-tube form, the dr_c-cancellation difference identity, Upsilon_Pi, Delta_Q / N_Q-1 collapse, similarity-preservation, numeric prefactor). Paper alignment is `aligned`. However two of the load-bearing checks are tautological by construction (F1, F2) — they verify name choices rather than physics — and the Mathematica file is a transliteration of the SymPy file (F3), reproducing the same tautologies and even a "STAGE 144" mislabel. Verdict: `findings`, not `clean`. Not `UNFIXABLE` (each finding has a mechanical fix) and not `CRITICAL_DOWNSTREAM` (the underlying identities are still correct; the script simply needs to verify them honestly — no downstream stage's value would change).

## Self-test notes

Mentally walked the F1/F2 fixes: (1) variable independence — in the F2 fix `BW` depends on `eps_kappa, eps_gamma, r_c`; substituting `eps_kappa -> eps*deps_k` makes it a polynomial in `eps`, so `D[BW_pert, eps] /. eps -> 0` is nonzero (linear in `deps_g - deps_k` and `1+r_c`) — the residual against `(1+r_c)*(deps_g - deps_k)/9` simplifies to 0 honestly. (2) trivial-case pre-check — for F1 fix, substituting `dln_gamma0 = 0` and `drc = 0` makes `depsg_branch` reduce to zero only if `dgamma0 = 0`, consistent. (3) paper round-trip — proposed fixes use the same constants `(1+r_c)/9`, `12 L_W^2 / (pi^2 a^2 (1+r_c))`, etc. as the notes; no new constants introduced. F3's "banner mislabel" is a separate cosmetic fix that does not affect the math.
