---
unit_id: 013
batch: I.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-25T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 3
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: []
  paper_appendix: present
---

# Audit unit 013 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_013.tex`
- notes: (none)
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part01.tex` (row 48)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.txt`

## What the paper claims

Stage 013 is titled "Mouth-Taylor master map" and the `\stagefield{Output}` line states: "Stage~013 exports the mouth-Taylor primitive map \eqref{eq:stage007-z0n0}--\eqref{eq:stage007-z4}." That output bundles four explicit identities:

- `z_0 = ∂_s(Q/Delta)` and `n_0 = ∂_s(P^2/Delta^2)` (eq stage007-z0n0).
- `z_2 = ∂_s((Q*S_2 - H_port*Delta)/Delta^2)` (eq stage007-z2).
- `z_4 = ∂_s((Q*(S_2^2 - Delta) - H_port*S_2*Delta)/Delta^3)` (eq stage007-z4).

The "Gate feed" paragraph adds two narrative deliverables: (a) the first-order prefactor slope `Xi_load = n_0/N_0 + z_0/D_0` (eq stage013-xi), and (b) the statement that `z_0, z_2, z_4` feed the even gates `K_1, H_even` (whose explicit linear combinations are not given in this card — they appear in stage 014: `K_1 = -z_2 - z_0/9`, `H_even = -z_4 + (2/3) z_2 - z_0/27`). The card closes with "The audit checks the one-sided Taylor derivatives and their dependence on the primitive bottlenecks." The part-appendix row (line 48 of `stage_appendix_part01.tex`) summarises the deliverable as "Mouth-local Taylor map for projected source and flux corrections" with status `\StatusExactClosure{}`. No `notes/stages/moving_throat_pde_stage013_*.md` files exist.

## What the script claims to verify

The SymPy script (and its paired Mathematica file) audit five clusters of claims: (M1) the one-sided Taylor projection lemma — integrals of `exp(-u)` and `u*exp(-u)` over `[0, oo)` recover the leading two Taylor coefficients of `X(u) = X_0 + ell*u*X_1 + ...`, with mutated-coefficient assertions confirming the checks are non-trivial; (M2) the closed-form moment identity `∫₀^∞ u² e^{-u} du = 2`; (M3) the chain-rule expansion of `Z(ell) = (Q - H_port*ell^2)/(Delta - S_2*ell^2 + ell^4)` at orders ell^0, ell^2, ell^4 (the SymPy script asserts the literal polynomials; the Mathematica script independently derives them via `Series` of a primitive `Zsource`); (M4) the analogous chain-rule expansion of `N(ell) = (P - G_w*ell^2)^2/(Delta - S_2*ell^2 + ell^4)^2` yielding `n_0, n_2, n_4`; (M5) substantive coefficient probes — `∂Xi/∂P' = 2/P` and `∂(deltaP_2)/∂G_w' = -2P/(D_0*Delta^2)`, plus a nonzero `∂(deltaP_4)/∂G_w'`; (M6) a "mechanism sieve" — that the 2x2 system `K_1|_{Sx=Hx=0} = 0, H_even|_{Sx=Hx=0} = 0` has only the trivial solution `(Q', Delta') = 0` (and the dual for `S_2, H_port`).

## Paper ↔ script cross-check

| Paper deliverable | Script coverage | Status |
|---|---|---|
| `z_0 = ∂_s(Q/Delta)` (eq stage007-z0n0) | SymPy line 58 literal `(Delta*q1 - Q*d1)/Delta^2`; Mathematica M3 independent series derivation | match |
| `n_0 = ∂_s(P^2/Delta^2)` (eq stage007-z0n0) | SymPy line 80 literal `2*P*(Delta*p1 - P*d1)/Delta^3`; Mathematica M4 independent series derivation | match |
| `z_2 = ∂_s((Q*S_2 - H_port*Delta)/Delta^2)` (eq stage007-z2) | SymPy line 59; Mathematica M3 | match |
| `z_4 = ∂_s((Q(S_2^2 - Delta) - H_port*S_2*Delta)/Delta^3)` (eq stage007-z4) | SymPy lines 60-69; Mathematica M3 | match |
| `Xi_load = n_0/N_0 + z_0/D_0` (eq stage013-xi) | SymPy line 100 instantiates `n_0/N_0 → 2*p1/P - 2*d1/Delta` (silently assumes `N_0 = P^2/Delta^2`); only one partial derivative `∂Xi/∂P'` is tested, not the closed-form equality | partial |
| "z_0, z_2, z_4 feed even gates K_1, H_even" (narrative) | SymPy lines 101-102 use literal coefficients (-1/9, 2/3, -1/27); the only tests on K_1, H_even are sieve-determinant / linear-solve checks (lines 146-151) which are insensitive to the specific coefficient values | partial |
| "audit checks the one-sided Taylor derivatives" (narrative) | M1, M2 closed-form moments and series recovery; non-trivial | match |
| `deltaP_2`, `deltaP_4` formulas and their `G_w'` dependence (lines 106-124) | sympy/mathematica only — no paper-side text mentions `deltaP_2` or `deltaP_4` or the specific coefficient `-2P/(D_0*Delta^2)` | extra |
| mechanism sieve via 2x2 determinants on K_1/H_even (lines 126-151) | sympy/mathematica only — paper card does not mention a sieve | extra |

Setting `paper_alignment: partial` — the four equation-numbered deliverables match cleanly, but the narrative `Xi_load` deliverable is only weakly exercised, the K_1/H_even coefficients are essentially untested, and significant script content (`deltaP_2`, `deltaP_4`, sieve) has no paper counterpart.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 26 | `series(Xproj, ell, 0, 2).removeO() - (X_0 + ell*X_1) == 0` | one-sided Taylor lemma (narrative) | yes |
| A2 | sympy | 27-30 | mutant of A1 must be nonzero | one-sided Taylor lemma | yes |
| A3 | sympy | 34 | `∫u·u·e^{-u} du - 2 == 0` | one-sided Taylor lemma | yes |
| A4 | sympy | 35 | `series(Xproj_W2, ell, 0, 2).removeO() - (X_0 + ell*mu1*X_1) == 0` | one-sided Taylor lemma | yes |
| A5 | sympy | 36-39 | mutant of A4 must be nonzero | one-sided Taylor lemma | yes |
| A6 | sympy | 104 | `diff(Xi, Px) - 2/P == 0` | partial — Xi_load eq | partial (one term) |
| A7 | sympy | 122 | `diff(deltaP2_der, Gx) + 2*P/(D0*Delta^2) == 0` | none (script-only) | n/a |
| A8 | sympy | 123-124 | `diff(deltaP4_der, Gx) != 0` | none (script-only) | n/a |
| A9 | sympy | 144 | `diff(Xi, Px) != 0` | partial — Xi_load | partial |
| A10 | sympy | 145 | `diff(deltaP4_der, Gx) != 0` (duplicate of A8) | none | n/a |
| A11 | sympy | 146 | `qd_matrix.det() != 0` | partial — K_1, H_even narrative | weak |
| A12 | sympy | 147 | `sh_matrix.det() != 0` | partial — K_1, H_even narrative | weak |
| A13 | sympy | 148-149 | `qd_only == [{Qx:0, Dx:0}]` | partial — K_1, H_even narrative | weak |
| A14 | sympy | 150-151 | `sh_only == [{Sx:0, Hx:0}]` | partial — K_1, H_even narrative | weak |
| M1 | mathematica | 41-44 | `Normal[Series[Xproj, {ell,0,1}]] - (X_0 + ell X_1) == 0` | one-sided Taylor lemma | yes |
| M2a | mathematica | 46-49 | `∫u^2 e^{-u} du - 2 == 0` | one-sided Taylor lemma | yes |
| M2b | mathematica | 50-54 | series recovery for `u·e^{-u}` weight | one-sided Taylor lemma | yes |
| M3a | mathematica | 81 | `z0 - z0Expected == 0` (z0 derived independently via Series) | `z_0` paper deliverable | yes |
| M3b | mathematica | 82 | `z2 - z2Expected == 0` | `z_2` paper deliverable | yes |
| M3c | mathematica | 83 | `z4 - z4Expected == 0` | `z_4` paper deliverable | yes |
| M4a | mathematica | 104 | `n0 - n0Expected == 0` | `n_0` paper deliverable | yes |
| M4b | mathematica | 105 | `n2 - n2Expected == 0` | none (script-only n_2) | yes-to-script |
| M4c | mathematica | 106 | `n4 - n4Expected == 0` | none (script-only n_4) | yes-to-script |
| M5a | mathematica | 141 | `D[Xi, Px] - 2/P == 0` | partial — Xi_load | partial |
| M5b | mathematica | 142 | `D[deltaP2Der, Gx] + 2P/(D0 Delta^2) == 0` | none (script-only) | n/a |
| M5c | mathematica | 143 | `D[deltaP4Der, Gx] != 0` | none (script-only) | n/a |
| M6 | mathematica | 156-166 | sieve determinants nonzero; solves give only trivial pair | partial — K_1, H_even narrative | weak |

## Findings

### F1 — paper_misalignment

**Severity:** low
**Subtype:** `paper_missing_script_claim`
**Files:**
- paper: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_013.tex` (whole card)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:106-124, 126-151`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl:126-166`

**What's wrong:**
The script verifies content the paper card does not mention. Concretely:

1. `deltaP_2` and `deltaP_4` are defined in SymPy lines 106-119 and tested via the substantive identity `∂(deltaP_2)/∂G_w' = -2P/(D_0*Delta^2)` and a nonzero `∂(deltaP_4)/∂G_w'`. The paper card never introduces `deltaP_2` or `deltaP_4`, never quotes the coefficient `-2P/(D_0*Delta^2)`, and the part-appendix line only describes "source and flux corrections" without naming these symbols or their formulas.
2. The "mechanism sieve" (lines 126-151) — that the 2x2 systems on `(K_1, H_even)|_{Sx=Hx=0}` and `(K_1, H_even)|_{Qx=Dx=0}` have only trivial solutions — has no paper-side mention.
3. Both `n_2` and `n_4` are derived by the Mathematica script and asserted by SymPy, but the paper card only states `n_0 = ∂_s(P^2/Delta^2)`; `n_2, n_4` are absent.

**Why this matters:**
The script is doing more than the paper card sanctions. The paper-side appendix row "Mouth-local Taylor map for projected source and flux corrections" is broad enough to plausibly cover this work, but the stage card itself does not enumerate these deliverables, so a reader of the paper cannot reconstruct what the audit transcript verifies. Two resolutions are possible: (a) expand the paper card to declare `deltaP_2, deltaP_4, n_2, n_4` and the mechanism sieve as deliverables of stage 013, or (b) trim the script to only the four equation-numbered deliverables and move the sieve / `deltaP` machinery to whatever stage actually owns them (likely 014 given the K_1, H_even definitions live there). This is a paper-vs-script scope question, not a math error.

**Required change:**
See `## Resolve before fix_loop` block in the directive.

**Verification:**
Once the user picks (a) or (b), a follow-up directive instructs Codex to either edit the paper card (option a) or trim the scripts (option b).

### F2 — insufficient_verification

**Severity:** medium
**Files:**
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:100, 104, 144`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl:113-118, 141`

**What's wrong:**
The paper card declares `Xi_load = n_0/N_0 + z_0/D_0` (eq stage013-xi) as a primary deliverable. The script defines

  `Xi = simplify((2*p1/P - 2*d1/Delta + q1/(D0*Delta) - Q*d1/(D0*Delta**2)).subs(subs_der) / mu1)`

which silently instantiates `n_0/N_0 → 2*p1/P - 2*d1/Delta` and `z_0/D_0 → q1/(D0*Delta) - Q*d1/(D0*Delta^2)`. The substitution for `n_0/N_0` requires `N_0 = P^2/Delta^2`; this identification is never (a) anchored to an upstream stage, (b) stated as a comment, or (c) made explicit anywhere in the paper card. More importantly, the only test on Xi is `assert_zero("dXi/dPprime", sp.diff(Xi, Px) - 2/P)` (SymPy line 104; Mathematica M5a line 141). That assertion probes a single partial derivative of one summand and cannot detect (i) a wrong sign on `z_0/D_0`, (ii) a missing factor, (iii) an incorrect identification of `N_0`, or (iv) an absent `q1/(D0*Delta) - Q*d1/(D0*Delta^2)` block. The closed-form equality `Xi == n_0/N_0 + z_0/D_0` (with explicit `n_0, z_0` from the paper-side definitions) is never asserted.

**Why this matters:**
`Xi_load` is one of the named outputs in the stage card's "Gate feed" paragraph and feeds downstream into the prefactor slope. A "PASS" transcript here purchases too much: a reader would assume the paper's `Xi_load = n_0/N_0 + z_0/D_0` has been verified, but only a single derivative coefficient `∂Xi/∂P' = 2/P` is actually exercised. If the developer had written `Xi = 2*p1/P - 2*d1/Delta` (dropping the `z_0/D_0` half entirely), the existing assertion would still pass — the `z_0/D_0` part carries no `p1` and so contributes 0 to `∂Xi/∂P'`. The check is anchored to the wrong piece of Xi.

**Required change:**
In SymPy, add an explicit closed-form assertion. Concretely, before line 104, define
```
z0_form = (Delta*q1 - Q*d1)/Delta**2     # = paper's z_0
n0_form = 2*P*(Delta*p1 - P*d1)/Delta**3  # = paper's n_0
N0_form = P**2/Delta**2                   # explicit instantiation
Xi_paper = (n0_form/N0_form + z0_form/D0).subs(subs_der) / mu1
assert_zero("Xi matches paper closed form", sp.simplify(Xi - Xi_paper))
```
This (a) makes the `N_0 = P^2/Delta^2` identification explicit, (b) tests the full Xi structure not just one partial, and (c) adds zero hardcoded literals — both pieces use paper-side definitions. Mirror in the Mathematica script with a matching `assertZero["Xi matches paper closed form", FullSimplify[Xi - XiPaper]]`. Do **not** delete the existing `dXi/dPprime` check — keep it as the existing substantive coefficient probe.

**Verification:**
Verifier confirms the new `assert_zero("Xi matches paper closed form", ...)` line appears in both scripts, and that the saved outputs show a corresponding `residual = 0` print line for each.

### F3 — insufficient_verification

**Severity:** medium
**Files:**
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_sympy_audit.py:101-102, 146-151`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage013_projected_maxwell_mouth_taylor_master_mathematica_audit.wl:119-124, 156-166`

**What's wrong:**
The script defines
```
K1 = simplify((-(z2 + z0/9)).subs(subs_der) / mu1)         # line 101
H_even = simplify(((-z4 + Rational(2,3)*z2 - z0/27)).subs(subs_der) / mu1)  # line 102
```
with hardcoded coefficients `1/9`, `2/3`, `-1/27`. These coefficients are the load-bearing physics of the gate definitions (they encode the projected weights of the Taylor coefficients onto the even-mode gates), and they match the paper's stage 014 definitions exactly (`K_1 = -z_2 - z_0/9`, `H_even = -z_4 + (2/3) z_2 - z_0/27`). Stage 013 does not state these formulas in its own card; it only narratively says "z_0, z_2, z_4 feed the even gates K_1 and H_even."

The script's only tests on K_1 / H_even are the sieve checks at lines 146-151 (and Mathematica 156-166): `assert_nonzero` on `qd_matrix.det()` and `sh_matrix.det()`, plus `qd_only == [{Qx:0, Dx:0}]` and `sh_only == [{Sx:0, Hx:0}]`. By inspection, these tests are insensitive to the specific values of the coefficients `1/9, 2/3, -1/27`. For example, with `Sx=Hx=0`:

  `K_1|_{Sx=Hx=0} = 2*Q*S2*Dx/Delta^3 - Qx/(9*Delta) + Q*Dx/(9*Delta^2)`

If one replaced `1/9` with `1/8`, the resulting `qd_matrix` would still have a generically nonzero determinant and the linear solve `qd_matrix * (Qx, Dx)^T == 0` would still yield the unique trivial solution `(0, 0)`. The same applies to `2/3` and `-1/27` in `H_even`. Thus the script reports "PASS" for any nonzero coefficient triple, not specifically for `(1/9, 2/3, -1/27)`. The literal values that are the actual physics carry-forward to stage 014 are never actively tested.

**Why this matters:**
A reader sees that K_1, H_even are referenced and that their derived sieve has the expected trivial-solution structure, and infers that the gate definitions are correct. They aren't — they're correct only because the developer typed in the right numbers. A single-character typo in `K1 = -(z2 + z0/9)` → `K1 = -(z2 + z0/91)` would still PASS every assertion in the current script (the sieve determinant remains a generically nonzero rational, and the linear solve still returns only the trivial pair). The gate coefficients are the load-bearing input to stage 014, so a silent error here would propagate undetected.

**Required change:**
Add coefficient-specific assertions in both engines that fail if `1/9`, `2/3`, or `-1/27` is altered. The minimal-scope fix is to substitute concrete values for the source-side primitives and check that K_1, H_even reduce to the expected numbers. Specifically, in SymPy after line 102, add
```
# Coefficient sanity: with q1=1, d1=0, s1=0, h1=0 (only Q' nonzero),
# K_1 should reduce to -z_0/9 = -1/(9*Delta), H_even to -z_0/27 = -1/(27*Delta).
K1_only_q = K1.subs(subs_der).subs({Qx:1, Dx:0, Sx:0, Hx:0, Px:0, Gx:0}) * mu1
H_only_q = H_even.subs(subs_der).subs({Qx:1, Dx:0, Sx:0, Hx:0, Px:0, Gx:0}) * mu1
assert_zero("K1 = -z0/9 when only Q' nonzero", K1_only_q - (-1)/(9*Delta))
assert_zero("H_even = -z0/27 when only Q' nonzero", H_only_q - (-1)/(27*Delta))

# Coefficient sanity: with s1=1 only, K_1 should reduce to -∂z2/∂s1 weight = -Q/Delta^2
# and H_even should reduce to (2/3)·(-Q/Delta^2) + 0 = -2*Q/(3*Delta^2).
K1_only_s = K1.subs(subs_der).subs({Qx:0, Dx:0, Sx:1, Hx:0, Px:0, Gx:0}) * mu1
H_only_s = H_even.subs(subs_der).subs({Qx:0, Dx:0, Sx:1, Hx:0, Px:0, Gx:0}) * mu1
assert_zero("K1 = -Q/Delta^2 when only S2' nonzero", K1_only_s - (-Q/Delta**2))
assert_zero("H_even = -2*Q/(3*Delta^2) when only S2' nonzero", H_only_s - (-Q*sp.Rational(2,3)/Delta**2 - (-2*Q*S2*1/Delta**3 + 0) ))
```

Actually the last assertion is fragile because `-z_4 + (2/3) z_2 - z_0/27` with `s1=1, q1=d1=h1=0` gives:

  z_0 = 0; z_2 = (Delta·Q·1)/Delta^3 = Q/Delta^2; z_4 = (-Delta^2·Hport·1 + 2·Delta·Q·S2·1)/Delta^4 = -Hport/Delta^2 + 2·Q·S2/Delta^3.

So `H_even = -z_4 + (2/3) z_2 - z_0/27 = Hport/Delta^2 - 2·Q·S2/Delta^3 + (2/3)·Q/Delta^2 - 0`. Hence the assertion should be:
```
assert_zero(
    "H_even when only S2' nonzero",
    H_only_s - (Hport/Delta**2 - 2*Q*S2/Delta**3 + sp.Rational(2,3)*Q/Delta**2)
)
```
The Codex directive includes the trap-aware form. Mirror these in the Mathematica script with `assertZero["K1 = -1/(9 Delta) when only Qx nonzero", K1 /. ... - (-1/(9 Delta))]` etc.

**Verification:**
Verifier confirms the four new assertions appear (two for K_1, two for H_even, in each engine) and that perturbing `1/9 → 1/8` in K_1 would now make `K1 = -z0/9 when only Q' nonzero` fail. (The verifier doesn't perturb; this is a paper round-trip check the auditor did mentally.)

## Independent-derivation check (Mathematica)

The Mathematica script is **not** a transliteration of the SymPy script — for the load-bearing M3 and M4 claims (the `z_n, n_n` chain-rule coefficients), the Mathematica script independently derives the coefficients by Taylor expansion of a primitive:

  SymPy lines 58-69: `z0 = (Delta*q1 - Q*d1)/Delta**2` (literal polynomial).
  Mathematica lines 66-70: `Zsource[x_] := (qFun[t] - hFun[t] x^2)/(dFun[t] - sFun[t] x^2 + x^4); Zexpansion = Normal[Series[Zsource[ell], {ell, 0, 4}]]; z0 = timeDerivativeAtZero[Coefficient[Zexpansion, ell, 0]];`

The Mathematica script then compares its independently derived `z0, z2, z4, n0, n2, n4` against the literal `Expected` forms (which match the SymPy literals). This is the textbook second-engine cross-check pattern. For the M1, M2 lemmas the two engines do compute the same integrals through structurally similar Mathematica/Sympy code, but that is unavoidable for closed-form moments and is not a transliteration concern. For M5 and M6 both engines apply the same `subsDer` substitution to the same `(K1, H_even, Xi, deltaP2, deltaP4)` formulas — this is parallel use of the same algebraic definitions, not a derivation, so no transliteration finding here either; the substantive findings are F2 and F3.

## Engine cross-check

Both engines run to PASS. Final-line residuals match:

- SymPy: "STATUS: PASS" (no per-assertion residuals in the saved transcript, but the script's `assert_zero` would raise on any nonzero residual; the bare "PASS" implies all 14 assertions returned 0/nonzero as expected).
- Mathematica: 19 `OK ... residual = 0` lines for the `assertZero` checks; 3 `OK ... residual = <expr>` lines for the `assertNonzero` checks where the expr is the (nonzero) symbolic residual. The nonzero residuals (sieve determinants and the `deltaP4` Gprime dependence) are well-formed rational expressions in the expected symbols, no surprise terms.

No engine disagreement.

## Verdict justification

The four equation-numbered paper deliverables (z_0, z_2, z_4, n_0) are now well-anchored: SymPy asserts literal polynomials, Mathematica independently derives them via series expansion of a primitive, and the two engines agree. The one-sided Taylor lemma (M1, M2) is real algebra, non-tautological, and cross-engine verified. The Mathematica script's existence and structure resolve the prior F1 (`missing_mathematica`). The current SymPy file's tautological-independence assertions from the v1 report have been removed (line 86-90 of v1 are gone), resolving v1's F2. The hardcoded chain-rule literals now have explicit anchor comments pointing to the Mathematica independent derivation (lines 49-57, 71-79), substantively resolving v1's F3. What remains is a layer of script-side scope (deltaP_2, deltaP_4, sieve) that the paper card doesn't sanction (F1, paper_misalignment — user-resolution), a weak test on the named `Xi_load` deliverable (F2, insufficient_verification — Codex-fixable), and a complete absence of coefficient-specific tests on K_1 and H_even literals (F3, insufficient_verification — Codex-fixable). None are stop-cold. Verdict: `findings`, three findings, paper_alignment `partial`.

## Self-test notes

I checked the variable-independence trap (audit prompt step 1) on the remaining assertions: every `sp.diff(EXPR, VAR)` in lines 104, 122, 123-124, 144, 145 has VAR among EXPR's symbols, so no tautological-derivative false-positive. I checked parity (step 2): integrals are over `[0, ∞)`, not symmetric, parity trap inapplicable. I checked the trivial-case substitution for F3's proposed assertions (step 3): with `s1=1, q1=d1=h1=0`, `z_2 = Q/Delta^2` and `z_4 = -Hport/Delta^2 + 2·Q·S2/Delta^3`, giving `H_even = Hport/Delta^2 - 2·Q·S2/Delta^3 + (2/3)·Q/Delta^2`, which matches the directive's expected RHS — the proposed assertion would PASS for the correct script and FAIL if `2/3` were perturbed to (say) `2/5`. Path specifications (step 4): all targeted paths exist; no new file creation needed.
