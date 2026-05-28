---
unit_id: 163
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage163_off_family_normal_coordinate.md]
  paper_appendix: present
---

# Audit unit 163 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_163.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage163_off_family_normal_coordinate.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex`
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.txt`

## What the paper claims

Stage 163 introduces the single first-order off-family normal coordinate
`delta_perp := delta g - g_-'(r_*) delta r`,
identifies it with the normalized parent compensation defect via `delta F = 4 sqrt(1+r_*^2) delta_perp` (notes Section 1; appendix Eq. `app-part04-deltaF-deltaperp`),
shows the compensation-ratio drift `delta R_q = -delta_perp/sqrt(1+r_*^2)` (notes Section 2),
gives the exact microscopic parent-variable closed form
`delta_perp = g_* dln(g_q K_s/(g_s lam)) + 1/(4 sqrt(1+r_*^2)) dln(K_s K_q / lam^2)` (notes Section 3; appendix Eq. `app-part04-deltaperp-microscopic`),
and routes outlet defects (delta C, delta E_2, delta E_4, Delta_Q) and the mouth-bias tangent/normal split through delta_perp (notes Sections 4-5). The paper card's purpose line states the audit target is the verification output for these identities; the checklist requires (a) deviations taken about the renormalized canonical point, (b) even-preservation gate, (c) tangent motion gives delta_perp = 0.

## What the script claims to verify

The SymPy and Mathematica scripts each set up symbols (r, g, dg, dr, microscopic dln_*, sigma_*, dkapW, dgamW, Sigma0, Sstar, dSigma0, dS), construct F = 1 + r^2 - 4(g-r)^2, R = (g-r)^2/(1+r^2), the lower-branch g_minus = r - sqrt(1+r^2)/2, and delta_perp = dg - g_-'(r) dr. They assert (i) delta F = 4 sqrt(1+r^2) delta_perp on the branch, (ii) delta R = -delta_perp/sqrt(1+r^2), (iii) tangent motion (dg = g_-'(r) dr) zeros delta F and delta R, (iv) the microscopic dln-expansion of delta_g - g_-' delta_r equals the closed-form delta_perp formula, (v) delta C = 4 sigma_* delta F / (1+r^2), (vi) the Stage 239 tangent/normal split of delta Pi factorizes correctly. They additionally print delta E_2, delta E_4, Delta_Q in their natural delta_C form and report Family-1 numerical coefficients.

## Paper <-> script cross-check

| Paper deliverable | Script check | Status |
|---|---|---|
| delta F = 4 sqrt(1+r_*^2) delta_perp | `expect_zero("delta F - 4 s delta_perp", dF - 4 s delta_perp)` | match |
| delta R_q = -delta_perp/sqrt(1+r_*^2) | `expect_zero("delta R + delta_perp/s", dR + delta_perp/s)` | match |
| Tangent motion gives delta_perp = 0 (checklist item c) | `dF.subs(dg, gp*dr) == 0`, `dR.subs(dg, gp*dr) == 0` | match (consequences of (i),(ii)) |
| Microscopic closed form for delta_perp | `expect_zero("microscopic delta_perp identity", delta_perp_micro - delta_perp_expected)` | match |
| delta C linkage to delta F | `expect_zero("delta C - 4 sigma_star deltaF/(1+r^2)", ...)` | match |
| Outlet defects delta E2, delta E4, Delta_Q linear in delta_perp | Defined and printed (not asserted as a zero residual) | partial — formula structure is shown but no explicit residual assertion against notes-form coefficient list 16/sqrt(1+r_*^2), 80/sqrt(1+r_*^2), 16/sqrt(1+r_*^2). The forms are equivalent by construction since delta C = 16 sigma_*/sqrt(1+r_*^2) delta_perp, so this is acceptable. |
| Tangent/normal split of delta Pi | `expect_zero("delta Pi tangent/normal split", dPi - dPi_expected)` | match |
| Family-1 numerical coefficients | Printed but not asserted | informational |

Paper alignment: aligned. Every load-bearing identity in the paper card and notes is covered by a non-trivial script assertion.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 60 | `expect_zero("delta F - 4 s delta_perp", dF - 4*s*delta_perp)` | delta F identity | yes |
| A2 | sympy | 61 | `expect_zero("delta R + delta_perp/s", dR + delta_perp/s)` | delta R_q identity | yes |
| A3 | sympy | 62 | tangent motion zeros delta F | checklist (c) | partial (corollary of A1) |
| A4 | sympy | 63 | tangent motion zeros delta R | checklist (c) | partial (corollary of A2) |
| A5 | sympy | 83 | microscopic delta_perp form | microscopic closed form | yes |
| A6 | sympy | 103-106 | delta C - 4 sigma_* dF/(1+r^2) | outlet linkage | yes |
| A7 | sympy | 117 | delta Pi tangent/normal split | mouth-bias decomposition | yes |
| B1 | mma | 41 | `expectZero["delta F - 4 s delta_perp", dF - 4*s*deltaPerp]` | delta F identity | yes |
| B2 | mma | 42 | `expectZero["delta R + delta_perp/s", dR + deltaPerp/s]` | delta R_q identity | yes |
| B3 | mma | 54 | microscopic delta_perp form | microscopic closed form | yes |
| B4 | mma | 69 | delta C - 4 sigma_* dF/(1+r^2) | outlet linkage | yes |
| B5 | mma | 78 | delta Pi tangent/normal split | mouth-bias decomposition | yes |

Both engines' assertion sets are substantively non-tautological. Each can fail under coefficient perturbations of g_-', delta r, delta g, microscopic-coefficient mismatches, or sign errors. The "tangent motion zeros" checks are corollaries of A1/A2 but still useful as substitution sanity probes.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage163_off_family_normal_coordinate_sympy_audit.py:42-119`
- `/var/projects/toy_throat_pde/mathematica/moving_throat_pde_stage163_off_family_normal_coordinate_mathematica_audit.wl:28-79`

**What's wrong:**

The Mathematica script is a near-line-for-line transliteration of the SymPy script. The variable choreography, computation steps, and assertion order are identical with only syntactic conversion. Quoted correspondences:

SymPy line 45-58:
```
r, g = sp.symbols("r g", positive=True, real=True)
s = sp.sqrt(1 + r**2)
gminus = r - s / 2
dg, dr = sp.symbols("dg dr", real=True)
gp = sp.diff(gminus, r)
F = 1 + r**2 - 4 * (g - r) ** 2
R = (g - r) ** 2 / (1 + r**2)
delta_perp = dg - gp * dr
dF = sp.diff(F, g).subs(g, gminus) * dg + sp.diff(F, r).subs(g, gminus) * dr
dR = sp.diff(R, g).subs(g, gminus) * dg + sp.diff(R, r).subs(g, gminus) * dr
```

Mathematica lines 28-39:
```
Clear[r, g, dg, dr];
$Assumptions = Element[{r, g, dg, dr}, Reals] && r > 0 && g > 0;
s = Sqrt[1 + r^2];
gMinus = r - s/2;
gPrime = D[gMinus, r];
fComp = 1 + r^2 - 4*(g - r)^2;
rComp = (g - r)^2/(1 + r^2);
deltaPerp = dg - gPrime*dr;
dF = (D[fComp, g] /. g -> gMinus)*dg + (D[fComp, r] /. g -> gMinus)*dr;
dR = (D[rComp, g] /. g -> gMinus)*dg + (D[rComp, r] /. g -> gMinus)*dr;
```

Same definition order, same construction method (`D[..., g] /. g -> gMinus`, mirroring `sp.diff(..., g).subs(g, gminus)`), same chosen residuals `dF - 4*s*deltaPerp` and `dR + deltaPerp/s`, same `deltaPerpExpected` form constructed in exactly the same way:

SymPy lines 79-82:
```
delta_perp_expected = sp.simplify(
    gminus * (dln_gq - dln_gs - dln_lam + dln_Ks)
    + (dln_Ks + dln_Kq - 2 * dln_lam) / (4 * s)
)
```

Mathematica lines 51-53:
```
deltaPerpExpected = FullSimplify[
  gMinus*(dlnGq - dlnGs - dlnLam + dlnKs) + (dlnKs + dlnKq - 2*dlnLam)/(4*s)
];
```

The mouth-bias section similarly mirrors `dPi = Expand[(1 - rStar*sStar)*dSigma0 - sigma0*(rStar*dS + sStar*dRFromPerp)]` (mma line 76) of `dPi = sp.expand((1 - Rstar*Sstar)*dSigma0 - Sigma0*(Rstar*dS + Sstar*dR_from_perp))` (sympy line 115).

The second-engine policy requires Mathematica to re-derive the identities independently (e.g., directly from F = 0 implicit-function definition of g_-, or from a power-series expansion in delta r about the lower branch, or by computing delta_perp via the gradient of F normalized by its magnitude). Instead Mathematica simply repeats the SymPy choreography in Mathematica syntax.

**Why this matters:**

A pure transliteration provides no cross-engine independence. If the SymPy script has a hidden algebraic error in how `dF` is constructed (e.g., wrong substitution order, wrong sign convention) the Mathematica script will reproduce the same error and "PASS" alongside. The whole point of a dual-engine audit is to catch exactly that.

**Required change:**

Refactor the Mathematica script so that at least two of its load-bearing identities are derived by an independent route from the SymPy script. Specifically:

1. Replace the `dF = (D[fComp,g] /. g -> gMinus)*dg + (D[fComp,r] /. g -> gMinus)*dr` total-derivative construction with the implicit-function route: compute `dgImp = -D[fComp, r]/D[fComp, g] /. g -> gMinus` and assert `Simplify[dgImp - gPrime] == 0`. Then verify `delta_perp` via `Series` expansion of `fComp[g, r] - fComp[gMinus[r0] + epsilon dg, r0 + epsilon dr]` to first order in epsilon, reading off the coefficient.

2. For the microscopic formula, replace the explicit `deltaG = gMinus*(dlnGq - dlnGs + 1/2*dlnKs - 1/2*dlnKq)` shortcut with a derivation that starts from `gExpr = gq*Sqrt[Ks]/(gs*Sqrt[Kq])`, computes `dgExpr = D[Log[gExpr], gq]*dlnGq*gq + ...` symbolically via the chain rule (e.g., via `Total[D[Log[gExpr], #]*Symbol["dln" <> SymbolName[#]]*# & /@ {gq, gs, Ks, Kq}]`), and then derives `delta_g` via `gExpr * (chain-rule sum)`. The script should not directly write the linearized form by hand and merely match it.

The remaining checks (delta C, delta Pi tangent/normal split) may stay parallel since their algebraic structure is straightforward, but the two checks above carry the bulk of the stage's content and must be independently derived.

**Verification:**

After patch, the Mathematica script contains at least two assertion blocks whose construction is structurally different from the SymPy script's, where the difference is not merely renaming of intermediate variables. The new `expectZero` lines should mention the cross-derivation route (e.g., `implicit-function route for g_-'`, `series-expansion route for delta_perp`).

## Independent-derivation check (Mathematica)

The Mathematica script is a transliteration of the SymPy script. See finding F1 above. Same residual choices (`dF - 4 s deltaPerp`, `dR + deltaPerp/s`), same explicit construction of `dF` and `dR` via differentiation followed by substitution, same hand-written form of `delta_g` and `delta_r` to which the Mathematica script then `FullSimplify`s the same expected closed form. No use of `Series`, `Solve` on `fComp == 0`, implicit-function theorem, or any other route distinct from the SymPy derivation.

## Engine cross-check

Both engines exit 0. Residuals all reduce to 0 in both. Sample concordance:

- SymPy output line 13: `delta F - 4 s delta_perp = 0` <-> Mathematica output line 13: same.
- SymPy line 14: `delta R + delta_perp/s = 0` <-> Mathematica line 15: same.
- SymPy line 18: `microscopic delta_perp identity = 0` <-> Mathematica line 18: same.
- SymPy line 24: `delta C - 4 sigma_star deltaF/(1+r^2) = 0` <-> Mathematica line 25: same.
- SymPy line 25: `delta Pi tangent/normal split = 0` <-> Mathematica line 27: same.

Family-1 numerical coefficients agree to all printed digits:
- `4 sqrt(1+r_*^2) ~ 8.1596676522425...` in both.
- `16 / sqrt(1+r_*^2) ~ 7.843456710202...` in both.
- `1/(4 sqrt(1+r_*^2)) ~ 0.1225540110969...` in both.

Engines agree, but the agreement is weakened by the transliteration finding: when both engines run the same derivation, agreement is structurally guaranteed and does not constitute genuine cross-validation.

## Verdict justification

The script's assertions exercise every load-bearing identity from the paper card's checklist and the notes' Sections 1-5: the delta F linkage, the delta R_q linkage, the microscopic delta_perp closed form, the outlet linkage delta C = 4 sigma_* delta F/(1+r^2), and the mouth-bias tangent/normal split. Each assertion is non-tautological (would fail under coefficient or sign perturbation). The microscopic closed-form check is the strongest: it derives `delta_perp_micro` from the chain-rule expansion of `delta r` and `delta g`, then matches it to the notes-quoted closed form. Numerical Family-1 coefficients carry forward correctly and are printed (not asserted, which is fine — they are downstream readbacks).

The script's banner says "STAGE 146" rather than "STAGE 163"; this is a label hygiene issue but does not affect the math.

The one substantive finding is that the Mathematica script is a transliteration of the SymPy script (F1, medium severity). The agreement between engines is therefore not an independent cross-check. Verdict: findings.

Attacks tried that did not break the script:
- Re-derived the microscopic closed form by hand from the chain-rule expansion of `delta r` and `delta g` against `delta_perp = delta_g - g_-'(r) delta_r`. Used the identity `r - r^2/(2 sqrt(1+r^2)) = gminus + 1/(2 sqrt(1+r^2))` to reduce. Result matched the script's `delta_perp_expected` form.
- Checked sign of `delta R = -delta_perp/sqrt(1+r^2)`: linearizing `R = (g-r)^2/(1+r^2)` at the lower branch gives `delta R = 2(g_- - r) delta_g / (1+r^2) - [(2(g_- - r)(-1))(1+r^2) - (g_- - r)^2 (2r)]/(1+r^2)^2 * delta_r`. On the lower branch `g_- - r = -s/2`, so `delta R = -s delta_g/(1+r^2) + [s(1+r^2) - (s^2/4)(2r)]/(1+r^2)^2 * delta_r`. Reducing: first term = -delta_g/s; second term = delta_r/s - (s^2)(r/2)/(1+r^2)^2 = delta_r/s - r/(2s^2). So delta R = -(delta_g - delta_r)/s ... wait this matches `-(delta_g - g_-' delta_r)/s = -delta_perp/s` only with the right g_-' simplification. Verified: g_-'(r) = 1 - r/(2s), and `r/(2 s^2) = (r/(2s))/s = (1 - g_-'(r))/s`; the algebra collapses correctly. Sign is correct.
- Checked delta F sign: `delta F = -8(g-r) delta g + (2r + 8(g-r)) delta r` at g = g_-: g_- - r = -s/2, so `delta F = -8(-s/2) delta g + (2r + 8(-s/2)) delta r = 4 s delta g + (2r - 4s) delta r`. With `g_-'(r) = 1 - r/(2s)`, `(2r - 4s) = -4(s - r/2) = -4s(1 - r/(2s)) = -4s g_-'(r)`. So `delta F = 4 s delta g - 4 s g_-'(r) delta r = 4 s delta_perp`. Sign correct.
- Tried to find a missing branch: the script only verifies the lower branch `g_-(r) = r - s/2`. Notes Section 1 also references an implicit upper branch by `mathcal F = 0`. However, paper card and notes Section 1 explicitly select the lower branch as the canonical one ("On the lower branch, g_-(r) = r - sqrt(1+r^2)/2"), and the entire downstream chain depends on the lower branch. Not a missing-branch finding.
- Tried to find `paper_misalignment`: appendix Eq `app-part04-deltaC-deltag` writes `delta C = 16 sigma_*/sqrt(1+r_F1^2) delta g`, whereas notes Section 4 and the script write `delta C = 16 sigma_*/sqrt(1+r_*^2) delta_perp`. These are not in conflict: the appendix expression is the natural form when delta_r = 0 is imposed (which then makes delta_perp = delta g - 0 = delta g). The appendix is summarizing in a frame where the even gate has already collapsed delta_r contributions; the stage 163 script chooses delta_perp as the proper general scalar. Both align — no finding.

## Self-test notes

Checked: (a) microscopic delta_perp algebraic derivation by hand — matches; (b) sign of delta F and delta R via direct linearization — matches; (c) Family-1 numerical coefficient cross-check between engines — matches to printed precision; (d) potential missing-branch concern (lower vs upper compensation branch) — paper explicitly selects lower, no finding; (e) potential paper_misalignment between appendix delta-g form and notes/script delta_perp form — reconciled via even-gate substitution, no finding. The transliteration finding (F1) is the only script-side defect.
