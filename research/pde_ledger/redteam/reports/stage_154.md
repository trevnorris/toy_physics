---
unit_id: 154
batch: IV.6
auditor_model: claude-opus-4-7-1m
audit_date: 2026-05-27T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
paper_alignment: partial
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files: [moving_throat_pde_stage154_coevolving_core_mouth_map.md]
  paper_appendix: present
---

# Audit unit 154 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_154.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage154_coevolving_core_mouth_map.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (row at line 32; `\input{stages/stage_154}` at line 1342)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.txt`

## What the paper claims

The stage 154 card (line 16) quotes the bottom-line `Output`: "Nonlinear fixed-point map $\Sigma\propto e^{-\Phi_{\Sigma_0}[\Sigma]}$, with $\mathcal R[\Sigma]=(\mathfrak g[\Sigma]-\mathfrak r_{F1})^2/(1+\mathfrak r_{F1}^2)$." The notes file extends this into four deliverables: (1) the exact Family-1 core ratio law `R[Σ] = (g[Σ]-r_{F1})^2/(1+r_{F1}^2)` together with the integral kernels `T_s`, `T_q` defining `Phi_{Σ_0}`; (2) compensation equivalence `g = g_* ⟺ R = 1/4` with `g_* = r_{F1} - sqrt(1+r_{F1}^2)/2`; (3) first-order defect transport `δR = -δg/sqrt(1+r_{F1}^2) + δg^2/(1+r_{F1}^2)`; (4) local slope/bias identity `Π[Σ] = Σ_0[1 - R[Σ] S[Σ]]` and (with Stage 242 susceptibility closure) `Σ_0 = (20/9) T̂_m^2`. The card's `Checks` list also names two part-IV-level checks (even-preservation; tangent motion → `δ_⊥ = 0`) that belong to later units in the 154–163 block per the appendix row at part04 line 32.

## What the script claims to verify

The SymPy script (docstring at lines 5–11) lists four checks: (1) exact Family-1 ratio law `R(g) = (g-r)^2/(1+r^2)`; (2) compensation equivalence `R(g_*) = 1/4` with `g_* = r - sqrt(1+r^2)/2`; (3) the exact shifted-R expansion around `g_*`; (4) a linearized slope identity `dΠ = (1-R_* S_*) dΣ_0 - Σ_0 (R_* dS + S_* dR)`. The Mathematica script asserts the same three `expectZero` identities with the same expected forms. Banner labels in both scripts mis-print "STAGE 137" instead of stage 154 (sympy line 29, mathematica line 26). Neither script defines or exercises the integral kernels `T_s[Σ]`, `T_q[Σ]`, nor the fixed-point operator `Σ = e^{-Φ}/∫e^{-Φ}`.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| (notes §1) `R[Σ] = (g-r_{F1})^2/(1+r_{F1}^2)` | `R = (g - r)**2 / (1 + r**2)` (sympy line 32; wl line 31) used as a *definition*, plus the consistency check via `R(g_*) - 1/4` and the shifted-R expansion | `match` (definitional+derived consistency) |
| (notes §1, card line 16) Fixed-point ansatz `Σ ∝ e^{-Φ_{Σ_0}[Σ]}` with kernels `T_s`, `T_q` | none — kernels never defined, integral operator never instantiated | `missing` (likely deferred to later units in 154-163 block, but card lists this as the Output) |
| (notes §2) `g = g_* ⟺ R = 1/4`, with `g_* = r_{F1} - sqrt(1+r_{F1}^2)/2` | `expect_zero("R(g_star) - 1/4", ...)` (sympy line 36; wl line 35) | `match` |
| (notes §3) Exact `R(g_* + δg) = 1/4 - δg/sqrt(1+r^2) + δg^2/(1+r^2)` | `expect_zero("exact shifted R formula", R_shift - R_shift_expected)` (sympy line 40; wl line 39) | `match` |
| (notes §4) Linearized `dΠ` identity | `expect_zero("dPi identity", dPi - dPi_expected)` (sympy line 57; wl line 51) — uses placeholder symbols rather than evaluating `Φ'_{Σ_0}[Σ](0)`; the identity is a pure algebraic expansion of `(Σ_0+dΣ_0)(1-(R+dR)(S+dS))` | `partial` — confirms the polynomial linearization but does not anchor `Π` to `Φ'(0)` or to `1 - R·S` via differentiation of the kernel form |
| (notes §4, last paragraph) `Σ_0 = (20/9) T̂_m^2` susceptibility closure | none | `missing` (notes explicitly refer to Stage 242 carry-forward; reasonable to omit here) |
| (card Checks) even-preservation; `δ_⊥ = 0` on tangent motion | none | `missing` (part-IV level; appendix row at line 32 routes these to stages 155-163) |

Dominant pattern: the algebraic R-law and its first-order transport are faithfully exercised; the fixed-point operator and the kernel-anchoring of `Π` are not. `paper_alignment: partial`.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 36 | `expect_zero("R(g_star) - 1/4", R.subs(g, g_star) - 1/4)` | notes §2 compensation equivalence | yes |
| A2 | sympy | 40 | `expect_zero("exact shifted R formula", R_shift - R_shift_expected)` | notes §3 defect transport | yes |
| A3 | sympy | 57 | `expect_zero("dPi identity", dPi - dPi_expected)` | notes §4 slope identity (algebraic linearization only) | partial — checks `(Σ_0+dΣ_0)(1-(R+dR)(S+dS))` linearization, not `Π = Φ'_{Σ_0}(0)` via the kernel definition |
| A4 | mathematica | 35 | `expectZero["R(g_star) - 1/4", (rFun /. g -> gStar) - 1/4]` | notes §2 | yes |
| A5 | mathematica | 39 | `expectZero["exact shifted R formula", rShift - rShiftExpected]` | notes §3 | yes |
| A6 | mathematica | 51 | `expectZero["dPi identity", dPi - dPiExpected]` | notes §4 (algebraic linearization only) | partial — same scope as A3 |

A1/A4 are non-trivial: substituting `g_* = r - sqrt(1+r^2)/2` into `(g-r)^2/(1+r^2)` gives `((sqrt(1+r^2)/2)^2)/(1+r^2) = 1/4`. The identity could have failed (e.g., if the sign or normalization of `g_*` were wrong). Non-tautological. A2/A5 (shifted-R) are also non-trivial: a sign error in either the expected `-dg/sqrt(1+r^2)` term or the `dg^2/(1+r^2)` term would fail. A3/A6 sit closer to tautology — they verify that hand-dropping the cross-terms `dΣ_0·dR`, `dΣ_0·dS`, `dR·dS` from a polynomial product reproduces the textbook linearization formula. The "expected" expression mirrors the manual dropping rules, so the residual *must* be zero by algebra; the check has very low discriminating power.

## Findings

### F1 — mathematica_transliteration

**Severity:** medium
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl:26-51`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage154_coevolving_core_mouth_sympy_audit.py:29-57`

**What's wrong:**
The Mathematica script is a line-by-line transliteration of the SymPy script — same variable choreography, same intermediate steps, same hand-written `*_expected` forms, even the same incorrect banner text. The .wl file does not independently re-derive the identities (no use of `Series[]`, no algebraic manipulation distinct from the .py file's flow); it simply re-types the same substitutions in Mathematica syntax. This violates the two-engine policy: both engines must derive the result independently from the physical premises.

Quoted corresponding sections:

Sympy lines 31–40:
```
g, r, dg = sp.symbols("g r dg", real=True)
R = (g - r)**2 / (1 + r**2)
...
g_star = r - sp.sqrt(1 + r**2) / 2
expect_zero("R(g_star) - 1/4", R.subs(g, g_star) - sp.Rational(1, 4))
R_shift = sp.expand(R.subs(g, g_star + dg))
R_shift_expected = sp.expand(sp.Rational(1, 4) - dg / sp.sqrt(1 + r**2) + dg**2 / (1 + r**2))
expect_zero("exact shifted R formula", R_shift - R_shift_expected)
```

Mathematica lines 28–39:
```
Clear[g, r, dg];
$Assumptions = Element[{g, r, dg}, Reals];
rFun = (g - r)^2/(1 + r^2);
...
gStar = r - Sqrt[1 + r^2]/2;
expectZero["R(g_star) - 1/4", (rFun /. g -> gStar) - 1/4];
rShift = Expand[rFun /. g -> gStar + dg];
rShiftExpected = Expand[1/4 - dg/Sqrt[1 + r^2] + dg^2/(1 + r^2)];
expectZero["exact shifted R formula", rShift - rShiftExpected];
```

Identical: variable names (mapped `g_star ↔ gStar`, `R_shift ↔ rShift`, `R_shift_expected ↔ rShiftExpected`), identical algebraic order, identical hand-written expected expression. The same `1/4 - dg/sqrt(1+r^2) + dg^2/(1+r^2)` form appears in both scripts — neither obtains it via an independent expansion (e.g., `Series[rFun /. g -> gStar + dg, {dg, 0, 2}] // Normal` in Mathematica), they both copy the closed-form answer in by hand.

Sympy lines 45–57 vs Mathematica lines 46–51 follow the same pattern for `Pi`/`piExpr`:
- both define `(Σ_0+dΣ_0)(1-(R+dR)(S+dS))`,
- both drop the exact same four cross-terms by symbol substitution `{dSigma0*dR → 0, ...}` / `{dSigma0*dR -> 0, ...}`,
- both compare against the same hand-typed `dPi_expected = (1-R*S*)dΣ_0 - Σ_0(R*dS + S*dR)`.

A genuinely independent Mathematica derivation would (e.g.) `Series[piExpr, {dSigma0, 0, 1}, {dR, 0, 1}, {dS, 0, 1}] // Normal` and extract the linear part automatically, not replicate the SymPy substitution dictionary verbatim.

Additional smoking gun: both scripts mis-print `"STAGE 137 — EXACT CO-EVOLVING CORE-MOUTH MAP"` in the banner (sympy line 29, wl line 26). The shared typo confirms the .wl was authored by copying the .py rather than independently.

**Why this matters:**
A transliteration provides no genuine cross-engine check. Any algebraic mistake in the SymPy script (e.g., a wrong sign in `R_shift_expected`, a wrong coefficient on `dg^2`) would be faithfully copied into the .wl and both engines would report PASS. The whole point of the mathematica mirror is to catch precisely those errors via a second independent path.

**Required change:**
Rewrite the Mathematica audit so its derivations are independent. Concretely:
- Replace `rShift = Expand[rFun /. g -> gStar + dg]` and the hand-typed `rShiftExpected` with an independent `Series` expansion: e.g., `rShiftDerived = Normal[Series[rFun /. g -> gStar + dg, {dg, 0, 2}]];` and then `expectZero["shifted R consistency", rShiftDerived - (1/4 - dg/Sqrt[1+r^2] + dg^2/(1+r^2))]`. This makes the check `series(R, dg, 0, 2) == claimed_formula` rather than `Expand[manual_substitution] - manual_formula`.
- Replace the `dPi` block's manual cross-term dropping (`piLin = Expand[piExpr] /. {dSigma0*dR -> 0, dSigma0*dS -> 0, dR*dS -> 0, ...}`) with an independent multivariate linearization, e.g. `piLin = Normal[Series[piExpr, {dSigma0, 0, 1}, {dR, 0, 1}, {dS, 0, 1}]];` then keep the same `expectZero["dPi identity", dPi - dPiExpected]`. The expected form is the paper's claim; the linearization should now be computed by Mathematica's own series machinery, not by a substitution dictionary copied from the .py.
- Fix the banner string at line 26 from `"STAGE 137 — EXACT CO-EVOLVING CORE-MOUTH MAP"` to `"STAGE 154 — EXACT CO-EVOLVING CORE-MOUTH MAP"`.

**Verification:**
After the rewrite, the .wl should still print PASS for both expectZero checks, but the residuals must come from `Series[]`-driven computations rather than from substitution rules. The banner line at the top of the mathematica output transcript must read "STAGE 154". The verifier can confirm by reading the patched .wl and observing that no `dSigma0*dR -> 0` / `dR*dS -> 0` style hand-substitutions remain in the linearization block, and that `Series[...]` (or `Normal[Series[...]]`) appears.

## Independent-derivation check (Mathematica)

The Mathematica script is **not** an independent derivation. See F1 above for evidence. Three smoking guns:
1. Banner typo "STAGE 137" identical in both (.py line 29, .wl line 26).
2. Hand-typed `rShiftExpected = 1/4 - dg/Sqrt[1+r^2] + dg^2/(1+r^2)` in .wl line 38 is the same text the .py line 39 types. The .wl never independently expands `(gStar + dg - r)^2/(1+r^2)` to obtain it via `Series[]`.
3. The substitution dictionary `{dSigma0*dR -> 0, dSigma0*dS -> 0, dR*dS -> 0, dSigma0*dR*dS -> 0}` at .wl line 47 is character-identical (up to Mathematica syntax) to the .py dict at line 49–53.

## Engine cross-check

Both engines print residual `0` for all three `expectZero` checks, and both Exit 0. They agree — but that agreement is trivial given F1 (the .wl is a copy of the .py), so the agreement provides no independent cross-validation. Once F1 is fixed and the .wl uses Mathematica's own `Series[]`, the agreement will be meaningful.

## Verdict justification

Verdict is `findings` (not `clean`): the SymPy script's three algebraic checks (A1, A2, A3) are correct and non-tautological enough (especially A1 and A2) to anchor the notes §2 and §3 deliverables; the dPi check (A3) is closer to tautology but acceptable as a polynomial-identity sanity test. However, the Mathematica mirror is a transliteration, so the second-engine cross-check is effectively absent. F1 must be addressed. There is no UNFIXABLE issue, no paper_misalignment requiring user resolution, and no downstream propagation risk (later stages depend on the *correctness* of the algebraic identities, which the SymPy script does verify; the lack of an independent Mathematica derivation is a process defect, not a math defect). Outputs are fresh (sympy 12:47 > script 11:56; mathematica 13:16 > script 11:56). `paper_alignment: partial` because the card's headline `Output` (the fixed-point operator with `T_s`, `T_q` kernels) is never instantiated, only the algebraic skeleton `R[Σ]` of it; the appendix row at line 32 makes clear the integral-operator content is shared across stages 154–163, so this is partial rather than `missing`. Attacks tried that failed: (i) checking whether `g_*` is sign-flipped — direct substitution shows `g_* = r - sqrt(1+r^2)/2` yields exactly `R = 1/4`, sign correct; (ii) checking whether the shifted-R expansion drops a sign on the linear `-dg/sqrt(1+r^2)` term — the sign is consistent with `(g_*+dg-r)^2 = (dg - sqrt(1+r^2)/2)^2`, correct; (iii) checking whether the dPi expansion misses a `Σ_0 dR dS` cross-term — the .py explicitly zeros `dR*dS` (.py line 51) so the residual is by construction the linear part.

## Self-test notes

Traps checked:
1. Variable independence — N/A here, no `sp.diff`/`D[]` derivatives appear in the script. The "linearization" is done by polynomial expansion and term-dropping, not by differentiation, so the "derivative-of-something-it-doesn't-depend-on" trap is not in scope.
2. Parity/symmetry on integrals — N/A; no integrals are evaluated in this stage's scripts. The kernels `T_s`, `T_q` are defined in the notes but never instantiated in code.
3. Trivial-case pre-check — For F1's proposed Mathematica fix using `Series[piExpr, {dSigma0, 0, 1}, {dR, 0, 1}, {dS, 0, 1}]`, the linear part of `(Σ_0+dΣ_0)(1-(R_*+dR)(S_*+dS))` is `Σ_0(1-R_* S_*) + dΣ_0(1-R_* S_*) - Σ_0(R_* dS + S_* dR)`, which matches `dPi_expected + Pi0`. So the new `expectZero["dPi identity", dPi - dPiExpected]` will continue to evaluate to `0` after the fix; it remains a true PASS, not a silent-pass artifact.
4. Path specifications — F1 is not a missing_script finding; the .wl already exists and is edited in place at `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage154_coevolving_core_mouth_mathematica_audit.wl`. No new file paths are introduced.
5. Paper round-trip — F1's directive does not change which identities are verified, only how the Mathematica engine derives them; the paper card and notes remain consistent with the (unchanged) `expectZero` targets.
