---
unit_id: 109
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
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
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage109_linearized_branch_selection.md
  paper_appendix: present
---

# Audit unit 109 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_109.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage109_linearized_branch_selection.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (only `\input{stages/stage_109}` line at L1252 — no prose row for this unit)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.txt`

## What the paper claims

The card designates the boxed `\begin{quote}` as the audit target: "First-order outgoing defect is `Delta_Q = 5b + a_0/3 + 9 a_5`" (stage_109.tex L16). The notes elaborate the derivation: linearize `chi_Q = 3(S*beta^5 + 9 Sigma_5)/(3S - Sigma_0)` around `S=1+eps s`, `beta=1+eps b`, `Sigma_0=eps a_0`, `Sigma_5=eps a_5`, obtaining the boxed expansion `chi_Q = 1 + eps(5b + a_0/3 + 9 a_5) + O(eps^2)`; the preservation condition `5b + a_0/3 + 9 a_5 = 0`; and the equivalent form `a_5 = -5b/9 - a_0/27`. The notes additionally state the immediate implications: (i) the `s` coefficient drops out, (ii) the `b` shift contributes linearly. The card's `Checks` field lists three secondary checks: (a) "pure scale and pure argument deformations separately"; (b) "Robin and standalone mixed-pole limits before imposing compensation"; (c) "compensated branch preserves even coefficients as well as odd normalization." The notes do not mention items (b) or (c).

## What the script claims to verify

The sympy script builds `chi = 3*(S*beta^5 + 9*Sigma5)/(3*S - Sigma0)` with the linearized ansatz; takes a series in `eps` to order 1; subtracts `expected = 1 + eps*(5*b + a0/3 + 9*a5)` and asserts the residual is zero; takes `coeff = (chi_series - 1)/eps` and asserts `d coeff/ds == 0` (scale cancellation); calls `sp.solve(sp.Eq(coeff, 0), a5)` to obtain `a5_sol` and asserts `coeff.subs(a5, a5_sol) == 0`. The Mathematica script mirrors this structure but adds one extra non-tautological check: `expectZero["a5 preservation condition + 5 b/9 + a0/27", a5Pres + 5*b/9 + a0/27]`, which anchors the Solve output to the paper's closed-form `a_5 = -5b/9 - a_0/27`.

## Paper ↔ script cross-check

| Paper deliverable | Script-side check | Status |
|---|---|---|
| Boxed audit target: `Delta_Q = 5b + a_0/3 + 9 a_5` (card L16) | sympy L34 `expect_zero('linearized chi law', chi_series - expected)` + mathematica L41 | match |
| Notes box: `chi_Q = 1 + eps(5b + a_0/3 + 9 a_5) + O(eps^2)` | same as above | match |
| Notes box: preservation `5b + a_0/3 + 9 a_5 = 0` | implicit in sympy L41 `sp.solve(sp.Eq(coeff, 0), a5)` and printed | partial (printed but not directly asserted in sympy; mathematica L49 asserts it) |
| Notes box: `a_5 = -5b/9 - a_0/27` | sympy: tautological subs at L43; mathematica: L49 anchors to paper closed-form (non-tautological) | partial (mathematica anchored, sympy tautological) |
| Notes implication (i): `s` cancels | sympy L38 `expect_zero('overall scale cancels', sp.diff(coeff, s))` + mathematica L45 | match |
| Card Check (a): pure scale AND pure argument deformations checked separately | scale: yes (sp.diff w.r.t. `s`); argument: no separate `sp.diff(coeff, b)` to anchor the `5b` coefficient | partial |
| Card Check (b): Robin and mixed-pole limits before compensation | (no script-side check) | missing |
| Card Check (c): compensated branch preserves even coefficients as well as odd normalization | only odd normalization (`coeff` vanishing) checked; even-coefficient preservation built into ansatz and not separately asserted | partial |

`paper_alignment` = `partial` — bottom-line audit target matches exactly; secondary card "Checks" (b)/(c)/part of (a) have no script-side counterpart, though the notes do not echo these checks.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 34 | `expect_zero('linearized chi law', chi_series - expected)` | boxed audit target / notes box | yes |
| A2 | sympy | 38 | `expect_zero('overall scale cancels', sp.diff(coeff, s))` | notes implication (i); card Check (a) part 1 | yes |
| A3 | sympy | 43 | `expect_zero('preservation substitution', coeff.subs(a5, a5_sol))` | notes preservation box | no (tautological — `a5_sol` constructed to make `coeff` vanish) |
| A4 | mathematica | 41 | `expectZero["linearized chi law", chiSeries - expected]` | boxed audit target / notes box | yes |
| A5 | mathematica | 45 | `expectZero["overall scale cancels", D[coeff, s]]` | notes implication (i) | yes |
| A6 | mathematica | 49 | `expectZero["a5 preservation condition + 5 b/9 + a0/27", a5Pres + 5*b/9 + a0/27]` | notes preservation closed-form | yes (anchors Solve output to paper) |
| A7 | mathematica | 50 | `expectZero["preservation substitution", coeff /. a5 -> a5Pres]` | notes preservation box | no (tautological — same as A3) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:41-43`

**What's wrong:**
The sympy "preservation substitution" check is tautological by construction. At L41 the script defines `a5_sol = sp.solve(sp.Eq(coeff, 0), a5)[0]` — i.e., `a5_sol` is defined as the value of `a5` that makes `coeff == 0`. The subsequent L43 assertion `expect_zero('preservation substitution', coeff.subs(a5, a5_sol))` is then guaranteed to be 0 by the way `a5_sol` was constructed; it cannot fail regardless of whether `coeff` has the physically correct form.

The Mathematica script has the analogous tautology at L50 but mitigates it with an anchored equality check at L49: `expectZero["a5 preservation condition + 5 b/9 + a0/27", a5Pres + 5*b/9 + a0/27]`. This compares the Solve output to the paper's closed-form `a_5 = -5b/9 - a_0/27` (notes box L70-72) and is non-tautological. The sympy script has no equivalent anchored check.

**Why this matters:**
Without anchoring the Solve output to the paper's closed-form, the sympy script never verifies the second of the notes' two boxed preservation expressions (`a_5 = -5b/9 - a_0/27`). A regression that miscomputed `coeff` (wrong coefficient on `b`, `a0`, or `a5`) would yield a wrong `a5_sol` but still pass the tautological L43 check. The wrong sign or factor would only be caught by A1; the second-channel cross-check is absent.

**Required change:**
Add an anchored equality check between `a5_sol` and the paper's closed-form. After L42 (the print) and before L43, insert:

```python
expected_a5_sol = -sp.Rational(5, 9)*b - sp.Rational(1, 27)*a0
expect_zero('a5 preservation closed-form', sp.simplify(a5_sol - expected_a5_sol))
```

Keep L43 as-is (it documents the substitution closes the loop, even if tautologically). The new check is non-tautological because the right-hand side is hard-coded from the paper notes and would catch a miscomputed `coeff`.

**Verification:**
After Codex applies, the sympy script must contain a new `expect_zero` call comparing `a5_sol` to `-5*b/9 - a0/27`. `redteam exec-sympy 109` must exit 0 and the new line should appear in the output between "a5 preservation condition" and "preservation substitution".

### F2 — mathematica_transliteration

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:31-50`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:25-43`

**What's wrong:**
The Mathematica script's algebraic flow mirrors the sympy script step-by-step rather than independently re-deriving the result. Corresponding sections:

(1) Linearized ansatz definitions, sympy L25-28 vs. mathematica L31-34: same five symbol introductions (`S`/`sNorm`, `beta`, `Sigma0`/`sigma0`, `Sigma5`/`sigma5`) with identical structural form `1 + eps*X` and `eps*X`. Same letter set.

(2) Build–series–subtract, sympy L29-34:
```
chi = sp.simplify(3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0))
chi_series = sp.series(chi, eps, 0, 2).removeO()
expected = 1 + eps*(5*b + a0/sp.Integer(3) + 9*a5)
expect_zero('linearized chi law', chi_series - expected)
```
vs. mathematica L36-41:
```
chiQ = FullSimplify[3*(sNorm*beta^5 + 9*sigma5)/(3*sNorm - sigma0), ...]
chiSeries = Expand[Normal[Series[chiQ, {eps, 0, 1}]]]
expected = 1 + eps*(5*b + a0/3 + 9*a5)
expectZero["linearized chi law", chiSeries - expected]
```
Same pattern: simplify → series → subtract expected → assert zero.

(3) Coefficient extraction and Solve, sympy L36-41:
```
coeff = sp.expand((chi_series - 1)/eps)
expect_zero('overall scale cancels', sp.diff(coeff, s))
a5_sol = sp.solve(sp.Eq(coeff, 0), a5)[0]
```
vs. mathematica L43-47:
```
coeff = FullSimplify[Expand[(chiSeries - 1)/eps], ...]
expectZero["overall scale cancels", D[coeff, s]]
a5Pres = FullSimplify[a5 /. First[Solve[coeff == 0, a5, Reals]], ...]
```
Same intermediate variable name `coeff`, same `(... - 1)/eps`, same `d/ds` check, same Solve-for-a5 mechanic.

The mathematica .wl differs only in (a) one extra anchored check at L49 and (b) cosmetic naming (`a5Pres` vs `a5_sol`, `sNorm` vs `S`). This is a port, not an independent derivation. Codified second-engine policy says both engines must derive from physical premises independently.

**Why this matters:**
A line-by-line port catches typos and Mathematica/SymPy implementation differences but not derivation mistakes that are baked into the shared algebraic choreography. If both engines pick the wrong branch of an intermediate Solve, the transliterated version will agree without independent witness.

**Required change:**
Reshape the mathematica derivation to use a different algebraic path. Suggested independent route: instead of series-expanding the full ratio, expand the numerator and denominator separately to first order in `eps`, then form the ratio using `1/(3 + eps x) ≈ 1/3 - eps x/9 + O(eps^2)` long-division style. Specifically replace L36-41 of the .wl with something like:

```mathematica
num = Expand[3*(sNorm*beta^5 + 9*sigma5)];
den = Expand[3*sNorm - sigma0];
numLin = Normal[Series[num, {eps, 0, 1}]];
denLin = Normal[Series[den, {eps, 0, 1}]];
invDen = Normal[Series[1/denLin, {eps, 0, 1}]];
chiSeries = Expand[Normal[Series[numLin*invDen, {eps, 0, 1}]]];
expected = 1 + eps*(5*b + a0/3 + 9*a5);
expectZero["linearized chi law", chiSeries - expected];
```

Then derive the preservation condition by setting `chiSeries == 1` and solving for `a5` independently (instead of using the `coeff` intermediate that mirrors the sympy script). Keep L49 (the anchored closed-form check) as the verifier — it remains the strongest non-tautological check.

**Verification:**
After Codex applies, the .wl should no longer share the `coeff = (chiSeries - 1)/eps` intermediate. The new path should still print "chi_Q series", "a5 preservation condition", and pass the anchored check at the corresponding line. `redteam exec-mathematica 109` must exit 0 with the new derivation visible in the output.

### F3 — paper_misalignment

**Subtype:** script_missing_paper_claim

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_109.tex:21-25`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:1-44`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:1-55`

**What's wrong:**
The paper card's `\stagefield{Checks}` field (L21-25) enumerates three checks the stage is supposed to perform:

> "Check pure scale and pure argument deformations separately. Check Robin and standalone mixed-pole limits before imposing compensation. Check that the compensated branch preserves the even coefficients as well as the odd normalization."

The scripts cover:
- "pure scale" via `d coeff/ds` (sympy L38, mathematica L45) — yes
- "pure argument deformation separately" — NO (no `d coeff/db` to isolate the `5b` coefficient; only implicit via the full equality check A1)
- "Robin and standalone mixed-pole limits before imposing compensation" — NO (the scripts do not introduce any Robin parameter or mixed-pole pole structure; they work purely with the linearized `(S, beta, Sigma_0, Sigma_5)` ansatz)
- "even coefficients as well as odd normalization" preservation — NO (only the odd normalization defect `coeff` vanishing is asserted; even-coefficient preservation is built into the ansatz definition and not separately verified)

However, the notes file (the authoritative derivation source per the audit policy) does NOT mention Robin/mixed-pole limits or even-coefficient preservation. The notes' deliverables are exactly: the expansion box, the preservation condition, the closed-form `a_5`. The scripts cover these notes-side deliverables fully. So the script aligns with the notes but only partially with the .tex `Checks` list.

This is a paper-side question, not a unilateral script bug: either the card's `Checks` field is boilerplate-y and should be trimmed to match what the stage actually proves, or the stage's scope should expand to include the Robin/mixed-pole and even-coefficient checks (which would likely require importing Robin/mixed-pole parameters from upstream stages — outside this stage's declared `Inputs`).

**Why this matters:**
If left as-is, the card promises three secondary checks of which only one (pure scale) is actually performed. A future reader auditing the card against the scripts will flag this. Resolution requires a paper-side decision.

## Independent-derivation check (Mathematica)

See finding F2 above. The .wl mirrors the .py step-by-step: same five symbol introductions (`S`/`sNorm`, `beta`, `Sigma0`, `Sigma5` with identical structural form), same series-of-the-ratio approach, same intermediate variable `coeff = (chi_series - 1)/eps`, same `D[coeff, s]` scale-cancellation check, same `Solve[coeff == 0, a5]` mechanic. The one substantive difference is the anchored closed-form check at .wl L49 (which sympy lacks — that's F1). Otherwise this is a transliteration.

## Engine cross-check

Both engines run to exit 0. Output residuals match symbol-for-symbol:

- chi_Q series: sympy prints `a0*eps/3 + 9*a5*eps + 5*b*eps + 1`; mathematica prints `1 + (a0*eps)/3 + 9*a5*eps + 5*b*eps`. Same expression.
- first-order defect coefficient: sympy `a0/3 + 9*a5 + 5*b`; mathematica `a0/3 + 9*a5 + 5*b`. Match.
- a5 preservation condition: sympy `-a0/27 - 5*b/9`; mathematica `-1/27*a0 - (5*b)/9`. Match.
- preservation substitution: both `0`. Match.

No engine disagreement.

## Verdict justification

The bottom-line audit target ("first-order outgoing defect `Delta_Q = 5b + a_0/3 + 9 a_5`") is verified by both engines via a non-tautological assertion (A1/A4) that exercises the actual linearization of the closed-form `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3S - Sigma_0)`. The scale-cancellation implication and the preservation closed-form are both verified (the closed-form anchored only on the Mathematica side — F1). Engine outputs agree. The findings are: (F1) a tautological substitution on the sympy side that should be supplemented with an anchored closed-form check; (F2) the .wl is a structural transliteration of the .py; (F3) the paper card's `Checks` field promises secondary checks the scripts do not perform, though the notes do not echo those promises — a paper-side resolution is needed. None of these defeat the bottom-line claim; verdict is `findings`, severity uniformly low. No `stop_cold` — F1/F2 are mechanically applicable Codex fixes and F3 is a user-routed paper question.

## Self-test notes

Checked the proposed F1 anchor: `a5_sol - (-5*b/9 - a0/27)` is identically `0` for the current `coeff = a0/3 + 9*a5 + 5*b` (since solving `5b + a0/3 + 9a5 = 0` for `a5` yields `a5 = -5b/9 - a0/27`); a regression in any of the three coefficients would propagate to `a5_sol` and break this anchor, so it is non-tautological. Checked the proposed F2 path: numerator/denominator separate series in `eps` is well-defined since the denominator `3 - eps a0 + 3 eps s` is `3 + O(eps)` (no pole at `eps=0`); the rational `1/denLin` expanded to order 1 is `1/3 - (eps/9)(-a0 + 3s) + O(eps^2)`, and multiplying by `numLin = 3 + eps(... )` reproduces the same expansion via a different algebraic route. Parity/independence not at issue (no integrals or derivatives w.r.t. dummy variables here).
