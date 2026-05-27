---
unit_id: 107
batch: IV.2
auditor_model: claude-opus-4-7[1m]
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
  notes_stage_files: [moving_throat_pde_stage107_general_dtn_deformation.md]
  paper_appendix: present
---

# Audit unit 107 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_107.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage107_general_dtn_deformation.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows mentioning Stages 107--113 / `\input{stages/stage_107}`)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.txt`

## What the paper claims

Stage 107 introduces the first explicit isotropic moving-throat DtN deformation family
`Lambda_2^def(z) = S*Lambda_2^out(beta z) + Sigma_0 + Sigma_2 z^2 + Sigma_4 z^4 + i Sigma_5 z^5` and
derives the exact map from deformation data `(S, beta, Sigma_0, Sigma_2, Sigma_4, Sigma_5)` to
the retarded quadrupole normalization `chi_Q`. The paper card quote is the bottom line:
"Deformed branch `S Lambda_2^out(beta z) + Sigma_0 + Sigma_2 z^2 + Sigma_4 z^4 + i Sigma_5 z^5`
gives explicit `chi_Q`."  The notes spell out three "boxed" deliverables:
(1) the normalized expansion `Y_2^def(z) = 1 - (L_2/L_0)z^2 + (L_2^2/L_0^2 - L_4/L_0)z^4 - i(L_5/L_0)z^5`,
(2) the canonical-even matching closed forms `Sigma_2 = -(3 S beta^2 - 3 S + Sigma_0)/9` and
`Sigma_4 = -(3 S beta^4 - 3 S + Sigma_0)/27`, and
(3) the post-even-match odd normalization `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)`.
The paper card "Checks" list also asks for "Check pure scale and pure argument deformations
separately" and "Check Robin and standalone mixed-pole limits before imposing compensation",
but those particular subchecks live further downstream (Stage 108+); Stage 107 is the algebra
backbone.

## What the script claims to verify

Both scripts construct `Lambda_def`, extract `(L_0, L_2, L_4, L_5)`, define the candidate
`Y_formula`, and assert it agrees with the SymPy `series` / Mathematica `Series` of
`L_0/Lambda_def` through order `z^5`. Both then solve the two even-matching equations
`-L_2/L_0 = 1/9`, `L_2^2/L_0^2 - L_4/L_0 = 4/81` for `(Sigma_2, Sigma_4)`, substitute into
`chi_Q = (-L_5/L_0)/(1/27)`, and assert the resulting expression equals the paper's boxed
`3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)`.  Mathematica additionally asserts that the solved
`Sigma_2` and `Sigma_4` equal the paper's boxed closed forms; SymPy only prints them.

## Paper -> script cross-check

| Paper-side deliverable | Sympy check | Mathematica check | Status |
|---|---|---|---|
| Box 1: `Y_2^def(z) = 1 - (L_2/L_0)z^2 + (L_2^2/L_0^2 - L_4/L_0)z^4 - i(L_5/L_0)z^5` | line 44 (`expect_zero('normalized expansion direct-formula', ...)`) | line 52 (`expectZero["normalized expansion direct-formula", ...]`) | match |
| Box 2a: `Sigma_2 = -(3 S beta^2 - 3 S + Sigma_0)/9` | line 58 prints only, no assertion | line 68 (`expectZero["Sigma2 exact formula", ...]`) | partial (Mathematica covers; SymPy only prints) |
| Box 2b: `Sigma_4 = -(3 S beta^4 - 3 S + Sigma_0)/27` | line 59 prints only, no assertion | line 69 (`expectZero["Sigma4 exact formula", ...]`) | partial (Mathematica covers; SymPy only prints) |
| Box 3: `chi_Q = 3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)` | lines 63-66 (`expect_zero(... chi_even - 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0))`) | lines 73-76 (`expectZero["chi_Q - ...", ...]`) | match |

Note on Box 3: `chi_Q = (-L_5/L_0)/(1/27)` does not algebraically depend on `Sigma_2` or
`Sigma_4`, so the substitution `chi_Q.subs(sol)` is a no-op and the assertion is structurally
the same as `(-L_5/L_0)/(1/27) == 3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)`. This still
exercises a non-trivial algebraic identity (the `-27 L_5/L_0` simplification), so it is not
tautological -- but the label "chi_Q under canonical-even matching" overstates what the
substitution did. This is informational, not a finding, because the paper itself states the
formula in the same way: "with [even matching] imposed, the exact odd normalization becomes
chi_Q = ..." -- which is mathematically equivalent to "chi_Q is independent of
(Sigma_2, Sigma_4) and equals this".

Set `paper_alignment: partial` because Box 2a/2b are only covered by one of the two engines.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 44 | `expect_zero('normalized expansion direct-formula', sp.expand(Y_direct - Y_formula))` | Box 1 (`Y_2^def` low-`z` form) | yes |
| A2 | sympy | 55-56 | `if len(sol) != 1: raise AssertionError(...)` | structural (unique even-match solution) | yes (uniqueness) |
| A3 | sympy | 63-66 | `expect_zero('chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0)', chi_even - 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0))` | Box 3 (`chi_Q` closed form) | yes |
| A4 | mathematica | 52 | `expectZero["normalized expansion direct-formula", yDirect - yFormula]` | Box 1 | yes |
| A5 | mathematica | 63 | `If[Length[sol] =!= 1, fail[...]]` | uniqueness | yes |
| A6 | mathematica | 68 | `expectZero["Sigma2 exact formula", (sigma2 /. sol) + (3*sNorm*beta^2 - 3*sNorm + sigma0)/9]` | Box 2a (`Sigma_2` closed form) | yes |
| A7 | mathematica | 69 | `expectZero["Sigma4 exact formula", (sigma4 /. sol) + (3*sNorm*beta^4 - 3*sNorm + sigma0)/27]` | Box 2b (`Sigma_4` closed form) | yes |
| A8 | mathematica | 73-76 | `expectZero["chi_Q - 3(sNorm beta^5 + 9 sigma5)/(3 sNorm - sigma0)", chiEven - 3*(sNorm*beta^5 + 9*sigma5)/(3*sNorm - sigma0)]` | Box 3 | yes |

The SymPy script has no analog of A6 / A7: it prints `Sigma2_evenmatch` and
`Sigma4_evenmatch` at lines 58-59 but never compares them to the paper's boxed forms.

## Findings

### F1 — insufficient_verification

**Severity:** medium

**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py:54-71`

**What's wrong:**

Notes box (Box 2a/2b) verbatim:
```
\Sigma_2 = -\frac{3S\beta^2 - 3S + \Sigma_0}{9},
\qquad
\Sigma_4 = -\frac{3S\beta^4 - 3S + \Sigma_0}{27}.
```
SymPy script lines 58-59:
```python
print('Sigma2_evenmatch =', sp.simplify(sol[Sigma2]))
print('Sigma4_evenmatch =', sp.simplify(sol[Sigma4]))
```
These two boxed identities are computed and printed but never wrapped in an
`expect_zero(...)` (or any other) assertion. The only assertion downstream of the solve is
`expect_zero('chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0)', chi_even - 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0))` at lines 63-66 -- but `chi_Q` does not depend on
`Sigma_2` or `Sigma_4`, so this assertion is independent of the closed-form values of
`Sigma_2` and `Sigma_4` returned by `sp.solve`. A regression that returned a wrong
`Sigma_2`/`Sigma_4` from `solve` would not fail any SymPy assertion in this script.

The Mathematica script does check both at lines 68-69 (`expectZero["Sigma2 exact formula", ...]` and `expectZero["Sigma4 exact formula", ...]`), so the second-engine check exists -- but the
two-engine policy expects both engines to lock in the paper's boxed deliverables, not just one.

**Why this matters:**

Two of the paper's three boxed Stage 107 deliverables are silent in SymPy. If `sp.solve`
returned a different branch (e.g. a SymPy version change in how `solve` orders or selects
roots), the script could still pass while reporting nonsense for `Sigma2_evenmatch` /
`Sigma4_evenmatch`. The transcript visually "shows" the right answer, but no assertion
guards it.

**Required change:**

After lines 58-59, add two `expect_zero` checks that assert the SymPy-side solved
`Sigma_2`/`Sigma_4` match the paper's boxed closed forms. Concretely (between the existing
`print('Sigma4_evenmatch ...')` line and the `chi_even = ...` line, i.e. between lines 59
and 61):

```python
expect_zero(
    'Sigma2 exact formula',
    sol[Sigma2] - (-(3*S*beta**2 - 3*S + Sigma0))/sp.Integer(9),
)
expect_zero(
    'Sigma4 exact formula',
    sol[Sigma4] - (-(3*S*beta**4 - 3*S + Sigma0))/sp.Integer(27),
)
```

Note the explicit `sp.Integer(9)` / `sp.Integer(27)` so the division stays symbolic. The
paper's notes write the formulas as `-(3 S beta^2 - 3 S + Sigma_0)/9` and
`-(3 S beta^4 - 3 S + Sigma_0)/27`; the script's printed values
(`-S*beta**2/3 + S/3 - Sigma_0/9`, `-S*beta**4/9 + S/9 - Sigma_0/27`) are the same numbers
once expanded, so `expect_zero(..., sol[Sigma2] - <paper form>)` will pass when the math is
right.

**Verification:**

After Codex applies, re-running the SymPy script should print two new lines:
```
Sigma2 exact formula = 0
Sigma4 exact formula = 0
```
and the existing `chi_Q - 3(S beta^5 + 9 Sigma5)/(3S - Sigma0) = 0` line should still appear.
The script must still exit 0.

## Independent-derivation check (Mathematica)

The `.wl` script follows the same algorithmic shape as the `.py`: build `Lambda_def`, extract
`L_0..L_5` via coefficient extraction, define `yFormula` symbolically, compare against
`Series[l0/lambdaDef, ...]`, solve the same two even-matching equations, substitute, and
assert the final `chi_Q` matches the paper's boxed form. Lines 33-34 of the `.wl` and lines
27-30 of the `.py` are structurally parallel, and the `yFormula` definitions at `.wl:42` and
`.py:43` are character-for-character the same identity.

However, this is the natural way to do the derivation: there is only one obvious algorithm
once the deformation family is written down, and both scripts compute the same coefficients
in the same order because the math dictates that order. The Mathematica script does use its
own primitives (`Coefficient[..., z, 2]`, `Series[..., {z, 0, 5}]`, `Solve`) rather than
calling out to the SymPy results. It also adds two assertions (`Sigma2 exact formula`,
`Sigma4 exact formula`) the SymPy script lacks, which would not appear in a mechanical
transliteration. I do not flag `mathematica_transliteration`: the parallel structure is
explained by the mathematics, not by porting choices.

## Engine cross-check

Both engines report the same final closed forms:

SymPy output:
```
L0 = -3*S + Sigma0
L2 = S*beta**2/3 + Sigma2
L4 = S*beta**4/9 + Sigma4
L5 = S*beta**5/9 + Sigma5
chi_Q = 3*(S*beta**5 + 9*Sigma5)/(3*S - Sigma0)
Sigma2_evenmatch = -S*beta**2/3 + S/3 - Sigma0/9
Sigma4_evenmatch = -S*beta**4/9 + S/9 - Sigma0/27
```

Mathematica output:
```
L0 = sigma0 - 3*sNorm
L2 = sigma2 + (beta^2*sNorm)/3
L4 = sigma4 + (beta^4*sNorm)/9
L5 = sigma5 + (beta^5*sNorm)/9
chi_Q = (-3*(9*sigma5 + beta^5*sNorm))/(sigma0 - 3*sNorm)        # == 3*(sNorm*beta^5 + 9*sigma5)/(3*sNorm - sigma0)
Sigma2_evenmatch = (-sigma0 + 3*sNorm - 3*beta^2*sNorm)/9        # == -S*beta^2/3 + S/3 - Sigma_0/9
Sigma4_evenmatch = (4*sigma0^2 - 24*sigma0*sNorm + ...)/(-81*sigma0 + 243*sNorm)
                                                                 # FullSimplify lands on the same value;
                                                                 # `expectZero["Sigma4 exact formula", ...]` PASSed
```

Both `chi_Q - 3(... beta^5 + 9 Sigma_5)/(3 ... - Sigma_0) = 0` checks pass in both engines.
Mathematica's `Sigma4_evenmatch` is left in an unsimplified form but the subsequent
`expectZero["Sigma4 exact formula", ...]` reduces it to 0. Engines agree.

## Verdict justification

The math is correct and both engines reach the same closed forms for `Y_2^def`, the boxed
even-matching `Sigma_2`/`Sigma_4`, and the final `chi_Q`. Mathematica explicitly anchors all
three boxed paper deliverables; SymPy anchors `Y_2^def` and `chi_Q` but only prints (does not
assert) the boxed `Sigma_2` and `Sigma_4` forms. That is a coverage gap, not a math error,
so the verdict is `findings` with one medium-severity `insufficient_verification` item; no
stop-cold flag.

Attacks attempted that failed: (1) check whether the `chi_Q` assertion is tautological after
the `Sigma_2`/`Sigma_4` substitution -- it is structurally equivalent to the unsubstituted
form because `chi_Q` is `Sigma_2`/`Sigma_4`-independent, but the simplification
`(-L_5/L_0)/(1/27) -> 3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)` is still non-trivial algebra;
(2) check whether `sp.im(coeff(Lambda_def, z, 5))` could pick up spurious real parts -- it
cannot, because all symbols are declared `real=True` and the `z^5` term is built from `I *
(real)`; (3) verify that the boxed `chi_Q` formula matches the notes verbatim modulo
factor-of-3 conventions -- it does (`3(S beta^5 + 9 Sigma_5)/(3 S - Sigma_0)`); (4) compare
`L_0..L_5` against the notes' explicit values -- exact match.

## Self-test notes

- Variable independence: the proposed new SymPy assertions for `Sigma_2`/`Sigma_4` reference
  `sol[Sigma2]` and `sol[Sigma4]`, both of which actually depend on `S`, `beta`, `Sigma_0`
  via the solved expressions. No trivial-zero-derivative trap.
- Symmetry/parity: not applicable here (no integrals; pure algebra over six scalar symbols).
- Trivial-case pre-check: substituting `S=1, beta=1, Sigma_0=0` into the proposed paper-form
  RHS gives `Sigma_2 = -(3 - 3 + 0)/9 = 0` and `Sigma_4 = -(3 - 3 + 0)/27 = 0`; substituting
  the same into the script's `Sigma2_evenmatch = -1/3 + 1/3 - 0 = 0` and `Sigma4_evenmatch = -1/9 + 1/9 - 0 = 0`. Residual zero, as required.
- Paper round-trip: the proposed fix only adds assertions; it introduces no new constants and
  no new derivations. No new `paper_misalignment` risk.
- Paths: the only fix target is the SymPy `.py` under `scripts/`; the `.wl` already has the
  Mathematica side covered.
