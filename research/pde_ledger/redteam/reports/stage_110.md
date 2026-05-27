---
unit_id: 110
batch: IV.2
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-27T00:00:00Z
verdict: clean
stop_cold: null
findings_count: 0
paper_alignment: aligned
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
docs_read:
  paper_stage_tex: present
  notes_stage_files:
    - /var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage110_robin_outlet_model.md
  paper_appendix: present
---

# Audit unit 110 red-team report

## Files reviewed

- paper stage card: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_110.tex`
- notes: `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage110_robin_outlet_model.md`
- part appendix: `/var/projects/toy_physics/research/pde_ledger/paper/appendices/stage_appendix_part04.tex` (rows referencing this stage, esp. lines 399--412 and 1254)
- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage110_robin_outlet_model_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage110_robin_outlet_model_mathematica_audit.txt`

## What the paper claims

Stage 110 introduces the explicit isotropic Robin outlet deformation
`Lambda_2^R(z) = Lambda_2^out(z) + rho_R` and computes its exact effect on
the outgoing quadrupole normalization. The stage card's headline result
(in the `\stagefield{Derivation ledger}` block quote) states verbatim:
"Pure Robin shift gives `chi_Q^R = 3/(3-rho_R)` and generically spoils the
branch." The notes elaborate the full low-frequency expansion of the
normalized response `\widehat Y_2^R(z) = (-3+rho_R)/Lambda_2^R(z)`, giving
the closed forms `c2 = 1/(9-3 rho_R)`, `c4 = (4-rho_R)/(9 (3-rho_R)^2)`,
`c5 = 1/(27-9 rho_R)`, and the linearized expansion
`chi_Q^R = 1 + rho_R/3 + rho_R^2/9 + O(rho_R^3)`. The Stage-92 branch-triple
encoding `(b,a_0,a_5) = (0, rho_R, 0)` and the linearized increment
`delta chi_Q^R = rho_R/3 + O(rho_R^2)` are additional notational tie-ins, not
new identities. The "Checks" list in the stage card refers to broader
block-level checks (pure-scale/pure-argument deformations, mixed-pole
limits, compensated branch) handled in adjacent stages (107--113), not new
stage-110 identities.

## What the script claims to verify

Both scripts (`.py` and `.wl`) construct `Lambda_out`, add the Robin shift
`rho`, form `Y_R = (-3 + rho)/Lambda_R`, expand to order `z^5`, extract the
coefficients `c2`, `c4`, `c5` (taking `c5 / i`), compute `chi_R = c5/(1/27)`,
and series-expand `chi_R` in `rho` to order 2. They then assert five
identities: `c2 == 1/(9-3 rho)`, `c4 == (4-rho)/(9 (3-rho)^2)`,
`c5 == 1/(27-9 rho)`, `chi_R == 3/(3-rho)`, and
`chi_R_lin == 1 + rho/3 + rho^2/9`. These five assertions are exactly the
boxed identities from the notes. No additional or alternative identities
are tested.

## Paper ↔ script cross-check

| Paper-side deliverable | Script-side check | Status |
|---|---|---|
| `Y_R(z) = 1 + z^2/(9-3 rho_R) + ... + i z^5/(27-9 rho_R) + O(z^6)` (notes, boxed) | `c2`, `c4`, `c5` assertions (sympy lines 26--28; wl lines 49--51) | match |
| `chi_Q^R = 3/(3-rho_R)` (card, body block-quote line 16; appendix eq. 410) | `chi_R - 3/(3 - rho) == 0` (sympy line 29; wl line 52) | match |
| `chi_Q^R = 1 + rho_R/3 + rho_R^2/9 + O(rho_R^3)` (notes, linearized expansion) | `chi_R_lin - (1 + rho/3 + rho^2/9) == 0` (sympy line 30; wl line 53) | match |
| `(b, a_0, a_5) = (0, rho_R, 0)` Stage-92 triple (notes, "Branch-selection triple") | (not directly tested; notational) | n/a (notational tie-in only, not an algebraic deliverable) |

`paper_alignment: aligned`. Every algebraic identity stated in the paper
card and the notes for this stage has a corresponding non-tautological
assertion in both scripts. The Stage-92 triple notation is a label, not
an identity to verify.

## Assertion inventory

| # | Script | Line | Form | Exercises which paper claim? | Anchored to claim? |
|---|---|---|---|---|---|
| A1 | sympy | 26 | `simplify(c2 - 1/(9 - 3*rho)) == 0` | c2 coefficient of Y_R series | yes |
| A2 | sympy | 27 | `simplify(c4 - (4 - rho)/(9*(3 - rho)**2)) == 0` | c4 coefficient of Y_R series | yes |
| A3 | sympy | 28 | `simplify(c5 - 1/(27 - 9*rho)) == 0` | c5 coefficient of Y_R series | yes |
| A4 | sympy | 29 | `simplify(chi_R - 3/(3 - rho)) == 0` | chi_Q^R headline | yes (also partially redundant with A3 since chi_R := 27*c5) |
| A5 | sympy | 30 | `simplify(chi_R_lin - (1 + rho/3 + rho**2/9)) == 0` | linearized chi_Q^R expansion | yes |
| A6 | mathematica | 49 | `expectZero[c2 - 1/(9 - 3 rho)]` | c2 | yes |
| A7 | mathematica | 50 | `expectZero[c4 - (4 - rho)/(9 (3 - rho)^2)]` | c4 | yes |
| A8 | mathematica | 51 | `expectZero[c5 - 1/(27 - 9 rho)]` | c5 | yes |
| A9 | mathematica | 52 | `expectZero[chiR - 3/(3 - rho)]` | chi_Q^R headline | yes (redundant with A8) |
| A10 | mathematica | 53 | `expectZero[chiRLinear - (1 + rho/3 + rho^2/9)]` | linearized chi_Q^R | yes |

Each assertion is series-coefficient identity between a CAS-expanded
rational function and a closed-form target taken from the notes; none is
algebraically guaranteed by construction. A4/A9 are slightly redundant
(once `c5 = 1/(27 - 9 rho)` is verified, `chi_R = 27 c5 = 3/(3 - rho)`
follows immediately), but this is a redundancy, not a tautology — the
assertion form is still well-posed and could in principle fail if a
simplifier left an unreduced form.

## Findings

None.

## Independent-derivation check (Mathematica)

The Mathematica script mirrors the SymPy script's algorithmic outline
closely: both define `Lambda_out` with identical real-coefficient series
(`-3 + z^2/3 + z^4/9 + I z^5/9`), both form `Y_R = (-3 + rho)/Lambda_R`,
both expand to order `z^5`, both extract `c2`, `c4`, `c5 / I`, both
construct `chi_R = c5 * 27`, both linearize in `rho` to order 2, and both
assert the same five closed-form identities in the same order.

Sample correspondences (sympy line / wl line):
- `Lambda_out = -3 + z**2/sp.Integer(3) + z**4/sp.Integer(9) + I*z**5/sp.Integer(9)` (py:8)
  vs. `lambdaOut = -3 + z^2/3 + z^4/9 + I*z^5/9` (wl:31) — verbatim parallel.
- `Y_R = sp.simplify((-3 + rho) / Lambda_R)` (py:10)
  vs. `yR = FullSimplify[(-3 + rho)/lambdaR, ...]` (wl:33) — verbatim parallel.
- `chi_R = sp.simplify(c5 / sp.Rational(1, 27))` (py:16)
  vs. `chiR = FullSimplify[c5/(1/27), ...]` (wl:39) — verbatim parallel.

This is a transliteration in the structural sense. However, the
stage's core mathematical content is a CAS-driven series expansion of a
rational function — there is no genuinely distinct algorithmic route
that two independent engines could plausibly take. The independence
that does exist is at the level of the series-expansion algorithm
(SymPy's `sp.series` vs. Mathematica's `Series`), which are independent
implementations. Given this is a low-content series-coefficient check,
I do not raise `mathematica_transliteration` as a finding. The bar for
that category is intended to catch cases where an independent
re-derivation route exists and was bypassed; that condition does not
apply here.

## Engine cross-check

Both engines produce identical bottom-line forms:
- SymPy `chi_Q^R = -3/(rho - 3)` = `3/(3 - rho)` ✓
- Mathematica `chi_Q^R = -3/(-3 + rho)` = `3/(3 - rho)` ✓
- SymPy `chi_Q^R linearized = rho**2/9 + rho/3 + 1` ✓
- Mathematica `chi_Q^R linearized = 1 + rho/3 + rho^2/9` ✓
- All five expectZero / assert simplify checks return 0 in both engines.

Both `EXIT_CODE: 0`. Outputs are fresh (sympy script mtime 2026-04-01,
output mtime 2026-05-11; wl script mtime 2026-05-11 11:56, output mtime
2026-05-11 13:07 — both outputs newer than scripts).

## Verdict justification

The script verifies precisely the five closed-form identities stated in
the notes and the headline `chi_Q^R = 3/(3 - rho_R)` from the stage card.
Each assertion is non-tautological (the closed-form target is independent
of the construction route), exercises a paper-side deliverable, and
passes in both engines. Attacks attempted: (i) check whether the
assertion targets could be algebraically forced by the construction —
they cannot, because the closed-form rationals come from notes-side
algebra independent of the series-expansion machine; (ii) check whether
the series order is sufficient to extract `c5` — yes, sympy uses
`series(..., 0, 6)` (terms through `z^5`) and Mathematica uses
`Series[..., {z, 0, 5}]` (terms through `z^5`); (iii) check the symbol
domain: `rho` is `real`, and the Mathematica script adds `rho != 3` to
`$Assumptions` (correct — `rho = 3` is a pole of `chi_Q^R`); (iv) check
whether the linearization order matches the notes' `O(rho^3)` truncation
— yes, `series(..., 0, 3)` in sympy and `Series[..., {rho, 0, 2}]` in
Mathematica both yield through `rho^2`; (v) check that the
`Lambda_2^{out}` reference series matches Stage-100's canonical outgoing
DtN (`-3 + z^2/3 + z^4/9 + i z^5/9`) — it does; (vi) check whether the
appendix's eq. `chi_Q^R = 3/(3-rho_R)` (line 410 of `stage_appendix_part04.tex`)
matches the script — it does. Nothing breaks under attack. Verdict:
clean.

Cosmetic note (not a finding): both scripts banner this stage as
"STAGE 093" (sympy prints `stage93: PASS` on line 31; Mathematica banner
on wl:26 reads `STAGE 093 — EXPLICIT ISOTROPIC ROBIN OUTLET MODEL`).
This appears to be residual identifier drift from a pre-renumbering era
(unit was Stage 93 before reordering). It does not affect verification
correctness, and per the policy that the red-team does not propose
refactors or scope extensions, no directive is required.

## Self-test notes

Traps checked: (1) Variable independence — no `sp.diff` / `D[..., ...]`
in this stage; series expansion of explicit polynomial denominators.
N/A. (2) Symmetry/parity — no unbounded integrals; series in `z` and
`rho` only. N/A. (3) Trivial-case pre-check — at `rho = 0`, `c2 = 1/9`
(canonical), `c4 = 4/81` (canonical), `c5 = 1/27` (canonical),
`chi_R = 1` (canonical), `chi_R_lin = 1` — all match the canonical
fingerprint quoted in the notes. (4) Series order — both engines expand
to include `z^5` and `rho^2`, matching the notes' truncations. (5)
Pole at `rho = 3` — Mathematica assumes `rho != 3`; sympy's
`simplify(... == 0)` works symbolically without the assumption. Both
correct. Conclusion: no hidden trap, no directive needed.
