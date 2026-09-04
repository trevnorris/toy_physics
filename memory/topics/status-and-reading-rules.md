---
title: Status and reading rules
type: topic
status: current
sources:
  - research/4d/paper/4d.tex
  - research/4d_1pn_full/paper/4d_1pn_full.tex
  - research/4d_2pn/paper/4d_2pn.tex
  - research/4d_em_fields/paper/4d_em_fields.tex
  - research/4d_plasma/paper/4d_plasma.tex
  - research/4d_2_5pn/paper/4d_2_5pn.tex
  - research/4d_4pn/paper/4d_4pn.tex
  - research/pde/paper/pde.tex
  - research/pde_audit/CLAIM_CHECK_MATRIX.md
  - research/pde_audit/notes/stage_v2_16_branch_freeze_no_refit_derivation.md
  - research/pde_audit/notes/stage_v2_26_program_status_after_audit.md
  - research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md
  - software/stage1_solver/STAGE1_VERDICT.md
last_updated: 2026-09-03
---

# Status and reading rules

## Use the weakest justified status

The framework paper defines the following claim ladder. When more than one
label could apply, use the weaker one.

| Status | Meaning in this repository |
|---|---|
| Exact | Follows from the declared parent action or definitions without a reduction. |
| Exact Within Closure | Follows exactly after a named finite ansatz, constitutive law, or other closure is imposed. |
| Reduced / Controlled Reduction | Requires a stated scale hierarchy, projection, truncation, or asymptotic regime. |
| Effective Closure | A useful modeled relation that is not derived from the strict parent. |
| Numerically Located | Found for a specified branch, discretization, and frozen input packet. |
| Open | Not established, not selected, or still awaiting a realized physical branch. |

Source: `research/pde/paper/pde.tex` — table `tab:claim-status` and its
immediately following reading rule.

The `status: current` field at the top of a memory page means only that the
page reflects the current repository. It does **not** upgrade every scientific
statement on that page to an exact or established result.

## Legacy capsule tags

Some source capsules retain extraction-era lines such as
`lifecycle=current`, `evidence=derived`, and `memory_review=ai_draft`. They are
not a competing status ladder:

| Legacy tag | How to read it now |
|---|---|
| `lifecycle=current` | The capsule did not identify the statement as superseded. This says nothing about whether the statement is exact, conditional, or physically realized. |
| `evidence=derived` | The cited source presents a derivation. Classify its result as Exact, Exact Within Closure, or Reduced according to the assumptions stated in the prose; never promote it from this tag alone. |
| `evidence=provisional` | The statement is adopted, conditional, interpretive, or awaiting stronger support. Use the prose to choose the weaker current-ladder label. |
| `evidence=open` | Read as Open. |
| `memory_review=ai_draft` | An editorial note about the original AI extraction, not scientific evidence and not a freshness indicator. |

When an old inline tag, the surrounding prose, and a current topic page differ,
the explicit scope and current topic synthesis control. Follow their citation
to the original source before relying on a technical detail.

## Boundaries that must remain visible

### Projection is not reduction

A brane observable such as
$W(\mathbf x,t)=\int dw\,\mathcal W(w)f(\mathbf x,w,t)$ is a weighted
observation of the parent field. Replacing the parent field by a zero mode,
discarding mixed channels, or imposing quasi-static behavior is a later
reduction. Projected leakage makes the brane subsystem open; projection does
not erase it.

Sources:

- `research/4d/paper/4d.tex` — sections `sec:braneobs` and `sec:emreduction`
- `research/4d_plasma/paper/4d_plasma.tex` — equation `eq:W_matched`

The matched plasma weight $W=Z$ is a special calibrated equality, not an
identity between the observation kernel $W$ and the localization profile $Z$.
Likewise, a zero-mode limit can suppress $A_w,F_{\mu w},J^w$ contributions but
does not remove them from the parent theory.

### An exact identity is not its familiar limiting equation

The longitudinal projection identity is exact once its definitions are fixed.
The Poisson equation and inverse-square regime require quasi-static,
longitudinal, nearly constant-density, and localized-source assumptions plus
small correction terms. Any identification of brane velocity with test-body
gravitational acceleration is an additional constitutive step.

Source: `research/4d/paper/4d.tex` — equation
`eq:longitudinal-identity`, section `sec:poisson-regime`, and equation
`eq:poisson-approx`.

Similarly, integrating the Maxwell equation over $w$ is exact, while standard
brane Maxwell theory requires boundary, zero-mode, current, and gauge
conditions. It was included in the parent action; it has not thereby emerged
from native throat variables.

Source: `research/4d_em_fields/paper/4d_em_fields.tex` — equations
`eq:integrated_exact` through `eq:brane_maxwell_gf`.

### Algebraic PN matching is conditional

The 1PN--4PN matches are exact coefficient or Legendre-transform statements
inside declared response and chart closures. They do not prove that the
moving-throat PDE realizes those inputs. The 2.5PN outgoing quadrupole
normalization remains open, and 4PN inherits the same datum.

Sources:

- `research/4d_1pn_full/paper/4d_1pn_full.tex` — sections `sec:taxonomy` and `sec:kappaPV-closure`
- `research/4d_2pn/paper/4d_2pn.tex` — section `sec:status-open`
- `research/4d_2_5pn/paper/4d_2_5pn.tex` — section `sec:conditional-theorem-gap-normalization`
- `research/4d_4pn/paper/4d_4pn.tex` — section `sec:fixed-open-still-open`

### A script pass is not a physical result

Symbolic and numerical harnesses test the equations, fixtures, tolerances, and
contracts they encode. A successful exit can establish internal consistency;
it cannot establish a missing closure, physical phase, boundary class, or
realized branch. Read the generated report or verdict, not only the exit code.

Sources:

- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — claim classifications and named limitations
- `software/stage1_solver/STAGE1_VERDICT.md` — sections `The result — a robust MISS` and `Interpretation — what it means (and does NOT mean)`

## Honest branch testing

Branch tests must be target-blind:

1. declare the parent/effective action, material and boundary class, source
   law, support, stability gates, and extraction definitions;
2. solve the branch;
3. freeze those choices before revealing the target residual;
4. export all compared observables from that same branch;
5. do not refit after the comparison.

Sources:

- `research/pde_audit/notes/stage_v2_16_branch_freeze_no_refit_derivation.md` — heading `Why the no-refit rule is mandatory`
- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — headings `Frozen Inputs` and `Non-Rescue Rules`

A miss then rejects the frozen branch or tested family under its assumptions.
It does not overturn exact parent identities, and it does not prove that every
possible completion fails. Conversely, a pass by an effective or manufactured
fixture would validate that fixture only, not the existence of a physical
branch.

## Recommended reading order

For any claim, ask in order:

1. What is the exact parent identity or declared action?
2. Which projection or controlled reduction is being applied?
3. Which effective or constitutive closure supplies the remaining inputs?
4. Has a self-consistent branch actually been located and frozen?
5. What remains open after the comparison?

That sequence is the compact replacement for the old Atlas evidence levels.
The topic page should give the current synthesis; its cited source remains
authoritative for technical detail.
