---
title: Projection and reduction
type: topic
status: current
sources:
  - research/4d/paper/4d.tex
  - research/4d_em_fields/paper/4d_em_fields.tex
  - research/4d_plasma/paper/4d_plasma.tex
  - research/pde/paper/pde.tex
  - software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md
last_updated: 2026-09-03
---

# Projection and reduction

## The central firewall

Projection and dimensional reduction are different operations.

- The normalized kernel $W(w)$ defines what a brane observer measures by
  averaging a bulk field across the transverse direction.
- The profile $Z(w)$ localizes the Maxwell action and sets the zero-mode
  normalization $Z_{\rm int}=\int Z(w)\,dw$.
- A dimensional reduction additionally imposes a mode ansatz and integrates
  over $w$ to obtain effective fields and couplings.

Using $W$ where $Z$ is required, or treating a projected field as an automatic
zero mode, changes the question being answered.

Sources:

- `research/4d/paper/4d.tex` — labels `sec:braneobs`, `sec:emreduction`, and `eq:qeff-canonical`
- `research/4d_plasma/paper/4d_plasma.tex` — labels `sec:projection` and `sec:localization`
- `research/pde/paper/pde.tex` — labels `sec:status-parent` and `eq:zero-mode-assumptions`

## Projection leaves an open subsystem

Projecting exact bulk continuity gives a brane balance law with a leakage term.
That term contains both transverse boundary flux and a $W'(w)j^w$ weighting
piece. It vanishes only under additional falloff/no-leak conditions. Projected
plasma energy and helicity similarly require transverse reservoirs or fluxes
to close their ledgers.

The brane should therefore be read as an open subsystem unless a source proves
the relevant transverse exchange negligible. A familiar three-dimensional
continuity, energy, or MHD equation is a controlled limit rather than the raw
projection.

Sources:

- `research/4d/paper/4d.tex` — labels `eq:proj-continuity`, `eq:Sleak-general`, and `eq:boundary-vanish`
- `research/4d_plasma/paper/4d_plasma.tex` — labels `eq:projected-continuity`, `eq:projected-energy-balance`, and `eq:projected-helicity-balance`

## Maxwell zero mode and mixed channels

Under $A_w=0$, $\partial_wA_\mu=0$, Lorenz gauge, and suppressed transverse
current, the localized Maxwell system reduces to a brane Maxwell equation with
effective coupling $\mu_0/Z_{\rm int}$. Source profiles must also be compatible
with $Z(w)$ or higher transverse modes are excited. Gaussian localization gives
a massless zero mode plus a gapped tower and a Coulomb potential with Yukawa
corrections.

This reduction suppresses $A_w$, $J^w$, $F_{\mu w}$ and the mixed invariants
$E_w,C_a$ in that regime; it does not erase them from the parent ontology.
Near a charged throat they may carry precisely the information that the pure
zero mode omits.

Sources:

- `research/4d/paper/4d.tex` — labels `eq:zero-mode`, `eq:brane-maxwell`, and `eq:source_profile_matching`
- `research/4d_em_fields/paper/4d_em_fields.tex` — labels `eq:brane-maxwell`, `eq:gaussian-spectrum`, and `eq:coulomb-yukawa`
- `research/pde/paper/pde.tex` — labels `eq:mixed-invariants`, `eq:zero-mode-assumptions`, and `eq:effective-maxwell`

## Relation to the brane-elastic fork

The Stage-0 brane-elastic action is classified as a `NEW_PARENT_ACTION`, not a
reduction or field redefinition of the localized Maxwell action. Its Cauchy,
rotational, and Cosserat constitutive options are postulated candidates. If a
future model keeps both the localized Maxwell carrier and a brane-elastic EM
carrier, it must explicitly resolve double counting.

Source:

- `software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md` — headings `VERDICT`, `Explicit action`, and `Constitutive-form menu and DOF audit`

## Reading rule

Before using a three-dimensional result, ask three questions: which kernel
defines the observable, which mode/profile ansatz performs the reduction, and
which boundary or leakage terms were dropped. Without those answers, a brane
formula is not yet a closed effective theory.

Related pages:

- `memory/topics/foundations.md`
- `memory/topics/charge-and-electromagnetism.md`
- `memory/topics/status-and-reading-rules.md`
