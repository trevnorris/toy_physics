---
title: U1 throat-body dynamics program
type: source
status: current
sources:
  - software/em_charge_attribute/reports/u1_body_dynamics.md
  - software/em_charge_attribute/reports/u1_body_dynamics_results.yaml
  - software/em_charge_attribute/run_u1_body_dynamics.sh
  - software/em_charge_attribute/run_u1_body_b2.sh
  - software/em_charge_attribute/run_u1_body_phaseC.sh
  - software/em_charge_attribute/u1_body_dynamics_sympy.py
  - software/em_charge_attribute/u1_body_dynamics_dual.wl
  - software/em_charge_attribute/u1_body_dynamics_compare.py
last_updated: 2026-09-03
---

# U1 throat-body dynamics program

## Purpose and current result

This program builds an effective collective-coordinate description of a
moving throat after retiring the independent brane polar field. It derives an
exterior co-moving family, classifies radial tails, solves five distinct
endpoint response laws, and reduces the field action to throat/body moments.
It is an effective-action calculation joined to open core traces, not a
nonlinear throat simulation.

The current report contains several phases, so its early “halt after Phase A”
language is historical within the same appended document. Read together, the
reported state is:

- **Phase A:** `COMPUTATION_VALID` and `U1_BASE_OK` on the declared positive
  coefficient stratum for all five endpoint cells.
- **Phase B1:** the indexed mechanics machinery runs, but all ten
  endpoint/ambient mechanics cells remain `UNRESOLVED` because required
  profile and open-root data are absent.
- **Phase B2:** `PASS_WITH_HONEST_OUTCOMES`, while return closure, the intake
  coefficient, and radiative residues remain unresolved.
- **Phase C:** not run.

Sources:

- `software/em_charge_attribute/reports/u1_body_dynamics.md` — headings `Verdict and Phase-A halt`, `Phase B1 — indexed mechanics remediation 3`, and `Phase B2 — Intake response and radiative residues`

## Established conditional structure

Phase A constructs co-moving continuity and momentum balances before reducing
them, derives the native control-surface force, solves six exterior channels,
and finds all six translated tails normalizable in its ambient-subtracted
scheme. It also constructs the linearized force-balance operator, translation
zero mode, five boundary/constraint response maps, and the reduced effective
Lagrangian coefficients.

The endpoint choices are genuinely distinct: no slip, free slip, permeable
texture, nonholonomic shear lock, and Robin/Rayleigh response. The scripts
exercise unstable, unresolved, and non-normalizable control branches, so
`U1_BASE_OK` is conditional rather than hardcoded.

Phase B1 carries a correction to the Phase-A brane-shear unit contribution,
changing $12\pi/5$ to $8\pi/3$ in the governed `U_XX -> G_VV` chain. The
tilt-profile rows remain explicitly unresolved. It also reports a passing G4
winding control, but no full mechanics signature is promoted from the missing
coefficients.

Sources:

- `software/em_charge_attribute/reports/u1_body_dynamics.md` — headings `Computed co-moving laws`, `far-field solve that decides G1`, `Endpoint response solves`, `Evaluated reduced moments and L_eff`, and `Phase-A amendment carried into B1`
- `software/em_charge_attribute/u1_body_dynamics_sympy.py` — Phase-A tail, balance, endpoint, and classification pipeline

## Evidence and operational notes

The Phase-A SymPy and Wolfram engines produce structured artifacts that the
comparator checks at the level of source-derived expressions, ODE residuals,
norms, endpoint coefficients, reduced monomials, and dependency sets. The
committed YAML is the aggregate machine-readable result and is large; agents
should query only the relevant phase or verdict fields rather than load it for
ordinary retrieval.

`run_u1_body_dynamics.sh` is a staged Phase-A/B1 runner. The B2 runner is a
heavy, environment-traced workflow with explicit resume contracts, while
`run_u1_body_phaseC.sh` dispatches the Phase-C stage runner. These are
specialized audit programs, not routine memory-update commands.

Sources:

- `software/em_charge_attribute/run_u1_body_dynamics.sh` — functions `run_stage1` and `run_stage2`
- `software/em_charge_attribute/run_u1_body_b2.sh` — stage/resume dispatch and isolated execution helpers
- `software/em_charge_attribute/u1_body_dynamics_compare.py` — functions `build_report`, `build_results`, and `main`

## Open questions

Core and surface functionals, return closure, several profile tangents, and the
open $M_h,c_E$ coefficients remain inputs rather than derived outputs. The
program does not yet determine a complete body mechanics tensor, radiation,
or Phase-C observable. Passing its executable gates shows consistency of the
encoded reductions and bookkeeping; it is not evidence that the toy model is
a physical U(1) theory.

Sources:

- `software/em_charge_attribute/reports/u1_body_dynamics.md` — headings `Frozen setup and honest scope`, `Derived positive content and honest exits`, and `Records, gates, and halt`
