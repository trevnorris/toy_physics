---
title: Moving-throat dynamics
type: topic
status: current
sources:
  - research/4d/paper/4d.tex
  - research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
  - research/pde/paper/pde.tex
  - research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md
  - research/pde_audit/notes/stage_v2_25_actual_branch_protocol_derivation.md
  - research/pde_audit/notes/stage_v2_26_program_status_after_audit.md
  - research/pde_audit/notes/stage_v2_28_physical_picture_and_ontology.md
  - research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md
  - research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md
  - software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md
  - software/stage1_solver/derivations/pathA_01_return_source_and_balance.md
  - software/stage1_solver/STAGE1_VERDICT.md
last_updated: 2026-09-03
---

# Moving-throat dynamics

## Current bottom line

The program has an exact matter/gauge parent, a distributed throat geometry,
and a consistent **linear effective** wall–support–gauge response framework.
It does not yet have a uniquely derived nonlinear autonomous throat action or
a physical exporter that solves the complete branch and produces all required
observables. Current reduced and frozen candidate families miss their target;
that rejects those candidates, not the framework's exact identities.

## Physical picture

The working defect is a finite-radius brane puncture and open conduit with a
mouth at $w=0$ and a finite bulk interior. Mouth data, interior support, radial
throughput, projected leakage, endpoint work, and outgoing radiation are
different pieces of the problem. A hard-capped pocket is not the operative
branch picture.

The distributed geometry is
$\Sigma(\mathbf X,t)=r-R(\Omega,w,t)$. The old $a(t),L(t)$ variables survive
as collective moments, while real $\ell=2$ wall/support modes become literal
degrees of freedom. Reduced quantities such as $D_0,C,P_0,N_2,N_4$ are
response readouts, not physical parts of the throat.

Sources:

- `research/pde/paper/pde.tex` — labels `eq:geometry-lift-levelset`, `eq:geometry-lift-a-moment`, `eq:geometry-lift-L-moment`, and `eq:geometry-lift-harmonic-decomp`
- `research/pde_audit/notes/stage_v2_28_physical_picture_and_ontology.md` — headings `Core ontology`, `Mouth versus interior`, `Open system and superfluid intake/output`, and `What the response coefficients do not describe`

## Three action levels

### 1. Strict parent

The strict parent is GNLS matter plus localized Maxwell, with geometry entering
the confinement potential. Varying confinement supplies a force on a declared
wall displacement; by itself it does not supply wall inertia, stiffness, or an
autonomous wall PDE.

### 2. Effective linear wall closure

After adopting the quadratic distributed wall action $S_\eta^{(2)}$, the
linear wall equation, modal split, boundary terms, positivity conditions, and
coupling to BdG/gauge response are consistent within that closure. Schur
elimination gives exact kernels for a selected finite mode list, not for the
unsolved nonlinear parent.

### 3. Promoted nonlinear throat action

A complete parent would add and freeze a nonlinear $S_\Sigma[R]$ or an
equivalent material action. Its constitutive functions and boundary laws are
not yet derived or uniquely selected. The Stage-0 brane-elastic Cauchy,
rotational, and Cosserat families are `NEW_PARENT_ACTION` candidates, not a
completed promotion. They also carry an unresolved double-counting question
if retained beside the existing localized Maxwell carrier.

Sources:

- `research/pde/paper/pde.tex` — labels `sec:status-parent`, `sec:linearized-pde`, `eq:linpde-wall-action`, and `sec:reduced-system`
- `research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md` — headings `Executive result`, `Variation of the confinement-only parent term`, and `Variation after adding the quadratic wall action`
- `software/stage1_solver/reports/pathA_23_stage0_action_and_contracts.md` — headings `VERDICT`, `Candidate in-plane constitutive menu`, and `Open questions handed to later stages`

## Return coupling

At quadratic order the declared confinement interaction fixes the reciprocal
matter return source

$$
S_\eta^{(\psi)}=-k_1\,\delta\rho,
$$

while the $k_2\rho_0\eta$ term renormalizes wall stiffness. The fixed
localization Maxwell action has no direct $\eta\to Z$ wall source, so gauge
return is matter-mediated and remains an on-shell functional until the coupled
matter/gauge problem is solved. Closing this return loop changes the realized
background and denominator-side response; it does not add a free numerator
knob in the stated quadratic derivation.

Source:

- `software/stage1_solver/derivations/pathA_01_return_source_and_balance.md` — headings `D1: Matter Return Source`, `D2: Gauge Return Source`, `D3: Static Self-Consistent Throat Balance`, and `D4: No New Numerator Magnitude`

## Branch workflow

An honest comparison follows this order:

1. choose the parent action, material/boundary class, source law, support
   basis, and stability gates;
2. solve one self-consistent branch;
3. freeze that branch and its extraction definitions before exposing target
   residuals;
4. export all observables from the same branch;
5. compare without changing support, normalization, boundary data, or
   moment-shape coefficients afterward.

Source:

- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — headings `Frozen Inputs`, `Required Packet`, and `Non-Rescue Rules`

## What has been tested

The PDE audit reports that reduced frozen candidate sweeps and manufactured
nonlinear candidates fail the target packet. The separately frozen Stage-1
effective-closure branch reports $R_{\rm norm}\approx-10.8$, with outgoing
transfer roughly seven orders of magnitude below the target. That discrepancy
is a robust miss for the tested branch under its frozen assumptions.

The branch deliberately omits the complete matter/gauge-to-wall return
coupling, so the result does not decide full Path A and does not license
post-result fitting. The earlier Stage-1 README remains useful for solver
architecture but is stale as an outcome summary.

Sources:

- `research/pde_audit/notes/stage_v2_25_actual_branch_protocol_derivation.md` — heading `Status after the target-blind simulation miss`
- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — heading `Current status statement`
- `software/stage1_solver/STAGE1_VERDICT.md` — headings `The result — a robust MISS`, `Interpretation — what it means (and does NOT mean)`, and `Open after Stage-1 (deferred)`

## Remaining completion gates

The physical exporter still needs a self-consistent stationary background,
nonlinear material/throat equations, open-boundary and source/port laws,
stability certificates, return closure, and same-branch extraction of all
response and outgoing coefficients. Density, phase/current, sound speed,
localization, and feedback must be evaluated on that same material branch
rather than supplied as unrelated inputs.

Sources:

- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — heading `What the audit says is missing`
- `research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md` — headings `What is already present`, `What is missing`, and `Why this can block the final PDE`

Related pages:

- `memory/topics/foundations.md`
- `memory/topics/quadrupole-normalization.md`
- `memory/topics/weak-axisymmetric-and-orbit-closure.md`
- `memory/scripts/stage1-solver.md`
