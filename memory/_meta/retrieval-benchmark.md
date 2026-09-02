# Phase 0 Retrieval Benchmark

Status: baseline before population
Benchmark version: 2
Last updated: 2026-09-01

## Purpose and protocol

These are retrieval tests, not canned answers. Each item records the question,
the capabilities a passing memory must provide, and the original source paths
and stable anchors the answer should ultimately cite. It deliberately does not
record the expected prose conclusion or numerical result.

For the migration gate:

1. Start at `memory/index.md` in a fresh session.
2. Answer from `memory/` topic, source, conflict, and script pages. Do not read
   `atlas/` or `graph/`.
3. A normal answer should require the index plus no more than four additional
   memory pages. Opening an original source to verify a citation is permitted
   during citation audit, but should not be necessary to discover the answer.
4. Cite repository-relative original sources and the anchors below. A memory
   page or legacy node alone is not sufficient provenance.
5. State evidence, memory-review, lifecycle, scope, and unresolved boundaries. Do not
   turn script output into a conclusion the interpretive source does not make.
6. Record pass/fail and missing navigation or content. Do not edit the benchmark
   question to match an incomplete memory.

The proposed memory destinations below guide initial population. They may be
renamed once if `index.md` keeps an explicit redirect/crosswalk and the same
question remains answerable.

## Retired cases

The following IDs are permanently retired because their sources are outside
the selected corpus; IDs are not reused:

- `B03` — v3 work frontier.
- `B04` — v3 interface-response revision.
- `B06` — v3 S10 Python comparator.
- `B07` — v3 S10 Wolfram audit.

Active cases for version 2 are `B01`, `B02`, `B05`, and `B08`–`B12`.

## B01 — Latest moving-throat treatment

**Question:** Where is the current framework-level treatment of the
moving-throat PDE, what does that treatment claim to provide, and which part of
branch realization remains open?

**Expected memory destination:** `memory/topics/moving-throat-dynamics.md` plus the
capsule for the moving-throat PDE paper.

**Expected original source targets:**

- `research/pde/paper/pde.tex` — `\label{sec:status-parent}`
- `research/pde/paper/pde.tex` — `\label{sec:how-to-use}`
- `research/pde/paper/pde.tex` — `\label{sec:discussion-gap}`
- `research/pde/paper/pde.tex` — `\label{app:verification-map}`

**Pass conditions:** The answer identifies the current source, distinguishes
framework/reduction content from realization, surfaces the open boundary, and
does not present the paper or its scripts as a solved full nonlinear branch.

## B02 — Projection is not reduction

**Question:** Why must brane projection not be described as brane reduction,
and what invalid downstream inference does that distinction prevent?

**Expected memory destination:**
`memory/topics/projection-and-reduction.md`, including the migrated status
firewall re-anchored to original papers.

**Expected original source targets:**

- `research/pde/paper/pde.tex` — heading
  `\subsection{Projection versus reduction and the mixed-sector firewall}`
- `research/pde/paper/pde.tex` — `\label{eq:projected-continuity}` and
  `\label{eq:zero-mode-assumptions}`
- `research/4d/paper/4d.tex` — `\label{sec:projection}`
- `research/4d/paper/4d.tex` — `\label{sec:emreduction}`
- `research/4d_plasma/paper/4d_plasma.tex` —
  `\label{sec:projection_defs}`

**Pass conditions:** The answer distinguishes the operations by role and
assumptions, preserves the mixed/open-system caveat, and cites original TeX
rather than `FIREWALL_PROJECTION_NOT_REDUCTION` or another generated Atlas
page.

## B05 — Post-Newtonian chain and its status boundaries

**Question:** Map the paper chain from the 1PN bridge and full 1PN assembly
through 2PN, 2.5PN, 3PN, and 4PN. What role does each paper play, and where do
the sources themselves place exact, conditional, local/tail, or open
boundaries?

**Expected memory destination:** `memory/topics/post-newtonian-ladder.md` plus
the relevant paper capsules.

**Expected original source targets:**

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` —
  `\label{sec:exec-summary}` and `\label{sec:1pn-interaction-alpha-K}`
- `research/4d_1pn_full/paper/4d_1pn_full.tex` — `\label{sec:claims}`,
  `\label{sec:nonclaims}`, and `\label{sec:full-eih}`
- `research/4d_2pn/paper/4d_2pn.tex` — `\label{sec:intro-claims}`,
  `\label{sec:intro-nonclaims}`, and `\label{sec:full-assembly}`
- `research/4d_2_5pn/paper/4d_2_5pn.tex` —
  `\label{sec:conditional-theorem-gap}` and
  `\label{sec:conditional-theorem-gap-statement}`
- `research/4d_3pn/paper/4d_3pn.tex` — `\label{sec:final-main}` and
  `\label{sec:fixed-open}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:local-final}`,
  `\label{sec:tail-bridge}`, `\label{sec:final-conditional-4pn-theorem}`,
  and `\label{sec:fixed-open}`

**Pass conditions:** The answer provides a navigable lineage, does not flatten
different theorem scopes into one “PN proved” statement, and makes the 4PN
local/tail distinction and conditional/open boundaries retrievable.

## B08 — Stage-1 verdict

**Question:** What did Stage-1 test, what verdict did the frozen target-blind
run record, how was robustness/reproducibility bounded, and what does the
result explicitly not decide?

**Expected memory destination:**
`memory/topics/stage1-branch-realization.md`, the Stage-1 verdict capsule, and
`memory/scripts/stage1-solver.md`.

**Expected original source targets:**

- `software/stage1_solver/STAGE1_VERDICT.md` — headings `## The test`,
  `## The result — a robust MISS`, `## Review (all 4 gates clean)`,
  `## Interpretation — what it means (and does NOT mean)`, and
  `## Open after Stage-1 (deferred)`
- `software/stage1_solver/README.md` — headings `## How to run`,
  `## What to expect`, and `### What the numbers mean — and what they do **not**`
- `software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_report.md`
  — report title and frozen-packet identifiers

**Pass conditions:** The answer reports the scoped test and verdict with its
frozen provenance, preserves the branch/closure limitation and deferred work,
and does not generalize the result to the full framework.

## B09 — Electric sign and charge-selection boundary

**Question:** In the current puncture-deflection electric-sign result, which
objects and conditional calculations are decided, and which boundary/variant,
magnitude, and ontology selections remain unresolved?

**Expected memory destination:** `memory/topics/emergent-charge.md` and
`memory/scripts/em-charge-attribute.md`.

**Expected original source targets:**

- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` —
  heading `## Concise result (body; ≤2 pages)` and subsections `Q-FIELD / Q-AMEND`,
  `Q-BC`, `Q-FORCE / Q-COMBINE`, `Q-MAG / hooks`, and `§4 landing`
- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` —
  headings `### B. Full ensemble definitions`, `### C. Q-BC evidence`, and
  `### E. Production checks`
- `software/em_charge_attribute/reports/emergent_em_construction.md` — headings
  `## The ±w throat is an endpoint, not a renamed Z2 charge`,
  `## One Maxwell field and the dual sign`, and
  `## Honest scope and remaining obligations`
- `software/em_charge_attribute/puncture_deflection_electric_sign_check.py` —
  functions `classify_bc`, `section4_adjudicate`, `truth_table`, and `main`
- `software/em_charge_attribute/puncture_deflection_electric_sign_check.wl` —
  emitted tags `SECTION4_LANDING` and `ENGINE_AGREE`

**Pass conditions:** The answer keeps computed conditional results separate
from physical selection, preserves the sign/charge ontology firewall, and
does not turn a Tier-A or unresolved-boundary result into earned universal
charge.

## B10 — Moving-throat magnetism software result

**Question:** What do the two construction routes in the current moving-throat
magnetism result compute and compare, what software family reproduces them,
and which sign, normalization, cone, and active-flux questions remain open?

**Expected memory destination:** `memory/topics/emergent-magnetism.md` and
`memory/scripts/em-charge-attribute.md`.

**Expected original source targets:**

- `software/em_charge_attribute/magnetism_moving_throat_result.md` — heading
  `## Concise result (body; ≤2 pages)` and subsections `Q-CURRENT and field identity`,
  `Q-BOOST (Route A)`, `Q-DIRECT (Route B, blind to Route A)`,
  `Q-COMPARE / Q-MAG`, and `SEALED §4 landing`
- `software/em_charge_attribute/magnetism_moving_throat_result.md` — headings
  `### C. Exhaustive sealed truth table`, `### D. Hooks`, and
  `### E. Atomic production-path mutation campaigns`
- `software/em_charge_attribute/magnetism_moving_throat_check.py` — functions
  `derive_source`, `derive_route_b`, `derive_route_a`, `compare_routes`,
  `section4_adjudicate`, and `main`
- `software/em_charge_attribute/magnetism_moving_throat_check.wl` — emitted tags
  `SECTION4_LANDING`, `SECTION4_ALL_R1`, and `ENGINE_AGREE`

**Pass conditions:** The answer describes the constructed objects, route
independence and literal evidence chain without importing a desired result; it
separates decided tensor/diagnostic content from unresolved physical
normalization, selection, and non-conservative contributions.

## B11 — Material closure and nonlinear handoff

**Question:** Which ingredients of the moving-throat PDE are already present,
which constitutive/material closure is still missing, and why can that gap
block a final nonlinear PDE even when the framework and audits are internally
consistent?

**Expected memory destination:**
`memory/topics/material-closure-and-nonlinear-handoff.md` plus the PDE and
PDE-audit capsules.

**Expected original source targets:**

- `research/pde/paper/pde.tex` — `\label{sec:status-parent}` and
  `\label{sec:discussion-gap}`
- `research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md` —
  headings `## 1. What is already present`, `## 2. What is missing`, and
  `## 3. Why this can block the final PDE`
- `software/stage1_solver/STAGE1_VERDICT.md` — headings
  `## Interpretation — what it means (and does NOT mean)` and
  `## Open after Stage-1 (deferred)`

**Pass conditions:** The answer distinguishes kinematic/geometric structure
from material closure, preserves the scoped Stage-1 result, and does not treat
an audit or numerical miss as a completed nonlinear theory.

## B12 — Scoped Stage-1 status revision

**Question:** Which build/status expectation in the Stage-1 README is revised
by the later frozen verdict, what is the exact scope of that revision, and
which README architecture and run instructions remain current?

**Expected memory destination:** `memory/conflicts.md`,
`memory/topics/stage1-branch-realization.md`, and the Stage-1 source capsules.

**Expected original source targets:**

- `software/stage1_solver/README.md` — headings `## Build status`,
  `## How to run`, `## What to expect`, and
  `### What the numbers mean — and what they do **not**`
- `software/stage1_solver/STAGE1_VERDICT.md` — headings `## The test`,
  `## The result — a robust MISS`, and
  `## Interpretation — what it means (and does NOT mean)`

**Pass conditions:** The answer treats the verdict as a section-scoped
revision rather than whole-source supersession and keeps architecture,
commands, and explicit nonclaims available.

## Recording a run

Append run records only after the memory has been populated:

Record one row for every active version-2 ID: `B01`, `B02`, `B05`, and
`B08`–`B12`.

```markdown
### Run YYYY-MM-DD — <commit>

| Item | Pass | Memory pages opened | Citation audit | Gap |
|---|---|---:|---|---|
| B01 | yes/no | <count and paths> | pass/fail | <short gap or none> |
```

A benchmark run is not a physics review. It tests whether the memory preserves
scope, evidence, provenance, precedence, and retrieval boundaries already
present in the sources.
