# Phase 0 Retrieval Benchmark

Status: baseline before population
Benchmark version: 1
Last updated: 2026-08-24

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

## B03 — Current v3 frontier and open gates

**Question:** What is the current v3 work frontier, which nearby step family is
closed, and which named requirements, defects, or gates still constrain the
next work?

**Expected memory destination:** `memory/topics/pde-ledger-v3-status.md` plus
the v3 governance and current-step capsules.

**Expected original source targets:**

- `research/pde_ledger_v3/CHARTER.md` — `{#two-halves}` and
  `{#falsification-standard}`
- `research/pde_ledger_v3/V3_STEP_PLAN.md` — `{#s11b}` and
  `{#s11b-split}`
- `research/pde_ledger_v3/steps/S11b_RUN_CHECKLIST.md` — heading
  `## ✅ Step record (step 12) DONE` and heading
  `## ✅ Card re-point (step 13) DONE`
- `research/pde_ledger_v3/steps/S11c_SCOPE.md` — headings `## The question`,
  `## The requirements`, and `## Scope boundaries (G14 — do not let C swallow neighbours)`
- `research/pde_ledger_v3/DEFECT_REGISTER.md` — heading
  `## C. Open negatives and gaps` and explicit anchors `{#c12}` through
  `{#c20}` where applicable
- `research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md` — heading `## Entries`

**Pass conditions:** The answer derives status from current governance and step
records rather than directory shape, separates done/in-progress/open, and does
not silently treat a dormant, named, or planned item as fixed.

## B04 — Explicit revision inside the v3 lineage

**Question:** Which later v3 step record explicitly revises or overturns an
earlier interface-response conclusion, what is the scope of that revision, and
what remains valid or open afterward?

**Expected memory destination:** `memory/conflicts.md` and
`memory/topics/pde-ledger-v3-status.md`.

**Expected original source targets:**

- `research/pde_ledger_v3/steps/S11bA_interface_response.md` — heading
  `## ⭐⭐⭐ The headline — the leak is FREQUENCY-DEPENDENT`
- `research/pde_ledger_v3/steps/S11bB_interface_assembly.md` — opening record
  boundary and heading
  `## ⭐⭐⭐ THE HEADLINE — the velocity-coupled leak costs an ENERGY RESERVOIR`
- `research/pde_ledger_v3/steps/S11bB_interface_assembly.md` — heading
  `### ⚠ Known limits — ⛔ recorded, not fixed`
- `research/pde_ledger_v3/steps/S11b_interface_coupling_law.md` — headings
  `## The result`, `## Static-or-instantaneous freezes, standing limits, and owed items`,
  and `## Provenance`

**Pass conditions:** The answer uses the explicit source relationship rather
than recency alone, limits supersession to the affected statement/scope,
retains the older position as lineage, and preserves the later record's known
limits.

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

## B06 — Representative Python comparator

**Question:** What does
`research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py` consume and
compare, what categories and artifacts does it emit, and what does its final
guard not establish?

**Expected memory destination:** `memory/scripts/pde-ledger-v3.md`.

**Expected original source targets:**

- `research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py` — module
  docstring, class `ComparatorInputError`, and functions `read_transcript`,
  `compare_records`, `compare_nullspace_basis`, `run`, and `main`
- `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out` — tags
  `CATEGORY: SUMMARY` and `FINAL_GUARD`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading
  `## Live comparator: measured scope, not a blanket verdict`

**Pass conditions:** The catalog entry names the input transcript grammar,
comparison/output categories, tracked output, and guard boundary. It attributes
interpretation to the step record and does not paraphrase an exit code as a
physical verdict.

## B07 — Representative Wolfram audit and paper-facing provenance

**Question:** What object family does
`research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
construct and emit, and how can an agent trace a paper-facing S10 statement
from card to step record, engine transcript, and comparator?

**Expected memory destination:** `memory/scripts/pde-ledger-v3.md` plus the S10
source-unit capsule.

**Expected original source targets:**

- `research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl`
  — functions `buildPackage`, `runPackage`, `runSpectrum`, `runModeSet`, and
  emitted tags `WL_S10_RUN_PAIRS`, `WL_S10_SKIPPED_PAIRS`, and
  `WL_S10_EMITTED_TAG_COUNT`
- `research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out`
  — corresponding `WL_S10_*` tags
- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` —
  `\label{step:S10}`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — headings
  `## What was supplied and what was computed`,
  `## Measured MAIN spectrum and count`, and
  `## Claims this step still does not establish`
- `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out` —
  `CATEGORY: SUMMARY`

**Pass conditions:** The answer distinguishes supplied premises from computed
objects, names literal output provenance, traces the interpretive hop, and
retains the step's explicit nonclaims and comparator-coverage limits.

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

## Recording a run

Append run records only after the memory has been populated:

```markdown
### Run YYYY-MM-DD — <commit>

| Item | Pass | Memory pages opened | Citation audit | Gap |
|---|---|---:|---|---|
| B01 | yes/no | <count and paths> | pass/fail | <short gap or none> |
```

A benchmark run is not a physics review. It tests whether the memory preserves
scope, evidence, provenance, precedence, and retrieval boundaries already
present in the sources.
