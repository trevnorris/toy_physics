# Moving-Throat PDE — Compressed Program Rewrite Plan

## Purpose

This plan defines how to produce a fresh compressed program document after the
full derivation review, the SymPy/Mathematica hardening pass, the theorem-scope
tightening, and the physics red-team pass.

The target is **not** a polished paper.
The target is a working master document that is:

- compact enough for fast session handoff,
- complete enough to support new derivation work without re-reading all stages,
- explicit about claim status and theorem scope,
- and stable enough that later sessions do not have to reconstruct the logic from
  the full note tree every time.

---

## Why Rewrite Instead of Patch

The current handoff file
[moving_throat_pde_handoff_full.md](/var/projects/toy_physics/notes/moving_throat_pde_handoff_full.md)
is useful, but it is no longer the right shape for the present state of the project.

Current problems:

1. it is too large for fast working use;
2. it now mixes at least two distinct jobs:
   - PDE engine handoff
   - reduced variable/invariant dictionary
3. it contains duplicate top-level structure rather than one clean narrative;
4. it predates the latest theorem-scope tightening and physics red-team wording;
5. it is not organized around the actual remaining risk:
   reduced-closure structure vs full PDE realization.

So the next document should be written as a **fresh master handoff**, not as a
piecemeal patch to the old one.

---

## Primary Deliverable

Recommended target file:

- [moving_throat_pde_program_compact.md](/var/projects/toy_physics/notes/moving_throat_pde_program_compact.md)

Recommended strategy:

- leave [moving_throat_pde_handoff_full.md](/var/projects/toy_physics/notes/moving_throat_pde_handoff_full.md) in place as a legacy source;
- write the new compact document from scratch;
- only retire or relabel the old handoff after the new one is complete and checked.

Optional supporting deliverable:

- `notes/moving_throat/review/COMPRESSED_PROGRAM_SOURCE_MAP.md`

That source map would list, for each section of the new compact document, the
stage files and review artifacts it was distilled from.

---

## Document Contract

The new compact document should obey these rules.

## 1. One role only

It is a **working master program document**.
It is not:

- a stage-by-stage history,
- a paper draft,
- a speculative essay,
- or a giant appendix dump.

## 2. One definition per object

Each important object should be defined once in its canonical section:

- parent fields and action
- moving-throat geometry variables
- reduced core/mouth variables
- grouped `P2` variables
- coherent/invariant variables

Later sections may cite those definitions, but should not redefine them unless a
different regime is being introduced explicitly.

## 3. Explicit claim-status tags

Every major statement should be labeled in one of these categories:

- `Exact`
- `Exact Within Closure`
- `Numerically Located`
- `Reduced / Controlled Reduction`
- `Effective Closure`
- `Open`

The status system already used in the current handoff is good and should be kept,
but applied more consistently.

## 4. No duplicate theorem ledgers

Each theorem or branch-selection result should appear once in its final form.
Status/consolidation sections may summarize it, but the master statement should
not be repeated in multiple incompatible versions.

## 5. Open problems separated from solved structure

Every section should end cleanly with:

- what is fixed,
- what is fixed only within closure,
- what is numerically located,
- and what remains open.

That is the main safeguard against later overclaim drift.

---

## Inputs That Must Be Treated As Canonical

The rewrite should be sourced from these current artifacts, not just from the old handoff.

Core inputs:

- [moving_throat_pde_handoff_full.md](/var/projects/toy_physics/notes/moving_throat_pde_handoff_full.md)
- [ASSUMPTION_LEDGER.md](/var/projects/toy_physics/notes/moving_throat/review/ASSUMPTION_LEDGER.md)
- [DEPENDENCY_LEDGER.md](/var/projects/toy_physics/notes/moving_throat/review/DEPENDENCY_LEDGER.md)
- [THEOREM_TIGHTENING_NOTES.md](/var/projects/toy_physics/notes/moving_throat/review/THEOREM_TIGHTENING_NOTES.md)
- [PROOF_HARDENING_PLAN.md](/var/projects/toy_physics/notes/moving_throat/review/PROOF_HARDENING_PLAN.md)
- [PHYSICS_RED_TEAM_STATUS.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_RED_TEAM_STATUS.md)

Physics audit inputs:

- [PHYSICS_FOUNDATIONS_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_FOUNDATIONS_AUDIT.md)
- [PHYSICS_MOUTH_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_MOUTH_AUDIT.md)
- [PHYSICS_COEVOLUTION_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_COEVOLUTION_AUDIT.md)
- [PHYSICS_REALIZATION_INVARIANT_AUDIT.md](/var/projects/toy_physics/notes/moving_throat/review/PHYSICS_REALIZATION_INVARIANT_AUDIT.md)

Hardening result inputs:

- current stage notes after theorem-scope tightening
- saved numerical stress outputs where they support scope statements
- current SymPy/Mathematica audit state only as verification background, not as
  main narrative content

---

## Recommended Final Structure

The new compact document should have one continuous narrative with these parts.

1. Purpose and status firewall
2. Parent theory and exact bulk equations
3. Moving-throat geometry and throat-mode decomposition
4. Reduced wall/support/outgoing engine
5. Family-1 core/mouth closure stack
6. Co-evolving mouth/core branch
7. Grouped real `P2` response and outgoing normalization bridge
8. Coherent local-kernel reduction and microscopic invariant structure
9. Final theorem ledger
10. Open realization gap
11. Minimal source anchors

That gives one master document instead of the current split personality:

- engine narrative
- then dictionary narrative

The dictionary content should be embedded where it is used, not appended as a
second quasi-independent document.

---

## Chunked Rewrite Workflow

The rewrite should be done in discrete chunks.
Each chunk should be drafted, reviewed against the ledgers/audits, and only then
folded into the new master document.

## Chunk 0 — Scaffold and Section Contract

### Goal

Freeze the outline and the document rules before writing math.

### Inputs

- current handoff structure
- theorem tightening notes
- physics red-team status

### Output

- top-level headings in the new compact file
- section-level status conventions
- notation firewall block

### Done criteria

- every later chunk has a reserved destination section;
- no section is trying to do two jobs at once.

---

## Chunk 1 — Parent PDE Core

### Goal

Write the exact parent theory and exact bulk equations in their cleanest final form.

### Scope

- bulk arena, fields, indices
- corrected charge ontology
- exact parent action
- exact GNLS equation
- exact Maxwell equation
- exact continuity/current/gauge-invariant identities

### Primary sources

- existing handoff Sections `1-4`
- stages `001-004` where wording was tightened
- foundations physics audit

### Exclusions

- no reduced branch claims yet
- no throat-specific heuristic narrative beyond what is needed to define the geometry coupling

### Done criteria

- a fresh session can reconstruct the exact parent PDE setup from this chunk alone;
- no reduced closure statement is mislabeled as exact parent theory.

---

## Chunk 2 — Moving-Throat Geometry and Reduced Wall Engine

### Goal

Compress the throat definition, wall action, mode decomposition, and finite-throat
support branch into one coherent reduced-engine block.

### Scope

- level-set throat definition
- stationary reference throat
- moving-wall confinement coupling
- harmonic decomposition and grouped real `P2` role
- finite-throat D/N support problem
- mouth DtN operator
- minimal distributed wall action
- axisymmetric and grouped reductions

### Primary sources

- existing handoff Sections `5-9`
- stage `001-002` note clarifications

### Done criteria

- the reduced wall/support machinery is defined once and only once;
- the grouped `P2` bridge is motivated without duplicating later response formulas.

---

## Chunk 3 — Conservative and Outgoing Reduced Operator Stack

### Goal

State the reduced conservative kernels and outgoing bridge that the PDE must supply.

### Scope

- conservative BdG, Maxwell, and mixed moments
- outgoing transfer moments
- grouped-lane operator
- normalized grouped response and isotropy
- outgoing `l=2` fingerprint
- outgoing prefactor and quadrupole-normalization condition

### Primary sources

- existing handoff Sections `10-12`
- stages `003-004`, `089`, `153`
- theorem tightening notes for outgoing scope

### Done criteria

- the reduced target is explicit;
- the outgoing bridge is clearly marked as reduced/coherent where appropriate;
- the remaining theorem gap is not overstated.

---

## Chunk 4 — Family-1 Core/Mouth Closure Stack

### Goal

Compress the positive-source, mouth-law, and explicit Family-1 gain derivation into
one clean branch stack.

### Scope

- positive-source theorem
- mouth boundary layer
- core-to-mouth gain map
- normalized gain family
- actual Family-1 mouth gains
- self-matched susceptibility closure
- branch-selection and singular-limit results

### Primary sources

- stages `108`, `112`, `120-127`
- mouth physics audit

### Done criteria

- the branch-selection logic is explicit and one-lane/positive-source scope is visible;
- no natural/equal-normalized branch is accidentally presented as a finite regular branch;
- the self-matched closure is clearly labeled as same-layer effective closure.

---

## Chunk 5 — Co-Evolving Mouth/Core Branch

### Goal

Write the co-evolving closure, the renormalized canonical branch, and the real
self-consistency cost in one place.

### Scope

- co-evolving fixed-point map
- first-order defect transport
- frozen-traction fixed point
- numerically located renormalized canonical branch
- status of the canonical branch after backreaction

### Primary sources

- stages `137-140`
- co-evolution physics audit
- numerical stress results for `137-139`

### Done criteria

- the old canonical point and the renormalized canonical point are clearly separated;
- the renormalization cost is stated once with current numbers;
- uniqueness and monotonicity are tied to the analyzed positive branch window.

---

## Chunk 6 — Grouped `P2` Linear Response and Anisotropy Bridge

### Goal

Summarize the linear grouped anisotropy bridge without re-deriving every stage.

### Scope

- grouped `P2` outlet map
- direct coefficients `delta kappa_W`, `delta gamma_W`
- even-consistency relation
- collapse to the obstruction pair `(K_A, G_A)`

### Primary sources

- stages `142-153`
- realization/invariant physics audit

### Done criteria

- the grouped bridge is concise but mathematically sufficient;
- it is clearly labeled as a linear-response statement on the compensated isotropic branch.

---

## Chunk 7 — Coherent Local-Kernel and Invariant Structure

### Goal

Compress the coherent tracking branch, microscopic slippage variables, monomial
invariants, and quotient theorem into one final reduced-invariant block.

### Scope

- coherent local-kernel variables
- microscopic grouped drift variables
- observable defect variables
- branch composites
- three monomial invariants
- similarity orbit
- exact finite quotient theorem

### Primary sources

- stages `020-031`
- stages `165-170`
- theorem tightening notes
- realization/invariant physics audit

### Done criteria

- first-order vs finite statements are cleanly separated;
- the “orbit” language is explicitly reduced-structural, not fundamental gauge language;
- the final open realization question is stated exactly once.

---

## Chunk 8 — Final Theorem Ledger and Open Problems

### Goal

Close the master document with the shortest possible theorem/open-problem ledger.

### Scope

- what is exact
- what is exact within closure
- what is numerically located
- what remains open
- what the actual remaining theorem gap is

### Primary sources

- dependency ledger
- theorem tightening notes
- physics red-team status

### Done criteria

- a fresh session can tell immediately what is actually solved and what is not;
- the final open problem is stated as branch realization, not vague “more work.”

---

## Per-Chunk Review Checklist

Every chunk should be checked against the same list before it is considered done.

1. Are all symbols defined before use?
2. Is every nontrivial claim tagged with the right status?
3. Is the scope branch-local / closure-local / numerical where needed?
4. Is any theorem stronger than what the proof currently supports?
5. Does the chunk duplicate a formula that should live earlier?
6. Are the remaining open questions separated from solved structure?
7. Could a fresh session use this chunk without reading the full stage tree?

---

## Integration Strategy

Recommended writing order:

1. Chunk 0
2. Chunk 1
3. Chunk 2
4. Chunk 3
5. Chunk 4
6. Chunk 5
7. Chunk 6
8. Chunk 7
9. Chunk 8

Recommended session strategy:

- do **one chunk at a time**;
- after drafting each chunk, compare it against:
  - assumption ledger
  - dependency ledger
  - theorem tightening notes
  - relevant physics audit
- then integrate it into the new master document before moving on.

Do **not** try to rewrite all of it in one pass.
The risk of duplicated notation and scope drift rises too quickly if the whole
document is drafted at once.

---

## Completion Criteria

The new compact program document is done when all of the following are true.

1. A fresh session can recover the parent PDE, the reduced branch stack, the
   grouped `P2` bridge, and the final invariant theorem from one file.
2. The file contains no duplicate top-level narrative blocks.
3. The claim-status tagging matches the current hardening/red-team state.
4. The final open problem is stated as:
   the full PDE branch realization question.
5. The document is shorter, clearer, and more stable than
   [moving_throat_pde_handoff_full.md](/var/projects/toy_physics/notes/moving_throat_pde_handoff_full.md).

---

## Recommended Next Step

Start with `Chunk 0` and `Chunk 1` together.

That gives the new master document:

- its structure,
- its status firewall,
- and its exact parent PDE block.

Once those are stable, the later reduced-closure chunks can be added without
reopening the foundations every time.
